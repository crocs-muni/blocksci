#include "coinjoin_module_py.hpp"

#include <blocksci/address/address.hpp>
#include <blocksci/chain/access.hpp>
#include <blocksci/chain/blockchain.hpp>
#include <blocksci/cluster/cluster.hpp>
#include <blocksci/heuristics/tx_identification.hpp>
#include <blocksci/scripts/script_range.hpp>
#include <queue>
#include <unordered_map>

#include "../external/json/single_include/nlohmann/json.hpp"
#include "caster_py.hpp"
#include "sequence.hpp"

struct CoinjoinNamespace {};

namespace py = pybind11;
using namespace blocksci;

using json = nlohmann::json;

void init_coinjoin_module(py::class_<Blockchain> &cl) {
    cl.def(
          "find_friends_who_dont_pay",
          [](Blockchain &chain, BlockHeight start, BlockHeight stop) {
              using MapType = std::vector<Transaction>;
              auto reduce_func = [](MapType &vec1, MapType &vec2) -> MapType & {
                  vec1.reserve(vec1.size() + vec2.size());
                  vec1.insert(vec1.end(), std::make_move_iterator(vec2.begin()), std::make_move_iterator(vec2.end()));
                  return vec1;
              };

              auto map_func = [](const Transaction &tx) -> MapType {
                  MapType result;
                  // if it is a ww2 coinjoin, then
                  if (blocksci::heuristics::isWasabi2CoinJoin(tx)) {
                      return {};
                  }

                  // go through all the outputs
                  for (const auto &output : tx.outputs()) {
                      if (!output.isSpent()) continue;
                      // ignore direct coinjoin remixes
                      if (blocksci::heuristics::isWasabi2CoinJoin(output.getSpendingTx().value())) {
                          continue;
                      }

                      // find an output that has all inputs from a coinjoin
                      bool all_inputs_from_cj = true;
                      for (auto input : output.getSpendingTx().value().inputs()) {
                          if (!blocksci::heuristics::isWasabi2CoinJoin(input.getSpentTx())) {
                              all_inputs_from_cj = false;
                              break;
                          }
                      }
                      if (!all_inputs_from_cj) {
                          continue;
                      }

                      // and check whether at least one of the outputs is mixed again
                      for (auto output2 : output.getSpendingTx().value().outputs()) {
                          if (!output2.isSpent()) continue;
                          if (blocksci::heuristics::isWasabi2CoinJoin(output2.getSpendingTx().value())) {
                              result.push_back(output.getSpendingTx().value());
                              break;
                          }
                      }
                  }
                  return result;
              };

              return chain[{start, stop}].mapReduce<MapType, decltype(map_func), decltype(reduce_func)>(map_func,
                                                                                                        reduce_func);
          },
          "Filter the blockchain to only include 'friends don't pay' transactions.", pybind11::arg("start"),
          pybind11::arg("stop"))
        .def(
            "filter_coinjoin_txes",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, std::string coinjoinType) {
                auto txes = chain[{start, stop}].filter([&](const Transaction &tx) {
                    return blocksci::heuristics::isCoinjoinOfGivenType(tx, coinjoinType);
                });

                // filter txes which are not connected to any other coinjoin tx
                std::unordered_set<Transaction> txSet;
                for (const auto &tx : txes) {
                    txSet.insert(tx);
                }

                std::vector<Transaction> result;
                result.push_back(txes[0]);

                for (const auto &tx : txSet) {
                    for (const auto &input : tx.inputs()) {
                        if (txSet.find(input.getSpentTx()) != txSet.end()) {
                            result.push_back(tx);
                            break;
                        }
                    }
                }
                return result;
            },
            "Filter coinjoin transactions", pybind11::arg("start"), pybind11::arg("stop"),
            pybind11::arg("coinjoin_type"))
        .def(
            "find_hw_sw_coinjoins",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop) {
                return chain[{start, stop}].filter([](const Transaction &tx) {
                    auto result = blocksci::heuristics::isLongDormantInRemixes(tx);

                    if (result == blocksci::heuristics::HWWalletRemixResult::False) {
                        return false;
                    }

                    if (result == blocksci::heuristics::HWWalletRemixResult::Trezor) {
                        return true;
                    }

                    return true;
                });
            },
            "Filter hw_sw coinjoin transactions", pybind11::arg("start"), pybind11::arg("stop"))

        .def(
            "find_traverses_between_coinjoins",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop) {
                using MapType = std::unordered_map<int,
                                                   std::unordered_map<uint64_t,  // <from coinjoin type> in uint64_t
                                                                      std::unordered_map<uint64_t,  // <to coinjoin
                                                                                                    // type> in uint64_t
                                                                                         uint64_t>>>;
                auto mapFunc = [](const Transaction &tx) -> MapType {
                    auto tag = blocksci::heuristics::getCoinjoinTag(tx);
                    if (tag == blocksci::heuristics::CoinJoinType::None) {
                        return {};
                    }

                    auto isDifferentCoinJoinThanTag = [&](blocksci::heuristics::CoinJoinType tag,
                                                          blocksci::heuristics::CoinJoinType newTag) {
                        return newTag != tag;
                    };
                    MapType result;

                    if (result.find(tx.block().height()) == result.end()) {
                        result[tx.block().height()] = {};
                    }

                    auto handleWhirlpoolTx0 = [&](const Transaction &tx, blocksci::heuristics::CoinJoinType oldTag) {
                        for (const auto &output : tx.outputs()) {
                            if (!output.isSpent()) continue;
                            auto spendingTx = output.getSpendingTx().value();
                            auto newTag = blocksci::heuristics::getCoinjoinTag(spendingTx);
                            if (newTag == blocksci::heuristics::CoinJoinType::None) {
                                continue;
                            }
                            if (result[spendingTx.block().height()].find(static_cast<uint64_t>(oldTag)) ==
                                result[spendingTx.block().height()].end()) {
                                result[spendingTx.block().height()][static_cast<uint64_t>(oldTag)] = {};
                            }

                            if (isDifferentCoinJoinThanTag(oldTag, newTag)) {
                                result[spendingTx.block().height()][static_cast<uint64_t>(oldTag)]
                                      [static_cast<uint64_t>(newTag)]++;
                            }
                        }
                    };

                    for (const auto &output : tx.outputs()) {
                        if (!output.isSpent()) continue;
                        auto spendingTx = output.getSpendingTx().value();
                        auto newTag = blocksci::heuristics::getCoinjoinTag(spendingTx);
                        if (newTag == blocksci::heuristics::CoinJoinType::None) {
                            handleWhirlpoolTx0(spendingTx, tag);
                            continue;
                        }

                        if (result[spendingTx.block().height()].find(static_cast<uint64_t>(tag)) ==
                            result[spendingTx.block().height()].end()) {
                            result[spendingTx.block().height()][static_cast<uint64_t>(tag)] = {};
                        }

                        if (isDifferentCoinJoinThanTag(tag, newTag)) {
                            result[spendingTx.block().height()][static_cast<uint64_t>(tag)]
                                  [static_cast<uint64_t>(newTag)]++;
                        }
                    }
                    return result;
                };

                auto reduceFunc = [](MapType &our, MapType &their) -> MapType & {
                    for (auto &[blockHeight, traverses] : their) {
                        if (our.find(blockHeight) == our.end()) {
                            our[blockHeight] = {};
                        }
                        auto &ourFrom = our[blockHeight];
                        for (auto &[fromCj, toCjs] : traverses) {
                            if (ourFrom.find(fromCj) == ourFrom.end()) {
                                ourFrom[fromCj] = {};
                            }
                            auto &ourTo = ourFrom[fromCj];
                            for (auto &[toCj, howMany] : toCjs) {
                                if (ourTo.find(toCj) == ourTo.end()) {
                                    ourTo[toCj] = 0;
                                }
                                ourTo[toCj] += howMany;
                            }
                        }
                    }
                    return our;
                };

                try {
                    return chain[{start, stop}].mapReduce<MapType, decltype(mapFunc), decltype(reduceFunc)>(mapFunc,
                                                                                                            reduceFunc);
                } catch (const std::exception &e) {
                    throw std::runtime_error(std::string("Error in mapReduce: ") + e.what());
                }
            },
            "Filter the blockchain to only include transactions that traverse between two coinjoins.",
            pybind11::arg("start"), pybind11::arg("stop"))

        .def(
            "get_coinjoin_consolidations",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, double inputOutputRatio,
               std::string coinjoinType, int maxHops) {
                using ResultType =
                    std::map<std::string, std::vector<Transaction>>;              // consolidation_type, [input_tx_hash]
                using MapType = std::vector<std::pair<Transaction, ResultType>>;  // tx_hash, ResultType

                auto globalVisited = std::unordered_set<uint256>();
                auto reduceFunc = [](MapType &vec1, MapType &vec2) -> MapType & {
                    vec1.reserve(vec1.size() + vec2.size());
                    vec1.insert(vec1.end(), std::make_move_iterator(vec2.begin()), std::make_move_iterator(vec2.end()));
                    return vec1;
                };

                auto mapFunc = [&](const Transaction &tx) -> MapType {
                    if (!blocksci::heuristics::isCoinjoinOfGivenType(tx, coinjoinType)) {
                        return {};
                    }

                    ResultType result;
                    result["certain"] = {};
                    result["possible"] = {};

                    std::queue<std::pair<const Transaction &, int>> bfsQueue;
                    std::unordered_set<uint256> visited;

                    bfsQueue.push({tx, 0});
                    visited.insert(tx.getHash());

                    while (!bfsQueue.empty()) {
                        auto [currentTx, depth] = bfsQueue.front();
                        bfsQueue.pop();

                        if (depth > maxHops) continue;

                        for (const auto &output : currentTx.outputs()) {
                            if (!output.isSpent()) continue;
                            auto spendingTx = output.getSpendingTx().value();

                            if (visited.count(spendingTx.getHash())) continue;
                            visited.insert(spendingTx.getHash());

                            if (blocksci::heuristics::isCoinjoinOfGivenType(spendingTx, coinjoinType)) {
                                // bfs_queue.push({spending_tx, depth + 1});
                                continue;
                            }

                            // if depth is 0, then check, if all the inputs are from the coinjoin
                            if (depth == 0) {
                                bool allInputsFromCj = true;
                                for (auto input : spendingTx.inputs()) {
                                    if (!blocksci::heuristics::isCoinjoinOfGivenType(input.getSpentTx(), coinjoinType)) {
                                        allInputsFromCj = false;
                                        break;
                                    }
                                }
                                if (!allInputsFromCj) {
                                    continue;
                                }
                            }

                            // Check if the spending tx is a consolidation tx
                            auto consolidationType =
                                blocksci::heuristics::getConsolidationType(spendingTx, inputOutputRatio);
                            if (consolidationType == blocksci::heuristics::ConsolidationType::Certain) {
                                result["certain"].push_back(spendingTx);
                                continue;
                            } else if (consolidationType == blocksci::heuristics::ConsolidationType::Possible) {
                                result["possible"].push_back(spendingTx);
                                continue;
                            }

                            // If it's not a consolidation, continue BFS
                            bfsQueue.push({spendingTx, depth + 1});
                        }
                    }

                    // sort "certain" and "possible" txes by total output value
                    auto sortingFn = [](const Transaction &tx1, const Transaction &tx2) {
                        auto sumFn = [](int64_t sum, const Output &output) { return sum + output.getValue(); };
                        return std::accumulate(tx1.outputs().begin(), tx1.outputs().end(), 0, sumFn) >
                               std::accumulate(tx2.outputs().begin(), tx2.outputs().end(), 0, sumFn);
                    };
                    std::sort(result["certain"].begin(), result["certain"].end(), sortingFn);

                    return {{tx, result}};
                };

                return chain[{start, stop}].mapReduce<MapType, decltype(mapFunc), decltype(reduceFunc)>(mapFunc, reduceFunc);
            },
            "Filter certain consolidation transactions", pybind11::arg("start"), pybind11::arg("stop"),
            pybind11::arg("input_output_ratio"), pybind11::arg("coinjoin_type"), pybind11::arg("hops"))
        .def(
            "compute_anonymity_degradation",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, int daysToConsider, std::string coinjoinType, std::optional<std::string> coinjoinSubType) {
                // CJTX and its anonymity sets
                using AnonymitySetsFuncType = std::unordered_map<Transaction, std::unordered_map<int64_t, int64_t>>;
                // For txes which consolidate inputs from multiple cjtxes. CJTX = Transaction, map = <value, count>
                // value = the anonymity set, count = how many times it appears in the consolidation tx (how much it
                // degrades the anonymity set)
                using PointingToTransactionsType =
                    std::unordered_map<Transaction, std::unordered_map<int64_t, int64_t>>;

                auto filteringFunc = [&](const Transaction &tx) -> bool {
                    return blocksci::heuristics::isCoinjoinOfGivenType(tx, coinjoinType, coinjoinSubType);
                };

                auto mapFunc = [&](const Transaction &tx) -> AnonymitySetsFuncType {
                    if (!filteringFunc(tx)) {
                        return {};
                    }

                    AnonymitySetsFuncType result;
                    auto anonymitySets = std::unordered_map<int64_t, int64_t>();

                    for (const auto &output : tx.outputs()) {
                        anonymitySets[output.getValue()]++;
                    }

                    result[tx] = anonymitySets;
                    return result;
                };

                auto reduceFunc = [](AnonymitySetsFuncType &map1,
                                     AnonymitySetsFuncType &map2) -> AnonymitySetsFuncType & {
                    for (const auto &[key, value] : map2) {
                        if (map1.find(key) == map1.end()) {
                            map1[key] = value;
                        } else {
                            for (const auto &[key2, value2] : value) {
                                if (map1[key].find(key2) == map1[key].end()) {
                                    map1[key][key2] = value2;
                                } else {
                                    map1[key][key2] += value2;
                                }
                            }
                        }
                    }
                    return map1;
                };

                auto initialAnonymitySets =
                    chain[{start, stop}].mapReduce<AnonymitySetsFuncType, decltype(mapFunc), decltype(reduceFunc)>(
                        mapFunc, reduceFunc);

                auto mapFunc2 = [&](const Transaction &tx) -> PointingToTransactionsType {
                    PointingToTransactionsType result;
                    auto coinJoinTag = blocksci::heuristics::getCoinjoinTag(tx);
                    if (coinJoinTag != blocksci::heuristics::CoinJoinType::None) {
                        return {};
                    }
                    for (const auto &input : tx.inputs()) {
                        auto inputTx = input.getSpentTx();
                        if (tx.block().timestamp() - inputTx.block().timestamp() > daysToConsider * 24 * 60 * 60) {
                            continue;
                        }

                        if (initialAnonymitySets.find(inputTx) == initialAnonymitySets.end()) {
                            continue;
                        }

                        if (result.find(inputTx) == result.end()) {
                            result[inputTx] = {};
                        }

                        if (result[inputTx].find(input.getValue()) == result[inputTx].end()) {
                            result[inputTx][input.getValue()] = 1;
                        } else {
                            result[inputTx][input.getValue()]++;
                        }
                    }
                    auto sum = 0;

                    for (const auto &[cjtx, anonymitySet] : result) {
                        for (const auto &[value, count] : anonymitySet) {
                            sum += count;
                        }
                    }

                    return sum > 1 ? result : PointingToTransactionsType{};
                };

                auto reduceFunc2 = [&](PointingToTransactionsType &map1,
                                       PointingToTransactionsType &map2) -> PointingToTransactionsType & {
                    for (const auto &[tx, anonymitySets] : map2) {
                        if (map1.find(tx) == map1.end()) {
                            map1[tx] = anonymitySets;
                        } else {
                            for (const auto &[value, count] : anonymitySets) {
                                if (map1[tx].find(value) == map1[tx].end()) {
                                    map1[tx][value] = count;
                                } else {
                                    map1[tx][value] += count;
                                }
                            }
                        }
                    }
                    return map1;
                };
                if (daysToConsider > 0) {
                    auto pointingToTransactions =
                        chain[{start, stop}]
                            .mapReduce<PointingToTransactionsType, decltype(mapFunc2), decltype(reduceFunc2)>(
                                mapFunc2, reduceFunc2);

                    for (const auto &[tx, anonymitySets] : pointingToTransactions) {
                        for (const auto &[value, count] : anonymitySets) {
                            initialAnonymitySets[tx][value] -= count;
                        }
                    }
                }

                std::unordered_map<Transaction, std::pair<double, int64_t>> result;
                for (const auto &[tx, anonymitySets] : initialAnonymitySets) {
                    double resultValue = 0;
                    int64_t totalCount = 0;
                    if (coinjoinType == "wasabi1") {
                        // take only standard denominations
                    }
                    for (const auto &[key, value] : anonymitySets) {
                        resultValue += lgamma(value + 1) / log(2);
                        totalCount += value;
                    }
                    result[tx] = std::make_pair(resultValue, totalCount);
                }

                return result;
            },
            "Compute anonymity degradation in coinjoins", pybind11::arg("start"), pybind11::arg("stop"),
            pybind11::arg("daysToConsider"), pybind11::arg("coinjoinType"), pybind11::arg("coinjoinSubType") = py::none());
    ;
}