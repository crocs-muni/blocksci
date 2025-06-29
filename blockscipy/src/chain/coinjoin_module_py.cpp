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
#include <optional>

struct CoinjoinNamespace {};

namespace py = pybind11;
using namespace blocksci;

using json = nlohmann::json;

std::unordered_set<Transaction> findLinkedCjTxes(
    int start, 
    int stop, 
    std::string coinjoinType,                                            
    Blockchain &chain, 
    std::optional<std::string> subtype = std::nullopt,
    std::optional<std::unordered_set<std::string>> falsePositives = std::nullopt
) {
        auto txes = chain[{start, stop}].filter([&](const Transaction &tx) {
            return blocksci::heuristics::isCoinjoinOfGivenType(tx, coinjoinType, subtype);
        });

        if (txes.empty()) {
            return {};
        }

        // filter txes which are not connected to any other coinjoin tx
        std::unordered_set<Transaction> txSet;
        for (const auto &tx : txes) {
            txSet.insert(tx);
        }

        std::unordered_set<Transaction> result;
        result.insert(txes[0]);

        for (const auto &tx : txSet) {
            if (falsePositives.has_value() && falsePositives.value().find(tx.getHash().GetHex()) != falsePositives.value().end()) {
                continue;
            }
            for (const auto &input : tx.inputs()) {
                if (txSet.find(input.getSpentTx()) != txSet.end()) {
                    result.insert(tx);
                    break;
                }
            }
        }
        return result;
}

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
            "find_kruw_paths_to_bybit",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, pybind11::set kruwTxIds, pybind11::set bybitAddresses,
               int maxHops) {
                std::unordered_set<Address> bybitAddressesSet = pybind11::cast<std::unordered_set<Address>>(bybitAddresses);
                std::unordered_set<Transaction> kruwTxSet = pybind11::cast<std::unordered_set<Transaction>>(kruwTxIds);
                std::cout << "KRUW txes: " << kruwTxSet.size() << std::endl;
                std::cout << "Bybit addresses: " << bybitAddressesSet.size() << std::endl;
                std::cout << "Max hops: " << maxHops << std::endl;
                
                
                // MapType = <tx_hash, <address, <<path>>>>
                using TxHashType = std::string;
                using TxPathType = std::vector<TxHashType>;
                using MapType = std::unordered_map<Transaction, std::unordered_map<Address, std::vector<TxPathType>>>;
                auto reduceFunc = [](MapType &our, MapType &their) -> MapType & {
                    for (auto &[tx, addresses] : their) {
                        if (our.find(tx) == our.end()) {
                            our[tx] = {};
                        }
                        auto &ourAddresses = our[tx];
                        for (auto &[address, paths] : addresses) {
                            if (ourAddresses.find(address) == ourAddresses.end()) {
                                ourAddresses[address] = {};
                            }

                            for (auto &path : paths) {
                                auto reversedPath = path;
                                std::reverse(reversedPath.begin(), reversedPath.end());
                                ourAddresses[address].push_back(reversedPath);
                            }

                        }
                    }
                    return our;
                };
                

                auto mapFunc = [&](const Transaction &tx) -> MapType {
                    if (kruwTxSet.find(tx) == kruwTxSet.end()) {
                        return {};
                    }
                    MapType result;

                    // do a bfs and check if any of the inputs come from the bybit addresses within maxHops
                    std::queue<std::pair<Transaction, TxPathType>> bfsQueue;
                    std::unordered_set<Transaction> visited;
                    bfsQueue.push({tx, {tx.getHash().GetHex()}});
                    while (!bfsQueue.empty()) {
                        auto [currentTx, path] = bfsQueue.front();
                        bfsQueue.pop();

                        if (path.size() > maxHops) {
                            continue;
                        }
                        for (const auto &input : currentTx.inputs()) {
                            // remove non-fresh inputs (coming from a coinjoin)
                            if (blocksci::heuristics::isWasabi2CoinJoin(input.getSpentTx())) {
                                continue;
                            }
                            auto a = input.getAddress();

                            if (bybitAddressesSet.find(a) != bybitAddressesSet.end()) {
                                if (result.find(tx) == result.end()) {
                                    result[tx] = {};
                                }
                                if (result[tx].find(a) == result[tx].end()) {
                                    result[tx][a] = {};
                                }
                                result[tx][a].push_back(path);
                            }

                            auto inputTx = input.getSpentTx();
                            if (visited.find(inputTx) != visited.end()) {
                                continue;
                            }

                            auto newPath = path;
                            newPath.push_back(inputTx.getHash().GetHex());
                            bfsQueue.push({inputTx, newPath});

                            visited.insert(inputTx);
                        }
                    }
                    
                    
                    return result;
                };

                return chain[{start, stop}].mapReduce<MapType, decltype(mapFunc), decltype(reduceFunc)>(mapFunc,
                                                                                                        reduceFunc);
            },
            pybind11::arg("start"), pybind11::arg("stop"),
            pybind11::arg("kruw_tx_ids"), pybind11::arg("bybit_addresses"), pybind11::arg("max_hops")
        )
        .def(
            "filter_coinjoin_txes",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, std::string coinjoinType) {
                return findLinkedCjTxes(start, stop, coinjoinType, chain);
                
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
            [](
                Blockchain &chain, 
                BlockHeight start, 
                BlockHeight stop, 
                int daysToConsider, 
                std::string coinjoinType, 
                std::optional<std::string> coinjoinSubType, 
                bool ignoreNonStandardDenominations,
                bool ignoreRemixes,
                std::optional<std::tuple<int64_t, int64_t>> ww2DenomsBucket,
		std::optional<std::unordered_set<std::string>> falseCoinjoins
            ) {
                std::unordered_map<std::string, std::unordered_set<Transaction>> cjsOfGivenType;
                cjsOfGivenType["wasabi1"] = findLinkedCjTxes(start, stop, "wasabi1", chain, std::nullopt, falseCoinjoins);
	    	cjsOfGivenType["wasabi2"] = findLinkedCjTxes(start, stop, "wasabi2", chain, std::nullopt, falseCoinjoins);
                cjsOfGivenType["whirlpool"] = findLinkedCjTxes(start, stop, "whirlpool", chain, std::nullopt, falseCoinjoins);
                if (coinjoinSubType.has_value()) {
                    cjsOfGivenType[coinjoinSubType.value()] = findLinkedCjTxes(start, stop, coinjoinType, chain, coinjoinSubType.value(), falseCoinjoins);
                }


                // CJTX and its anonymity sets
                using AnonymitySetsFuncType = std::unordered_map<Transaction, std::unordered_map<int64_t, int64_t>>;
                // For txes which consolidate inputs from multiple cjtxes. CJTX = Transaction, map = <value, count>
                // value = the anonymity set, count = how many times it appears in the consolidation tx (how much it
                // degrades the anonymity set)
                using PointingToTransactionsType =
                    std::unordered_map<Transaction, std::unordered_map<int64_t, int64_t>>;

                auto shouldIgnoreDenomination = [&](const int64_t &outputValue) -> bool {
                    if (ww2DenomsBucket.has_value()) {
                        auto [minDenom, maxDenom] = ww2DenomsBucket.value();
                        if (outputValue < minDenom || outputValue >= maxDenom) {
                            return true;
                        }
                    }
                    if (!ignoreNonStandardDenominations) {
                        return false;
                    }
                    if (coinjoinType == "whirlpool") {
                        return false;
                    }
                    // Check if the output is a standard denomination
                    if (coinjoinType == "wasabi1") {
                        return outputValue > 0.12 * 1e8 || outputValue < 0.08 * 1e8;
                    }
                    if (coinjoinType == "wasabi2") {
                        return CoinjoinUtils::ww2_denominations.find(outputValue) == CoinjoinUtils::ww2_denominations.end();
                    }

                    return false;
                    
                };

                auto getSizesOfAnonymitySetsPerOutput = [&](const Transaction &tx) -> AnonymitySetsFuncType {
                    if (coinjoinSubType.has_value()) {
                        if (cjsOfGivenType[coinjoinSubType.value()].find(tx) == cjsOfGivenType[coinjoinSubType.value()].end()) {
                            return {};
                        }
                    } else {
                        if (cjsOfGivenType[coinjoinType].find(tx) == cjsOfGivenType[coinjoinType].end()) {
                            return {};
                        }
                    }

                    AnonymitySetsFuncType result;
                    auto anonymitySets = std::unordered_map<int64_t, int64_t>();

                    auto sum = 0;
                    for (const auto &output : tx.outputs()) {
                        // if it is not yet spent we don't really care about it
                        // if it is spent, we need to check if it is a coinjoin
                        if (output.isSpent()) {
                            if (!output.getSpendingTx().has_value()) {
                                continue;
                            }
                            auto spendingTx = output.getSpendingTx().value();
                            if (ignoreRemixes && !(cjsOfGivenType["wasabi1"].find(spendingTx) == cjsOfGivenType["wasabi1"].end() &&
                                cjsOfGivenType["wasabi2"].find(spendingTx) == cjsOfGivenType["wasabi2"].end() &&
                                cjsOfGivenType["whirlpool"].find(spendingTx) == cjsOfGivenType["whirlpool"].end())) {
                                continue;
                            }

                        }                         
                        if (shouldIgnoreDenomination(output.getValue())) {
                            continue;
                        }
                        

                        if (anonymitySets.find(output.getValue()) == anonymitySets.end()) {
                            anonymitySets[output.getValue()] = 0;
                        } 
                        anonymitySets[output.getValue()]++;
                        sum++;
                    }

                    result[tx] = anonymitySets;

                    return result;
                };

                auto combineAnonymitySetSizes = [](AnonymitySetsFuncType &map1,
                                     AnonymitySetsFuncType &map2) -> AnonymitySetsFuncType & {
                    for (const auto &[key, value] : map2) {
                        if (map1.find(key) == map1.end()) {
                            map1[key] = value;
                            continue;
                        } 
                        for (const auto &[key2, value2] : value) {
                            if (map1[key].find(key2) == map1[key].end()) {
                                map1[key][key2] = value2;
                            } else {
                                map1[key][key2] += value2;
                            }
                        }
                    }
                    return map1;
                };

                auto initialAnonymitySets =
                    chain[{start, stop}].mapReduce<AnonymitySetsFuncType, decltype(getSizesOfAnonymitySetsPerOutput), decltype(combineAnonymitySetSizes)>(
                        getSizesOfAnonymitySetsPerOutput, combineAnonymitySetSizes);

                auto decreaseAnonymityIfConsolidated = [&](const Transaction &tx) -> PointingToTransactionsType {
                    PointingToTransactionsType result;
                    // if it _is_ coinjoin and we ignore remixes
                    if (!(cjsOfGivenType["wasabi1"].find(tx) == cjsOfGivenType["wasabi1"].end() &&
                        cjsOfGivenType["wasabi2"].find(tx) == cjsOfGivenType["wasabi2"].end() &&
                        cjsOfGivenType["whirlpool"].find(tx) == cjsOfGivenType["whirlpool"].end())) {
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

                        if (shouldIgnoreDenomination(input.getValue())) {
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

                auto combineDataAfterAnonymityLoss = [&](PointingToTransactionsType &map1,
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
                            .mapReduce<PointingToTransactionsType, decltype(decreaseAnonymityIfConsolidated), decltype(combineDataAfterAnonymityLoss)>(
                                decreaseAnonymityIfConsolidated, combineDataAfterAnonymityLoss);

                    for (const auto &[tx, anonymitySets] : pointingToTransactions) {
                        for (const auto &[value, count] : anonymitySets) {
                            initialAnonymitySets[tx][value] -= count;
                        }
                    }
                }

                std::unordered_map<Transaction, std::unordered_map<std::string, int64_t>> result;
                for (const auto &[tx, anonymitySets] : initialAnonymitySets) {
                    int64_t notIgnoredOutputs = 0;
                    int64_t totalCount = 0;
                    int64_t notRemixedOutputs = 0;
                    int64_t notRemixedNotIgnoredOutputs = 0;
                    for (const auto &output : tx.outputs()) {
                        int raised = 0;
                        if (!shouldIgnoreDenomination(output.getValue())) {
                            notIgnoredOutputs++;
                            raised++;
                        }
                        if (!output.isSpent()) {
                            notRemixedOutputs++;
                            raised++;
                        } else {
                            auto spendingTx = output.getSpendingTx().value();
                            if (cjsOfGivenType["wasabi1"].find(spendingTx) == cjsOfGivenType["wasabi1"].end() &&
                                cjsOfGivenType["wasabi2"].find(spendingTx) == cjsOfGivenType["wasabi2"].end() &&
                                cjsOfGivenType["whirlpool"].find(spendingTx) == cjsOfGivenType["whirlpool"].end()) {
                                notRemixedOutputs++;
                                raised++;
                            }
                        }
                        if (raised == 2) {
                            notRemixedNotIgnoredOutputs++;
                        }
                    }
                    
                    for (const auto &[key, value] : anonymitySets) {
                        // resultValue += lgamma(value + 1) / log(2);
                        totalCount += value;
                    }
                    // result[tx] = std::make_pair(resultValue, totalCount);
                    result[tx] = {
                        {"total_count", totalCount},
                        {"all_outputs", tx.outputCount()},
                        {"not_remixed_outputs", notRemixedOutputs},
                        {"not_ignored_outputs", notIgnoredOutputs},
                        {"not_remixed_not_ignored_outputs", notRemixedNotIgnoredOutputs},
                    };
                }

                return result;
            },
            "Compute anonymity degradation in coinjoins. Ignore remixes by default", pybind11::arg("start"), pybind11::arg("stop"),
            pybind11::arg("daysToConsider"), pybind11::arg("coinjoinType"), pybind11::arg("coinjoinSubType") = py::none(), pybind11::arg("ignoreNonStandardDenominations") = false,
            pybind11::arg("ignoreRemixes") = true, pybind11::arg("ww2DenomsBucket") = py::none(), pybind11::arg("falseCoinjoins") = py::none()
        );
    ;
}
