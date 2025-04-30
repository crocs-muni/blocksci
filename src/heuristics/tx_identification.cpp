//
//  tx_identification.cpp
//  blocksci
//
//  Created by Harry Kalodner on 12/1/17.
//

#include <algorithm>
#include <blocksci/chain/block.hpp>
#include <blocksci/chain/coinjoin_utils.hpp>
#include <blocksci/chain/input.hpp>
#include <blocksci/chain/output.hpp>
#include <blocksci/chain/transaction.hpp>
#include <blocksci/heuristics/tx_identification.hpp>
#include <blocksci/scripts/script_variant.hpp>
#include <iostream>
#include <numeric>
#include <optional>
#include <range/v3/range_for.hpp>
#include <unordered_map>
#include <unordered_set>

namespace blocksci {
namespace heuristics {


    // Peeling chains have one input and two outputs
    bool looksLikePeelingChain(const Transaction &tx) { return (tx.outputCount() == 2 && tx.inputCount() == 1); }

    // A transaction is considered a peeling chain if it has one input and two outputs,
    // and either the previous or one of the next transactions looks like a peeling chain.
    bool isPeelingChain(const Transaction &tx) {
        if (!looksLikePeelingChain(tx)) {
            return false;
        }
        // Check if past transaction is peeling chain
        if (looksLikePeelingChain(tx.inputs()[0].getSpentTx())) {
            return true;
        }
        // Check if any spending transaction is peeling chain
        RANGES_FOR(auto output, tx.outputs()) {
            if (output.isSpent() && looksLikePeelingChain(*output.getSpendingTx())) {
                return true;
            }
        }
        return false;
    }

    bool isCoinjoin(const Transaction &tx) {
        return isWhirlpoolCoinJoin(tx) || isWasabi2CoinJoin(tx) || isWasabi1CoinJoin(tx);
    }

    struct OutputBucket {
        int64_t currentValue;
        int64_t goalValue;

        explicit OutputBucket(int64_t goal) : currentValue(0), goalValue(goal) {}

        bool isFull() const { return currentValue >= goalValue; }

        void addInput(int64_t val) { currentValue += val; }

        void removeInput(int64_t val) { currentValue -= val; }

        int64_t remaining() const {
            if (goalValue > currentValue) {
                return goalValue - currentValue;
            } else {
                return 0;
            }
        }

        bool operator<(const OutputBucket &other) const { return remaining() < other.remaining(); }
    };

    CoinJoinResult _getSumCount(std::vector<int64_t> &values, std::vector<OutputBucket> buckets,
                                int64_t totalRemaining, int64_t valueLeft, size_t maxDepth, size_t &depth);
    CoinJoinResult getSumCount(std::vector<int64_t> &values, std::vector<int64_t> bucketGoals, size_t maxDepth);

    CoinJoinResult _getSumCount(std::vector<int64_t> &values, std::vector<OutputBucket> buckets,
                                int64_t totalRemaining, int64_t valueLeft, size_t maxDepth, size_t &depth) {
        if (totalRemaining > valueLeft) {
            return CoinJoinResult::False;
        }

        buckets.erase(std::remove_if(buckets.begin(), buckets.end(), [&](auto &bucket) { return bucket.isFull(); }),
                        buckets.end());
        if (buckets.size() == 0) {
            return CoinJoinResult::True;
        }

        if (values.size() == 0) {
            return CoinJoinResult::False;
        }

        depth++;

        if (maxDepth != 0 && depth > maxDepth) {
            return CoinJoinResult::Timeout;
        }

        std::sort(buckets.rbegin(), buckets.rend());

        int64_t lastValue = values.back();
        values.pop_back();
        valueLeft -= lastValue;
        for (auto &bucket : buckets) {
            int64_t remaining = totalRemaining - bucket.remaining();
            bucket.addInput(lastValue);
            remaining += bucket.remaining();
            auto res = _getSumCount(values, buckets, remaining, valueLeft, maxDepth, depth);
            if (res != CoinJoinResult::False) {
                return res;
            }
            bucket.removeInput(lastValue);
        }
        values.push_back(lastValue);
        return CoinJoinResult::False;
    }

    CoinJoinResult getSumCount(std::vector<int64_t> &values, std::vector<int64_t> bucketGoals, size_t maxDepth) {
        std::sort(values.begin(), values.end());
        std::sort(bucketGoals.rbegin(), bucketGoals.rend());

        int64_t valueLeft = 0;
        for (auto value : values) {
            valueLeft += value;
        }

        int64_t totalRemaining = 0;
        for (auto goal : bucketGoals) {
            totalRemaining += goal;
        }

        std::vector<OutputBucket> buckets;
        buckets.reserve(bucketGoals.size());
        for (auto goal : bucketGoals) {
            buckets.emplace_back(goal);
        }
        size_t depth = 0;
        return _getSumCount(values, buckets, totalRemaining, valueLeft, maxDepth, depth);
    }

    CoinJoinResult isCoinjoinExtra(const Transaction &tx, int64_t minBaseFee, double percentageFee,
                                    size_t maxDepth) {
        if (tx.inputCount() < 2 || tx.outputCount() < 3) {
            return CoinJoinResult::False;
        }

        uint16_t participantCount = (tx.outputCount() + 1) / 2;
        if (participantCount > tx.inputCount()) {
            return CoinJoinResult::False;
        }

        std::unordered_map<Address, int64_t> inputValues;
        RANGES_FOR(auto input, tx.inputs()) { inputValues[input.getAddress()] += input.getValue(); }

        if (participantCount > inputValues.size()) {
            return CoinJoinResult::False;
        }

        std::unordered_map<int64_t, std::unordered_set<Address>> outputValues;
        RANGES_FOR(auto output, tx.outputs()) { outputValues[output.getValue()].insert(output.getAddress()); }

        using pair_type = decltype(outputValues)::value_type;
        auto pr = std::max_element(
            std::begin(outputValues), std::end(outputValues),
            [](const pair_type &p1, const pair_type &p2) { return p1.second.size() < p2.second.size(); });

        if (pr->second.size() != participantCount) {
            return CoinJoinResult::False;
        }

        if (pr->first == 546 || pr->first == 2730) {
            return CoinJoinResult::False;
        }

        std::vector<int64_t> values;
        values.reserve(inputValues.size());
        for (auto &pair : inputValues) {
            values.push_back(pair.second);
        }

        int64_t goalValue = pr->first;

        int64_t maxPossibleFee = std::max(minBaseFee, static_cast<int64_t>(goalValue * percentageFee));

        std::vector<int64_t> bucketGoals;
        for (uint16_t i = 0; i < participantCount; i++) {
            bucketGoals.push_back(goalValue);
        }

        size_t j = 0;
        RANGES_FOR(auto output, tx.outputs()) {
            if (output.getValue() != goalValue) {
                bucketGoals[j] += output.getValue();
                j++;
            }
        }

        for (auto &goal : bucketGoals) {
            if (maxPossibleFee > goal) {
                goal = 0;
            } else {
                goal -= maxPossibleFee;
            }
        }

        return getSumCount(values, bucketGoals, maxDepth);
    }

    CoinJoinResult isPossibleCoinjoin(const Transaction &tx, int64_t minBaseFee, double percentageFee,
                                        size_t maxDepth) {
        if (tx.outputCount() == 1 || tx.inputCount() == 1) {
            return CoinJoinResult::False;
        }

        std::unordered_map<int64_t, uint16_t> outputValues;
        RANGES_FOR(auto output, tx.outputs()) { outputValues[output.getValue()]++; }

        using pair_type = decltype(outputValues)::value_type;
        auto pr = std::max_element(std::begin(outputValues), std::end(outputValues),
                                    [](const pair_type &p1, const pair_type &p2) { return p1.second < p2.second; });

        // There must be at least two outputs of equal value to create an anonymity set
        if (pr->second == 1) {
            return CoinJoinResult::False;
        }

        std::unordered_map<Address, int64_t> inputValues;
        RANGES_FOR(auto input, tx.inputs()) { inputValues[input.getAddress()] += input.getValue(); }

        if (inputValues.size() == 1) {
            return CoinJoinResult::False;
        }

        std::vector<Output> unknownOutputs;
        RANGES_FOR(auto output, tx.outputs()) {
            if (inputValues.find(output.getAddress()) == inputValues.end()) {
                unknownOutputs.push_back(output);
            }
        }

        if (unknownOutputs.size() <= 1) {
            return CoinJoinResult::False;
        }

        outputValues.clear();
        for (auto &output : unknownOutputs) {
            outputValues[output.getValue()]++;
        }
        pr = std::max_element(std::begin(outputValues), std::end(outputValues),
                                [](const pair_type &p1, const pair_type &p2) { return p1.second < p2.second; });
        // There must be at least two outputs of equal value to create an anonymity set
        if (pr->second == 1) {
            return CoinJoinResult::False;
        }

        std::vector<int64_t> values;
        values.reserve(inputValues.size());
        for (auto &pair : inputValues) {
            values.push_back(pair.second);
        }

        int64_t maxPossibleFee = std::max(minBaseFee, static_cast<int64_t>(pr->first * percentageFee));
        int64_t goalValue = 0;
        if (pr->first > goalValue) {
            goalValue = pr->first - maxPossibleFee;
        }

        std::vector<int64_t> bucketGoals = {goalValue, goalValue};

        return getSumCount(values, bucketGoals, maxDepth);
    }

    bool isDeanonTx(const Transaction &tx) {
        if (tx.isCoinbase()) {
            return false;
        }

        if (tx.outputCount() == 1) {
            return false;
        }

        std::unordered_set<AddressType::Enum> inputCounts;
        RANGES_FOR(auto input, tx.inputs()) { inputCounts.insert(input.getAddress().type); }

        if (inputCounts.size() != 1) {
            return false;
        }

        AddressType::Enum inputType = *inputCounts.begin();

        bool seenType = false;
        RANGES_FOR(auto output, tx.outputs()) {
            if (output.getType() == inputType) {
                if (seenType) {
                    return false;
                }
                seenType = true;
            }
        }

        return seenType;
    }

    namespace {
        ranges::optional<Address> getInsidePointer(const ranges::optional<Address> &address, DataAccess &access);

        ranges::optional<Address> getInsidePointer(const Address &pointer, DataAccess &access) {
            if (pointer.type == AddressType::Enum::SCRIPTHASH) {
                script::ScriptHash scriptHashAddress(pointer.scriptNum, access);
                return getInsidePointer(scriptHashAddress.getWrappedAddress(), access);
            } else if (pointer.type == AddressType::Enum::WITNESS_SCRIPTHASH) {
                script::WitnessScriptHash witnessScriptHashAddress(pointer.scriptNum, access);
                return getInsidePointer(witnessScriptHashAddress.getWrappedAddress(), access);
            } else {
                return pointer;
            }
        }

        ranges::optional<Address> getInsidePointer(const ranges::optional<Address> &address, DataAccess &access) {
            if (address) {
                return getInsidePointer(*address, access);
            } else {
                return ranges::nullopt;
            }
        }
    }  // namespace

    struct DetailedType {
        AddressType::Enum mainType;
        bool hasSubtype;
        AddressType::Enum subType;
        int i;
        int j;

        DetailedType(const Address &pointer, DataAccess &access)
            : mainType(pointer.type), hasSubtype(false), subType(AddressType::Enum::NONSTANDARD), i(0), j(0) {
            auto insidePointer = getInsidePointer(pointer, access);
            if (insidePointer) {
                subType = insidePointer->type;
                hasSubtype = true;
                if (subType == AddressType::Enum::MULTISIG) {
                    script::Multisig multisigAddress(insidePointer->scriptNum, access);
                    i = multisigAddress.getRequired();
                    j = static_cast<int>(multisigAddress.getTotal());
                }
            }
        }

        bool operator==(const DetailedType &other) const {
            if (mainType != other.mainType || subType != other.subType) {
                return false;
            }
            if (mainType == AddressType::Enum::SCRIPTHASH && (!hasSubtype || !other.hasSubtype)) {
                return false;
            }
            if (mainType == AddressType::Enum::WITNESS_SCRIPTHASH && (!hasSubtype || !other.hasSubtype)) {
                return false;
            }
            if (subType == AddressType::Enum::MULTISIG && (i != other.i || j != other.j)) {
                return false;
            }

            return true;
        }

        bool operator!=(const DetailedType &other) const { return !operator==(other); }
    };

    struct DetailedTypeHasher {
        size_t operator()(const DetailedType &b) const {
            std::size_t seed = 123945432;
            hash_combine(seed, b.mainType);
            if (b.hasSubtype) {
                hash_combine(seed, b.subType);
            }
            if (b.subType == AddressType::Enum::MULTISIG) {
                hash_combine(seed, b.i);
                hash_combine(seed, b.j);
            }
            return seed;
        }
    };

    bool isChangeOverTx(const Transaction &tx) {
        if (tx.isCoinbase()) {
            return false;
        }

        std::unordered_set<DetailedType, DetailedTypeHasher> outputTypes;
        RANGES_FOR(auto output, tx.outputs()) {
            outputTypes.insert(DetailedType{output.getAddress(), tx.getAccess()});
        }

        if (outputTypes.size() != 1) {
            return false;
        }

        std::unordered_set<DetailedType, DetailedTypeHasher> inputTypes;
        RANGES_FOR(auto input, tx.inputs()) { inputTypes.insert(DetailedType{input.getAddress(), tx.getAccess()}); }

        if (inputTypes.size() != 1) {
            return false;
        }

        return *outputTypes.begin() != *inputTypes.begin();
    }

    bool containsKeysetChange(const Transaction &tx) {
        if (tx.isCoinbase()) {
            return false;
        }

        std::unordered_set<Address> multisigOutputs;
        RANGES_FOR(auto output, tx.outputs()) {
            auto pointer = getInsidePointer(output.getAddress(), tx.getAccess());
            if (pointer && pointer->type == AddressType::Enum::MULTISIG) {
                multisigOutputs.insert(*pointer);
            }
        }

        if (multisigOutputs.size() == 0) {
            return false;
        }

        std::unordered_set<Address> multisigInputs;
        RANGES_FOR(auto input, tx.inputs()) {
            auto pointer = getInsidePointer(input.getAddress(), tx.getAccess());
            if (pointer && pointer->type == AddressType::Enum::MULTISIG) {
                if (multisigOutputs.find(*pointer) == multisigOutputs.end()) {
                    multisigInputs.insert(*pointer);
                }
            }
        }

        if (multisigInputs.size() == 0) {
            return false;
        }

        std::unordered_set<Address> containedOutputs;
        for (auto &pointer : multisigOutputs) {
            std::function<void(const blocksci::Address &)> visitFunc = [&](const blocksci::Address &add) {
                containedOutputs.insert(add);
            };
            pointer.getScript().visitPointers(visitFunc);
        }

        for (auto &pointer : multisigInputs) {
            bool foundMatch = false;
            std::function<void(const blocksci::Address &)> visitFunc = [&](const blocksci::Address &add) {
                if (containedOutputs.find(add) != containedOutputs.end()) {
                    foundMatch = true;
                    return;
                }
            };
            pointer.getScript().visitPointers(visitFunc);
            if (foundMatch) {
                return true;
            }
        }

        return false;
    }

    void findSubsets(std::vector<int>::iterator start, std::vector<int>::iterator end, int sum,
                        std::vector<int> &path, int target) {
        if (sum == target) {
            return;
        }

        if (sum > target || start == end || path.size() >= 4) {
            return;
        }

        path.push_back(*start);
        findSubsets(std::next(start), end, sum + *start, path, target);

        path.pop_back();
        findSubsets(std::next(start), end, sum, path, target);
    }

    void subsetSum(std::vector<int> &numbers, int target, std::vector<int> &path) {
        findSubsets(numbers.begin(), numbers.end(), 0, path, target);
    }

    /**
     * tx must have an input older than 6 months
     * tx must have an output that is spent within 3 days (spendingTx)
     * spendingTx must go to a Wasabi2 CoinJoin (to_cj_spendingTx)
     * to_cj_spendingTx has some value, I need to find the subset sum in coinjoin outputs that equals
     * to_cj_spendingTx input
     */
    HWWalletRemixResult isLongDormantInRemixes(const Transaction &tx) {
        // the tx has to be dormant for 6 months and is >= 1 BTC (out -> hww)
        if (tx.outputCount() != 1) {
            return HWWalletRemixResult::False;
        }

        auto output = tx.outputs()[0];

        uint32_t tx_timestamp = tx.block().timestamp();
        uint32_t oldest_input = tx_timestamp;

        for (const auto &input : tx.inputs()) {
            if (input.getSpentTx().block().timestamp() < oldest_input) {
                oldest_input = input.getSpentTx().block().timestamp();
            }
        }

        if (tx_timestamp - oldest_input < 15814800) {  // 6 months
            return HWWalletRemixResult::False;
        }

        if (output.getValue() <= 100000000) {
            return HWWalletRemixResult::False;
        }
        // then it has to move to a new within 3 days and get remixed (hww -> cj space)

        if (!output.isSpent()) {
            return HWWalletRemixResult::False;
        }

        auto spendingTx = output.getSpendingTx().value();

        if (output.block().timestamp() - tx.block().timestamp() > 259200) {  // should be less than 3 days
            return HWWalletRemixResult::False;
        }

        for (const auto &o : spendingTx.outputs()) {
            // then to a mix (cj space -> mix)
            auto to_cj_output = o;
            auto value = o.getValue();
            if (!to_cj_output.isSpent()) {
                continue;
            }

            auto to_cj_spendingTx = to_cj_output.getSpendingTx().value();

            if (!blocksci::heuristics::isWasabi2CoinJoin(to_cj_spendingTx)) {
                continue;
            }

            std::vector<int> subsets_values = {}, path = {};
            for (const auto &o : to_cj_spendingTx.outputs()) {
                subsets_values.push_back(o.getValue());
            }

            subsetSum(subsets_values, value - value * 0.003, path);

            if (path.size() == 0) {
                continue;
            }

            // then if output stays for 6 months it's a HW trezor suite wallet remix (starting from 06/23)
            if (to_cj_spendingTx.block().timestamp() > 1685570400) {
                if (!to_cj_spendingTx.outputs()[0].isSpent() 
                || to_cj_spendingTx.outputs()[0].getSpendingTx().value().block().timestamp() - to_cj_spendingTx.block().timestamp() < 15814800) {
                    return HWWalletRemixResult::Trezor;
                }
            }

            // if output goes in 3 days and stays for 6 months it's a SW wallet remix
            if (to_cj_spendingTx.block().timestamp() - to_cj_output.getSpendingTx().value().block().timestamp() <
                259200) {
                if (!to_cj_spendingTx.outputs()[0].isSpent() ||
                    to_cj_spendingTx.outputs()[0].getSpendingTx().value().block().timestamp() -
                            to_cj_spendingTx.block().timestamp() <
                        15814800) {
                    return HWWalletRemixResult::SW;
                }
            }
        }
        // otherwise it's just false
        return HWWalletRemixResult::False;
    }

    /**
     * WW1 detection based on detection of the coordinator fee in the outputs.
     * The coordinator fee was sent to one of 2 fixed addresses.
     * All outputs must be segwit, there must be at least 2 different output values 
     * and each must appear at least twice.
     */
    bool _isWasabi1StaticCoordinator(const Transaction &tx) {
        std::unordered_map<int64_t, int> outputValues;
        bool is_coord_output = false;

        for (auto output : tx.outputs()) {
            outputValues[output.getValue()]++;
            if (output.getType() != AddressType::Enum::WITNESS_PUBKEYHASH) {
                return false;
            }
            if (!is_coord_output) {
                auto address = output.getAddress().toString();
                if (address == "WitnessPubkeyAddress(bc1qs604c7jv6amk4cxqlnvuxv26hv3e48cds4m0ew)" ||
                    address == "WitnessPubkeyAddress(bc1qa24tsgchvuxsaccp8vrnkfd85hrcpafg20kmjw)") {
                    is_coord_output = true;
                }
            }
        }
        
        if (!is_coord_output) {
            return false;
        }

        // there are at least 2 different outputs and each appeared at least twice
        if (outputValues.size() < 2) {
            return false;
        }

        for (auto &output : outputValues) {
            if (output.second > 2) {
                return true;
            }
        }

        return false;
    }

    /**
     * After the removing of the static coordinators, we need to employ more "heuristical"
     * methods to detect Wasabi1 CoinJoins.
     * - The base denomination is fixed at 0.1 BTC +- 0.02 BTC
     * - There are at least 2 different output values
     * - The most frequent output value appears at least 10 times
     * - The most frequent output value is almost the base denomination
     * - There are more inputs than the most frequent output value
     * - All outputs are segwit
     */
    bool _isWasabi1PostStaticCoordinator(const Transaction &tx) {                                
        std::unordered_map<int64_t, int> outputValues;
        for (auto output : tx.outputs()) {
            // all inputs segwit only
            if (output.getType() != AddressType::Enum::WITNESS_PUBKEYHASH) {
                return false;
            }
            outputValues[output.getValue()]++;
        }

        auto uniqueOutputValues = std::count_if(std::begin(outputValues), std::end(outputValues),
                                                [](const std::pair<int64_t, int> &p) { return p.second == 1; });

        if (uniqueOutputValues < 2) {
            return false;
        }

        auto mostFrequentEqualOutputCount =
            std::max_element(std::begin(outputValues), std::end(outputValues),
                                [](const std::pair<int64_t, int> &p1, const std::pair<int64_t, int> &p2) {
                                    return p1.second < p2.second;
                                });

        if (mostFrequentEqualOutputCount->second < 10) {
            return false;
        }

        // and there are more inputs than most frequent equal outputs
        if (tx.inputCount() < mostFrequentEqualOutputCount->second) {
            return false;
        }

        // and the most frequent equal output is almost the base denomination 0.1 BTC +- 0.02 BTC

        if (mostFrequentEqualOutputCount->first < 0.08 * 1e8 || mostFrequentEqualOutputCount->first > 0.12 * 1e8) {
            return false;
        }

        return true;
    }

    /**
     * Wasabi1 CoinJoin detection heuristic.
     * If the transaction happened before a certain block, different heuristics must be used.
     */
    bool isWasabi1CoinJoin(const Transaction &tx) {
        if (tx.getBlockHeight() < blocksci::CoinjoinUtils::FirstWasabiBlock) {
            return false;
        }
        // before FirstWasabiNoCoordAddressBlock WW1 had different base denominations and fixed coordinators
        if (tx.getBlockHeight() < blocksci::CoinjoinUtils::FirstWasabiNoCoordAddressBlock) {
            return _isWasabi1StaticCoordinator(tx);
        }

        if (blocksci::heuristics::isWasabi2CoinJoin(tx)) {
            return false;
        }
        // after FirstWasabiNoCoordAddressBlock WW1 had fixed base denomination and no coordinators
        return _isWasabi1PostStaticCoordinator(tx);
    }

    /**
     * Wasabi2 CoinJoin detection heuristic.
     * Ported from Dumplings
     */
    bool isWasabi2CoinJoin(const Transaction &tx, std::optional<uint64_t> inputCount, bool debug) {
        // first ww2 coinjoin block

        auto blockHeight = tx.getBlockHeight();
        if (blockHeight < blocksci::CoinjoinUtils::FirstWasabi2Block) {
            return false;
        }
        auto inputLimit = blockHeight < 850237 ? 50 : 20;
        for (const auto &input : tx.inputs()) {
            if (input.getType() != AddressType::Enum::WITNESS_PUBKEYHASH &&
                input.getType() != AddressType::Enum::WITNESS_UNKNOWN) {
                if (debug) {
                    std::cout << "Input type: " << input.getType() << std::endl;
                }
                return false;
            }
        }
        for (const auto &output : tx.outputs()) {
            if (output.getType() != AddressType::Enum::WITNESS_PUBKEYHASH &&
                output.getType() != AddressType::Enum::WITNESS_UNKNOWN) {
                if (debug) {
                    std::cout << "Output type: " << output.getType() << std::endl;
                }
                return false;
            }
        }

        if ((inputCount.has_value() && tx.inputCount() != inputCount) ||
            (!inputCount.has_value() && tx.inputCount() < inputLimit)) {
            if (debug) {
                std::cout << "Input count: " << tx.inputCount() << std::endl;
            }
            return false;
        }

        // Inputs are ordered descending and there are more than 3 input addresses.
        long long prev_input_value = -1;
        std::unordered_set<Address> usedAddresses;
        std::unordered_set<std::string> usedWitnessUnknownAddresses;
        for (const auto &input : tx.inputs()) {
            if (prev_input_value != -1 && input.getValue() > prev_input_value) {
                if (debug) {
                    std::cout << "Input value: " << input.getValue() << std::endl;
                }
                return false;
            }

            if (input.getAddress().getType() == AddressType::Enum::WITNESS_UNKNOWN) {
                script::WitnessUnknown witnessUnknownAddress(input.getAddress().scriptNum, tx.getAccess());
                auto *data = witnessUnknownAddress.getDataForMe();
                std::string stringData;
                stringData.resize(32);
                std::copy(data->scriptData.begin(), data->scriptData.end(), stringData.begin());

                usedWitnessUnknownAddresses.insert(stringData);

            } else {
                usedAddresses.insert(input.getAddress());
            }
            prev_input_value = input.getValue();
        }

        if (usedAddresses.size() + usedWitnessUnknownAddresses.size() < 5) {
            if (debug) {
                std::cout << "Input addresses: " << usedAddresses.size() << " + " << usedWitnessUnknownAddresses.size()
                          << std::endl;
            }
            return false;
        }

        usedAddresses.clear();
        usedWitnessUnknownAddresses.clear();

        // Outputs are ordered descending + there are more than 80% of outputs that are WW2 denominations
        long long prev_output_value = -1;
        int count = 0;
        for (const auto &output : tx.outputs()) {
            if (prev_output_value != -1 && output.getValue() > prev_output_value) {
                if (debug) {
                    std::cout << "Output value: " << output.getValue() << std::endl;
                }
                return false;
            }

            if (output.getAddress().getType() == AddressType::Enum::WITNESS_UNKNOWN) {
                script::WitnessUnknown witnessUnknownAddress(output.getAddress().scriptNum, tx.getAccess());
                auto *data = witnessUnknownAddress.getDataForMe();
                std::string stringData;
                stringData.resize(32);
                std::copy(data->scriptData.begin(), data->scriptData.end(), stringData.begin());

                usedWitnessUnknownAddresses.insert(stringData);
            } else {
                usedAddresses.insert(output.getAddress());
            }

            prev_output_value = output.getValue();
            if (CoinjoinUtils::ww2_denominations.find(output.getValue()) !=
                CoinjoinUtils::ww2_denominations.end()) {
                count++;
            }
        }
        if (usedAddresses.size() + usedWitnessUnknownAddresses.size() < 5) {
            if (debug) {
                std::cout << "Output addresses: " << usedAddresses.size() << " + " << usedWitnessUnknownAddresses.size()
                          << std::endl;
            }
            return false;
        }

        if (debug) {
            std::cout << "Output count: " << count << ", " << tx.outputCount() * (blockHeight < 850237 ? 0.8 : 0.75)<< std::endl;
        }
        return count > (tx.outputCount() * (blockHeight < 850237 ? 0.8 : 0.75));
    }

    /**
     * Whirlpool CoinJoin detection heuristic, ported from Dumplings.
     * The heuristic is based on the following rules:
     * - The current pool size must be one of the following: 0.001 BTC, 0.01 BTC, 0.1 BTC, 0.5 BTC
     * - The number of inputs must be between 5 and 10
     * - The number of outputs must be between 5 and 10 (but the same as the number of inputs)
     * - All outputs must have the same value as the current pool size
     * - There must be at least one input with the current pool size
     * - All inputs must have the current pool size or be within 0.0011 BTC of the current pool size (fee)
     */
    bool isWhirlpoolCoinJoin(const Transaction &tx) {
        if (tx.getBlockHeight() < blocksci::CoinjoinUtils::FirstSamouraiBlock) {
            return false;
        }

        auto current_pool_size = tx.outputs()[0].getValue();

        if (current_pool_size != 50000000 && current_pool_size != 5000000 && current_pool_size != 1000000 &&
            current_pool_size != 100000) {
            return false;
        }
        auto input_count = tx.inputCount();
        auto output_count = tx.outputCount();

        if (input_count < 5 || input_count > 10) {
            return false;
        }

        if (output_count < 5 || output_count > 10) {
            return false;
        }

        if (input_count != output_count) {
            return false;
        }

        std::unordered_map<int64_t, int> outputValues;
        RANGES_FOR(auto output, tx.outputs()) { outputValues[output.getValue()]++; }

        if (outputValues.size() != 1) {
            return false;
        }

        if (outputValues.begin()->first != current_pool_size) {
            return false;
        }

        auto poolSizedInputCount = std::count_if(
            tx.inputs().begin(), tx.inputs().end(),
            [current_pool_size](const Input &input) { return input.getValue() == current_pool_size; });

        if (poolSizedInputCount < 1) {
            return false;
        }

        //tx.Inputs.Where(x => x.PrevOutput.Value != poolSize).All(x => x.PrevOutput.Value.Almost(poolSize, Money.Coins(0.0011m)));
        return std::all_of(tx.inputs().begin(), tx.inputs().end(), [current_pool_size](const Input &input) {
            return input.getValue() == current_pool_size 
                || std::abs(input.getValue() - current_pool_size) < 110000;
        });
    }

    ConsolidationType getConsolidationType(const Transaction &tx, double inputOutputRatio) {
        if (tx.inputCount() != 1 && tx.outputCount() == 1) {
            return ConsolidationType::Certain;
        }

        auto total_output_value =
            std::accumulate(tx.outputs().begin(), tx.outputs().end(), 0,
                            [](int64_t sum, const Output &output) { return sum + output.getValue(); });

        // Check if there is one huge output that has the majority of the output value, while the rest are small,
        // it's a consolidation
        if (tx.outputCount() > 1 && tx.outputCount() < (tx.inputCount() * inputOutputRatio)) {
            // Check if any input is larger than half of the total output value
            if (std::any_of(tx.inputs().begin(), tx.inputs().end(), [total_output_value](const Input &input) {
                    return input.getValue() > total_output_value * 0.5;
                })) {
                return ConsolidationType::None;
            }

            // Determine the maximum output value
            auto max_output_value =
                std::max_element(tx.outputs().begin(), tx.outputs().end(),
                                    [](const Output &a, const Output &b) { return a.getValue() < b.getValue(); });

            // If the largest output is more than 90% of the total output value, it's a possible consolidation
            if ((*max_output_value).getValue() > total_output_value * 0.9) {
                return ConsolidationType::Possible;
            }
        }

        return ConsolidationType::None;
    }

    bool isCoinjoinOfGivenType(const Transaction &tx, const std::string &type, std::optional<std::string> subtype) {
        if (type == "wasabi1") {
            return isWasabi1CoinJoin(tx);
        } else if (type == "wasabi2") {
            return isWasabi2CoinJoin(tx);
        } else if (type == "whirlpool") {
            auto isWhirlpool = isWhirlpoolCoinJoin(tx);
            if (!isWhirlpool) {
                return false;
            }
            if (!subtype.has_value()) {
                return isWhirlpool;
            }
            // 0.001 BTC, 0.01 BTC, 0.05 BTC, and 0.5 BTC (to satoshis)
            if (subtype.value() == "50m") {
                return tx.inputs()[0].getValue() == 50000000;
            } else if (subtype.value() == "5m") {
                return tx.inputs()[0].getValue() == 5000000;
            } else if (subtype.value() == "1m") {
                return tx.inputs()[0].getValue() == 1000000;
            } else if (subtype.value() == "100k") {
                return tx.inputs()[0].getValue() == 100000;
            } else {
                return false;
            }

        } else {
            return false;
        }
    }

    CoinJoinType getCoinjoinTag(const Transaction &tx) {
        if (blocksci::heuristics::isWasabi2CoinJoin(tx)) {
            // 850237 is July 1st 2024
            return tx.block().height() < 850237 ? blocksci::heuristics::CoinJoinType::WW2zkSNACKs : blocksci::heuristics::CoinJoinType::WW2PostzkSNACKs;
        } else if (blocksci::heuristics::isWhirlpoolCoinJoin(tx)) {
            return blocksci::heuristics::CoinJoinType::Whirlpool;
        } else if (blocksci::heuristics::isWasabi1CoinJoin(tx)) {
            return blocksci::heuristics::CoinJoinType::WW1;
        }  else {
            return blocksci::heuristics::CoinJoinType::None;
        }
    }
}  // namespace heuristics
}  // namespace blocksci
