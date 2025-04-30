//
//  tx_properties_py.hpp
//  blockscipy
//
//  Created by Malte Moeser on 8/26/19.
//

#ifndef tx_properties_py_h
#define tx_properties_py_h

#include "method_tags.hpp"

#include <blocksci/chain/algorithms.hpp>
#include <blocksci/chain/block.hpp>
#include <blocksci/heuristics/tx_identification.hpp>

#include <optional>
#include <pybind11/chrono.h>
#include <pybind11/operators.h>

#include "proxy.hpp"
#include "proxy_create.hpp"
#include "proxy_utils.hpp"

struct AddTransactionMethods {
    template <typename FuncApplication>
    void operator()(FuncApplication func) {
        using namespace blocksci;

        func(property_tag, "output_count", &Transaction::outputCount, "The number of outputs this transaction has");
        func(property_tag, "input_count", &Transaction::inputCount, "The number of inputs this transaction has");
        func(property_tag, "size_bytes", &Transaction::totalSize, "The size of this transaction in bytes");
        func(property_tag, "base_size", &Transaction::baseSize, "The size of the non-segwit data in bytes");
        func(property_tag, "total_size", &Transaction::totalSize, "The size all transaction data in bytes");
        func(property_tag, "virtual_size", &Transaction::virtualSize, "The weight of the transaction divided by 4");
        func(property_tag, "weight", &Transaction::weight, "Three times the base size plus the total size");
        func(property_tag, "locktime", &Transaction::locktime, "The locktime of this transasction");
        func(property_tag, "version", &Transaction::getVersion, "The version of this transaction");
        func(property_tag, "block_height", &Transaction::getBlockHeight, "The height of the block that this transaction was in");
        func(property_tag, "block_time", +[](const Transaction &tx) -> std::chrono::system_clock::time_point {
            return tx.block().getTime();
        }, "The time that the block containing this transaction arrived");
        func(property_tag, "observed_in_mempool", &Transaction::observedInMempool, "Returns whether this transaction was seen in the mempool by the recorder");
        func(property_tag, "time_seen", &Transaction::getTimeSeen, "If recorded by the mempool recorder, the time that this transaction was first seen by your node");
        func(property_tag, "timestamp_seen", &Transaction::getTimestampSeen, "If recorded by the mempool recorder, the time that this transaction was first seen by your node");
        func(property_tag, "block", &Transaction::block, "The block that this transaction was in");
        func(property_tag, "index", +[](const Transaction &tx) { return tx.txNum; }, "The internal index of this transaction");
        func(property_tag, "hash", &Transaction::getHash, "The 256-bit hash of this transaction");
        func(property_tag, "input_value", totalInputValue<Transaction &>, "The sum of the value of all of the inputs");
        func(property_tag, "output_value", totalOutputValue<Transaction &>, "The sum of the value of all of the outputs");
        func(property_tag, "fee", +[](const Transaction &tx) -> int64_t {
            return fee(tx);
        }, "The fee paid by this transaction");
        func(method_tag, "fee_per_byte", +[](const Transaction &tx, const std::string &sizeMeasure) -> int64_t {
            auto txFee = fee(tx);
            if (sizeMeasure == "total") {
                return txFee / tx.totalSize();
            } else if (sizeMeasure == "base") {
                return txFee / tx.baseSize();
            } else if(sizeMeasure == "weight") {
                return txFee / tx.weight();
            } else if(sizeMeasure == "virtual") {
                return txFee / tx.virtualSize();
            } else {
                throw std::invalid_argument{"Size measure must be one of total, base, weight, or virtual"};
            }
        }, "The (rounded) ratio of fee paid to size of this transaction (in byte). By default this uses virtual size, but passing total, base, weight, or virtual let's you choose a different size measure", pybind11::arg("size_measure") = "virtual");
        func(method_tag, "fee_per_kbyte", +[](const Transaction &tx, const std::string &sizeMeasure) -> int64_t {
            auto txFee = fee(tx);
            if (sizeMeasure == "total") {
                return txFee * 1000 / tx.totalSize();
            } else if (sizeMeasure == "base") {
                return txFee * 1000 / tx.baseSize();
            } else if(sizeMeasure == "weight") {
                return txFee * 1000 / tx.weight();
            } else if(sizeMeasure == "virtual") {
                return txFee * 1000 / tx.virtualSize();
            } else {
                throw std::invalid_argument{"Size measure must be one of total, base, weight, or virtual"};
            }
        }, "The (rounded) ratio of fee paid to size of this transaction (in kbyte). By default this uses virtual size, but passing total, base, weight, or virtual let's you choose a different size measure", pybind11::arg("size_measure") = "virtual");
        func(property_tag, "op_return", +[](const Transaction &tx) -> ranges::optional<Output> {
            return getOpReturn(tx);
        }, "If this transaction included a null data address, return its output. Otherwise return None");
        func(method_tag, "includes_output_of_type", includesOutputOfType, "Check whether the given transaction includes an output of the given address type", pybind11::arg("address_type"));
        func(property_tag, "is_coinbase", &Transaction::isCoinbase, "Return's true if this transaction is a Coinbase transaction");
        func(property_tag, "is_wasabi2_coinjoin", +[](const Transaction &tx) -> bool {
            return blocksci::heuristics::isWasabi2CoinJoin(tx);
        }, "Return's true if this transaction is a Wasabi CoinJoin transaction");
        func(property_tag, "is_ww2_cj_debug", +[](const Transaction &tx) -> bool {
            return blocksci::heuristics::isWasabi2CoinJoin(tx, std::nullopt, true);
        }, "Return's true if this transaction is a Wasabi CoinJoin transaction (debug mode)");
        func(property_tag, "is_whirlpool_coinjoin", +[](const Transaction &tx) -> bool {
            return blocksci::heuristics::isWhirlpoolCoinJoin(tx);
        }, "Return's true if this transaction is a Whirlpool CoinJoin transaction");
        func(property_tag, "is_wasabi1_coinjoin", +[](const Transaction &tx) -> bool {
            return blocksci::heuristics::isWasabi1CoinJoin(tx);
        }, "Return's true if this transaction is a WW1 CoinJoin transaction");
        func(method_tag, "is_wasabi2_coinjoin_with_input_count", +[](const Transaction &tx, uint64_t min_input_count) -> bool {
            return blocksci::heuristics::isWasabi2CoinJoin(tx, std::optional<uint64_t>{min_input_count}, false);
        }, "Return's true if this transaction is a CoinJoin transaction with given `min_input_count`", pybind11::arg("min_input_count") = 50); 
        func(property_tag, "mempool_space_link", +[](const Transaction &tx) -> std::string {
            return "https://mempool.space/tx/" + tx.getHash().GetHex();
        }, "A link to the mempool.space page for this transaction");
        func(property_tag, "coinjoin_tag", +[](const Transaction &tx) -> int64_t {
            if (blocksci::heuristics::isWasabi2CoinJoin(tx)) {
                // 850237 is July 1st 2024
                return tx.block().height() < 850237 ? static_cast<int64_t>(blocksci::heuristics::CoinJoinType::WW2zkSNACKs) : static_cast<int64_t>(blocksci::heuristics::CoinJoinType::WW2PostzkSNACKs);
            } else if (blocksci::heuristics::isWhirlpoolCoinJoin(tx)) {
                return static_cast<int64_t>(blocksci::heuristics::CoinJoinType::Whirlpool);
            } else if (blocksci::heuristics::isWasabi1CoinJoin(tx)) {
                return static_cast<int64_t>(blocksci::heuristics::CoinJoinType::WW1);
            }  else {
                return static_cast<int64_t>(blocksci::heuristics::CoinJoinType::None);
            }
        }, "The type of coinjoin this transaction is, if any");
    }
};


#endif /* tx_properties_py_h */
