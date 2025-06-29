//
//  tx_identification.hpp
//  blocksci
//
//  Created by Harry Kalodner on 12/1/17.
//

#ifndef tx_identification_hpp
#define tx_identification_hpp

#include <blocksci/blocksci_export.h>
#include <blocksci/chain/chain_fwd.hpp>
#include <blocksci/scripts/scripts_fwd.hpp>
#include <optional>
#include <string>

namespace blocksci {
    class DataAccess;
    namespace heuristics {


    // original BlockSci heuristics    
    enum class BLOCKSCI_EXPORT CoinJoinResult {True, False, Timeout};
    CoinJoinResult BLOCKSCI_EXPORT isPossibleCoinjoin(const Transaction &tx, int64_t minBaseFee, double percentageFee, size_t maxDepth);
    CoinJoinResult BLOCKSCI_EXPORT isCoinjoinExtra(const Transaction &tx, int64_t minBaseFee, double percentageFee, size_t maxDepth);

    bool BLOCKSCI_EXPORT isDeanonTx(const Transaction &tx);
    bool BLOCKSCI_EXPORT containsKeysetChange(const Transaction &tx);
    bool BLOCKSCI_EXPORT isChangeOverTx(const Transaction &tx);
    bool BLOCKSCI_EXPORT isPeelingChain(const Transaction &tx);


    enum class BLOCKSCI_EXPORT CoinJoinType { WW2zkSNACKs, WW2PostzkSNACKs, WW1, Whirlpool, None };
    bool BLOCKSCI_EXPORT isCoinjoin(const Transaction &tx);
    
    /**
     * Wasabi 2 CoinJoin detection, ported from Dumplings.
     * 2 addtional checks added:
     *  1. After June 30, 2024, lower input limit to 40 to account for emerging post-zksnacks coordinators
     *  2. We require at least 5 different input and 5 different output addresses in the transaction
     * 
     * @param tx The transaction to check
     * @param inputCount if provided, then the transaction must have this exact input count to be considered
     *                a Wasabi 2 CoinJoin. Used for collecting transactions around a specific timeframe.
     * 
     * @return true if the transaction is a Wasabi 2 CoinJoin, false otherwise
     */
    bool BLOCKSCI_EXPORT isWasabi2CoinJoin(const Transaction &tx, std::optional<uint64_t> inputCount = std::nullopt, bool debug = false);
    
    /** 
     * Wasabi 1 CoinJoin detection, ported from Dumplings.
     * 
     * @param tx The transaction to check
     * 
     * @return true if the transaction is a Wasabi 1 CoinJoin, false otherwise
     */
    bool BLOCKSCI_EXPORT isWasabi1CoinJoin(const Transaction &tx);

    /**
     * Whirlpool CoinJoin detection, ported from Dumplings.
     * 
     * @param tx The transaction to check
     * 
     * @return true if the transaction is a Whirlpool CoinJoin, false otherwise
     */
    bool BLOCKSCI_EXPORT isWhirlpoolCoinJoin(const Transaction &tx);
    /**
     * Given a string "wasabi2", "wasabi1", or "whirlpool", return whether 
     * the transaction is the coinjoin of the given type.
     * 
     * @param tx The transaction to check
     * @param type The type of coinjoin to check for (wasabi1, wasabi2, whirlpool)
     * 
     * @return true if the transaction is a coinjoin of the given type, false otherwise
     */
    bool BLOCKSCI_EXPORT isCoinjoinOfGivenType(const Transaction &tx, const std::string &type, std::optional<std::string> subtype = std::nullopt);

    /**
     * Given a transaction, return the type of coinjoin it is.
     * Separate values for "No CoinJoin" and "postzkSNACKs Wasabi 2" are provided.
     * 
     * @param tx The transaction to check
     * 
     * @return the type of coinjoin the transaction is.
     */
    CoinJoinType BLOCKSCI_EXPORT getCoinjoinTag(const Transaction &tx);
    

    // consolidation
    enum class BLOCKSCI_EXPORT ConsolidationType {
        Certain, Possible, BigRoundOutput, None
    };

    /**
     * Check if the given transaction is some sort of consolidation.
     * If the consolidation is certain (e.g. multiple inputs to 1 output), inputOutputRatio is ignored.
     * 
     * @param tx The transaction to check
     * @param inputOutputRatio The ratio of inputs to outputs to consider a transaction a consolidation
     * 
     * @return the type of consolidation the transaction is.
     */
    ConsolidationType BLOCKSCI_EXPORT getConsolidationType(const Transaction &tx, double inputOutputRatio = 2);

    enum class BLOCKSCI_EXPORT HWWalletRemixResult {Trezor, SW, False};

    /**
     * Hugely heuristical check for the type of wallet (hw or sw) that is mixing coins.
     * Based around the ideas of how the HW (trezor) and SW wallets operate, we can make 
     * assumptions about transactions around CoinJoins.
     * For example:
     * - transaction tx has oldest input older than 6 months, moves the funds to a single output
     * - the output is then an input to a mix
     * - after the mix, the funds are dormant for more than 6 months
     * - then we say it's most likely a HW wallet.
     */
    HWWalletRemixResult BLOCKSCI_EXPORT isLongDormantInRemixes(const Transaction &tx);
}} // namespace blocksci

#endif /* tx_identification_hpp */
