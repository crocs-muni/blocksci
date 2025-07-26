

#ifndef coinjoin_cluster_manager_hpp
#define coinjoin_cluster_manager_hpp

#include <blocksci/blocksci_export.h>

#include <blocksci/heuristics/tx_identification.hpp>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "../external/dset/dset.h"
#include "cluster.hpp"
#include "cluster_fwd.hpp"

namespace blocksci {
    namespace coinjoin_heuristics {
        struct ClusteringHeuristic;
    }

    struct BLOCKSCI_EXPORT AddressDisjointSets {
        DisjointSets disjoinSets;
        std::unordered_map<Address, uint32_t> addressIndices;

        AddressDisjointSets(uint32_t size, std::unordered_map<Address, uint32_t> addressIndices_)
            : disjoinSets(size), addressIndices(addressIndices_) {}

        uint32_t size() const { return disjoinSets.size(); }

        void link_addresses(const Address &address1, const Address &address2) {
            auto firstAddressIndex = addressIndices.at(address1);
            auto secondAddressIndex = addressIndices.at(address2);
            disjoinSets.unite(firstAddressIndex, secondAddressIndex);
        }

        void resolveAll();

        uint32_t find(uint32_t index) { return disjoinSets.find(index); }
    };

    class ClusterAccess;

    class BLOCKSCI_EXPORT CoinjoinClusterManager {
        std::unique_ptr<ClusterAccess> access;
        uint32_t clusterCount;

        friend class blocksci::Cluster;

       public:
        CoinjoinClusterManager(const std::string &baseDirectory, DataAccess &access);
        CoinjoinClusterManager(CoinjoinClusterManager &&other);
        CoinjoinClusterManager &operator=(CoinjoinClusterManager &&other);
        ~CoinjoinClusterManager();

        /**
         * Create a clustering around CoinJoin transactions of given type.
         *
         * First, find the CoinJoin transactions and collect addresses within 2 hops of them.
         * Then, create clusters using the provided clustering heuristic function.
         * Finally, serialize the cluster data to the output directory.
         *
         * @param chain BlockRange to cluster
         * @param clusteringFunc Clustering heuristic function
         * @param outputPath Path to output directory
         * @param overwrite Overwrite existing cluster data
         * @param coinjoinType Type of CoinJoin transactions to cluster (wasabi1, wasabi2, whirlpool, ashigaru)
         * @param maxHops Maximum number of hops to collect addresses around CoinJoin transactions
         * @return CoinjoinClusterManager instance
         */
        static CoinjoinClusterManager createClustering(
            BlockRange &chain, const blocksci::coinjoin_heuristics::ClusteringHeuristic &clusteringFunc,
            const std::string &outputPath, bool overwrite = false, std::string coinjoinType = "None", int maxHops = 2);

        Cluster getCluster(const Address &address) const;

        ranges::any_view<Cluster, ranges::category::random_access | ranges::category::sized> getClusters() const;

        ranges::any_view<TaggedCluster> taggedClusters(const std::unordered_map<Address, std::string> &tags) const;
    };

    using cluster_range = decltype(std::declval<CoinjoinClusterManager>().getClusters());
}  // namespace blocksci

#endif /* coinjoin_cluster_manager_hpp */
