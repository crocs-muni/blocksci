//
//  cluster_py.cpp
//  blocksci
//
//  Created by Harry Kalodner on 1/14/17.
//  Copyright © 2017 Harry Kalodner. All rights reserved.
//

#include "cluster_py.hpp"

#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <blocksci/chain/block.hpp>
#include <blocksci/chain/blockchain.hpp>
#include <blocksci/cluster/cluster_manager.hpp>
#include <blocksci/cluster/coinjoin_cluster_manager.hpp>
#include <blocksci/cluster/coinjoin_clustering_heuristics.hpp>
#include <blocksci/heuristics/change_address.hpp>
#include <range/v3/range_for.hpp>

#include "caster_py.hpp"
#include "cluster_properties_py.hpp"
#include "proxy_utils.hpp"
#include "ranges_py.hpp"
#include "self_apply_py.hpp"

namespace py = pybind11;
using namespace blocksci;

int64_t totalOutWithoutSelfChurn(const blocksci::Block &block, blocksci::ClusterManager &manager) {
    int64_t total = 0;
    RANGES_FOR(auto tx, block) {
        std::set<uint32_t> inputClusters;
        RANGES_FOR(auto input, tx.inputs()) {
            auto cluster = manager.getCluster(input.getAddress());
            if (cluster.getTypeEquivSize() < 30000) {
                inputClusters.insert(cluster.clusterNum);
            }
        }
        RANGES_FOR(auto output, tx.outputs()) {
            if ((!output.isSpent() || output.getSpendingTx()->getBlockHeight() - block.height() > 3) &&
                inputClusters.find(manager.getCluster(output.getAddress()).clusterNum) == inputClusters.end()) {
                total += output.getValue();
            }
        }
    }
    return total;
}

void init_cluster_manager(pybind11::module &s) {
    s.def("total_without_self_churn", totalOutWithoutSelfChurn);

    py::class_<ClusterManager>(s, "ClusterManager", "Class managing the cluster dat")
        .def(py::init(
            [](std::string arg, blocksci::Blockchain &chain) { return ClusterManager(arg, chain.getAccess()); }))
        .def_static(
            "create_clustering",
            [](const std::string &location, Blockchain &chain, BlockHeight start, BlockHeight stop,
               heuristics::ChangeHeuristic &heuristic, bool shouldOverwrite, bool ignoreCoinJoin) {
                py::scoped_ostream_redirect stream(std::cout, py::module::import("sys").attr("stdout"));
                if (stop == -1) {
                    stop = chain.size();
                }
                auto range = chain[{start, stop}];
                return ClusterManager::createClustering(range, heuristic, location, shouldOverwrite, ignoreCoinJoin);
            },
            py::arg("location"), py::arg("chain"), py::arg("start") = 0, py::arg("stop") = -1,
            py::arg("heuristic") = heuristics::ChangeHeuristic{heuristics::NoChange{}},
            py::arg("should_overwrite") = false, py::arg("ignore_coinjoin") = true)
        .def(
            "cluster_with_address",
            [](const ClusterManager &cm, const Address &address) -> Cluster { return cm.getCluster(address); },
            py::arg("address"), "Return the cluster containing the given address")
        .def(
            "clusters", [](ClusterManager &cm) -> Range<Cluster> { return {cm.getClusters()}; },
            "Get a list of all clusters (The list is lazy so there is no cost to calling this method)")
        .def(
            "tagged_clusters",
            [](ClusterManager &cm, const std::unordered_map<blocksci::Address, std::string> &tags)
                -> Iterator<TaggedCluster> { return cm.taggedClusters(tags); },
            py::arg("tagged_addresses"),
            "Given a dictionary of tags, return a list of TaggedCluster objects for any clusters containing tagged "
            "scripts");
}

void init_coinjoin_cluster_manager(pybind11::module &s) {
    py::class_<CoinjoinClusterManager>(s, "CoinjoinClusterManager", "Class managing the cluster dat for coinjoins")
        .def(py::init([](std::string arg, blocksci::Blockchain &chain) {
            return CoinjoinClusterManager(arg, chain.getAccess());
        }))
        .def(
            "cluster_with_address",
            [](const CoinjoinClusterManager &cm, const Address &address) -> Cluster { return cm.getCluster(address); },
            py::arg("address"), "Return the cluster containing the given address")
        .def(
            "clusters", [](CoinjoinClusterManager &cm) -> Range<Cluster> { return {cm.getClusters()}; },
            "Get a list of all clusters (The list is lazy so there is no cost to calling this method)")
        .def(
            "tagged_clusters",
            [](CoinjoinClusterManager &cm, const std::unordered_map<blocksci::Address, std::string> &tags)
                -> Iterator<TaggedCluster> { return cm.taggedClusters(tags); },
            py::arg("tagged_addresses"),
            "Given a dictionary of tags, return a list of TaggedCluster objects for any clusters containing tagged "
            "scripts")
        .def_static(
            "create_clustering",
            [](Blockchain &chain, BlockHeight start, BlockHeight stop, const std::string &outputPath, bool overwrite,
               std::string coinjoinType, const std::string &heuristic) {
                py::scoped_ostream_redirect stream(std::cout, py::module::import("sys").attr("stdout"));
                auto heuristicFunc = blocksci::coinjoin_heuristics::getClusteringHeuristic(heuristic);

                if (stop == -1) {
                    stop = chain.size();
                }
                auto range = chain[{start, stop}];
                return CoinjoinClusterManager::createClustering(range, heuristicFunc, outputPath, overwrite, coinjoinType);
            },
            py::arg("chain"), py::arg("start"), py::arg("stop"), py::arg("output_path"),
            py::arg("overwrite") = false, py::arg("coinjoin_type") = "None", py::arg("heuristic") = "None",
            "Creates a clustering of the blockchain using the given heuristic and saves it to the given output path.\n"
            "Possible heuristics are:\n"
            "- 'None': No heuristic is used\n"
            "- 'OneOutputConsolidation': Clusters are created based on transactions with multiple inputs and one output\n"
            "Possible coinjoin types are:\n"
            "- 'none': No coinjoin type is used\n"
            "- 'wasabi2': Wasabi CoinJoin type 2\n"
            "- 'wasabi1': Wasabi CoinJoin type 1\n"
            "- 'whirlpool': Whirlpool CoinJoin");
}

void init_cluster(py::class_<Cluster> &cl) {
    applyMethodsToSelf(cl, AddClusterMethods{});

    cl.def("__len__", &Cluster::getSize)
        .def("__hash__", [](const Cluster &cluster) { return cluster.clusterNum; })
        .def("txes", &Cluster::getTransactions, "Returns a list of all transactions involving this cluster")
        .def(
            "in_txes",
            [](Cluster &cluster) {
                pybind11::print("Warning: `in_txes` is deprecated. Use `input_txes` instead.");
                return cluster.getInputTransactions();
            },
            "Returns a list of all transaction where this cluster was an input")
        .def("input_txes", &Cluster::getInputTransactions,
             "Returns a list of all transaction where this cluster was an input")
        .def(
            "out_txes",
            [](Cluster &cluster) {
                pybind11::print("Warning: `out_txes` is deprecated. Use `output_txes` instead.");
                return cluster.getOutputTransactions();
            },
            "Returns a list of all transaction where this cluster was an output")
        .def("output_txes", &Cluster::getOutputTransactions,
             "Returns a list of all transaction where this cluster was an output");
}

void addClusterRangeMethods(RangeClasses<Cluster> &classes) {
    addAllRangeMethods(classes);
}
