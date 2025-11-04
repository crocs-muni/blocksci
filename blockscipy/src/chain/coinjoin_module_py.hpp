//
//  blockchain_py.hpp
//  blocksci
//
//  Created by Harry Kalodner on 4/30/18.
//

#ifndef coinjoin_module_py_hpp
#define coinjoin_module_py_hpp

#include <pybind11/pybind11.h>

#include <blocksci/chain/chain_fwd.hpp>

void init_coinjoin_module(pybind11::class_<blocksci::Blockchain> &cl);

#endif /* coinjoin_module_py_hpp */
