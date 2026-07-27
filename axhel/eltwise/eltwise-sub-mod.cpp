// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-sub-mod.hpp" 
#include "eltwise/eltwise-sub-mod-native.hpp"
#include "eltwise/eltwise-sub-mod-sve.hpp"
#include "axhel/util/compiler.hpp"
#include "axhel/util/debug.hpp"

namespace unipi {
namespace axhel {


    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        
        AXHEL_UNROLL(4)
        for(uint64_t i=0; i<n; ++i){
            const uint64_t op1_val = op1[i];
            const uint64_t op2_val = op2[i];

            uint64_t diff = op1_val - op2_val;

            diff += static_cast<uint64_t>(op1_val < op2_val) * mod;

            res[i] = diff;
        }
    }


    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod) {
        
        AXHEL_UNROLL(4)
        for (uint64_t i = 0; i<n; ++i) {
            const uint64_t op1_val = op1[i];

            uint64_t diff = op1_val - op2;

            diff += static_cast<uint64_t>(op1_val < op2) * mod;

            res[i] = diff;
        }        
    }


    void EltwiseSubMod(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseSubMod vector-vector -> SVE");
        EltwiseSubModSVE(res,op1,op2,n,mod);
        #else
        AXHEL_LOG("EltwiseSubMod vector-vector -> native");
        EltwiseSubModNative(res,op1,op2,n,mod);
        #endif
    }


    void EltwiseSubMod(uint64_t* res, const uint64_t* op1, const uint64_t op2, uint64_t n, uint64_t mod) {
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseSubMod vector-scalar -> SVE");
        EltwiseSubModSVE(res,op1,op2,n,mod);
        #else
        AXHEL_LOG("EltwiseSubMod vector-scalar -> native");
        EltwiseSubModNative(res,op1,op2,n,mod);
        #endif
    }
}
}