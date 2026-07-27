// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-add-mod.hpp" 
#include "eltwise/eltwise-add-mod-native.hpp"
#include "eltwise/eltwise-add-mod-sve.hpp"
#include "axhel/util/compiler.hpp"
#include "axhel/util/debug.hpp"

namespace unipi {
namespace axhel {


    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        
        AXHEL_UNROLL(4)
        for(uint64_t i=0; i<n; ++i){
            const uint64_t op1_val = op1[i];
            const uint64_t op2_val = op2[i];

            uint64_t sum = op1_val + op2_val;

            sum -= static_cast<uint64_t>(sum >= mod) * mod;

            res[i] = sum;
        }
    }


    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod) {
        
        AXHEL_UNROLL(4)
        for (uint64_t i = 0; i<n; ++i) {
            const uint64_t op1_val = op1[i];

            uint64_t sum = op1_val + op2;

            sum -= static_cast<uint64_t>(sum >= mod) * mod;

            res[i] = sum;
        }        
    }


    void EltwiseAddMod(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {       
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseAddMod vector-vector -> SVE");
        EltwiseAddModSVE(res,op1,op2,n,mod);
        #else
        AXHEL_LOG("EltwiseAddMod vector-vector -> native");
        EltwiseAddModNative(res,op1,op2,n,mod);
        #endif
    }


    void EltwiseAddMod(uint64_t* res, const uint64_t* op1, const uint64_t op2, uint64_t n, uint64_t mod) {
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseAddMod vector-scalar -> SVE");
        EltwiseAddModSVE(res,op1,op2,n,mod);
        #else
        AXHEL_LOG("EltwiseAddMod vector-scalar -> native");
        EltwiseAddModNative(res,op1,op2,n,mod);
        #endif
    }
}
}