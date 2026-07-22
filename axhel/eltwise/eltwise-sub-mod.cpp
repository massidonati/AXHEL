// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-sub-mod.hpp" 

#include "eltwise/eltwise-sub-mod-native.hpp"
#include "eltwise/eltwise-sub-mod-sve.hpp"

namespace unipi {
namespace axhel {


    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        
        //TODO: check

        #pragma unroll 4
        for(uint64_t i=0; i<n; ++i){
            const uint64_t op1_val = op1[i];
            const uint64_t op2_val = op2[i];

            uint64_t diff = op1_val - op2_val;

            diff += static_cast<uint64_t>(op1_val < op2_val) * mod;

            res[i] = diff;
            /*
            uint64_t op1_val = *op1;
            uint64_t op2_val = *op2;
            if(op1_val >= op2_val) {
                *res = op1_val - op2_val;
            } else {
                *res = op1_val + mod - op2_val;
            }
            ++op1;
            ++op2;
            ++res;
            */
        }
    }


    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod) {
        
        //TODO: check

        #pragma unroll 4
        for (uint64_t i = 0; i<n; ++i) {
            const uint64_t op1_val = op1[i];

            uint64_t diff = op1_val - op2;

            diff += static_cast<uint64_t>(op1_val < op2) * mod;

            res[i] = diff;

            /*
            uint64_t op1_val = *op1;

            if(op1_val >= op2) {
                *res = op1_val - op2;
            } else {
                *res = op1_val + mod - op2;
            }
            ++op1;
            ++res;
            */
        }        
    }


    void EltwiseSubMod(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        //TODO: check
       
        #ifdef AXHEL_HAS_SVE
        EltwiseSubModSVE(res,op1,op2,n,mod);
        #else
        EltwiseSubModNative(res,op1,op2,n,mod);
        #endif
    }


    void EltwiseSubMod(uint64_t* res, const uint64_t* op1, const uint64_t op2, uint64_t n, uint64_t mod) {
         //TODO: check

        #ifdef AXHEL_HAS_SVE
        EltwiseSubModSVE(res,op1,op2,n,mod);
        #else
        EltwiseSubModNative(res,op1,op2,n,mod);
        #endif
    }
}
}