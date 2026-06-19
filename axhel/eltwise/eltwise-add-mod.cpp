// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-add-mod.hpp" 

#include "eltwise/eltwise-add-mod-native.hpp"
#include "eltwise/eltwise-add-mod-sve.hpp"

namespace unipi {
namespace axhel {


    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        
        //TODO: check

        #pragma unroll 4
        for(uint64_t i=0; i<n; ++i){
            const uint64_t op1_val = op1[i];
            const uint64_t op2_val = op2[i];

            uint64_t sum = op1_val + op2_val;

            sum -= static_cast<uint64_t>(sum >= mod) * mod;

            res[i] = sum;

            /*
            uint64_t sum = *op1 + *op2;
            if(sum >= mod) {
                *res = sum - mod;
            } else {
                *res = sum;
            }
            ++op1;
            ++op2;
            ++res;
            */
        }
    }


    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod) {
        
        //TODO: check

        //uint64_t diff = mod - op2;

        #pragma unroll 4
        for (uint64_t i = 0; i<n; ++i) {
            const uint64_t op1_val = op1[i];

            uint64_t sum = op1_val + op2;

            sum -= static_cast<uint64_t>(sum >= mod) * mod;

            res[i] = sum;
            /*
            if (*op1 >= diff) {
                *res = *op1 - diff;
            } else {
                *res = *op1 + op2;
            }
            ++op1;
            ++res;
            */
        }        
    }


    void EltwiseAddMod(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        //TODO: check
       
        #if defined(AXHEL_HAS_SVE) || defined(AXHEL_HAS_SVE2)
        EltwiseAddModSVE(res,op1,op2,n,mod);
        #else
        EltwiseAddModNative(res,op1,op2,n,mod);
        #endif
    }


    void EltwiseAddMod(uint64_t* res, const uint64_t* op1, const uint64_t op2, uint64_t n, uint64_t mod) {
         //TODO: check

        #if defined(AXHEL_HAS_SVE) || defined(AXHEL_HAS_SVE2)
        EltwiseAddModSVE(res,op1,op2,n,mod);
        #else
        EltwiseAddModNative(res,op1,op2,n,mod);
        #endif
    }
}
}