// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <arm_sve.h>
#include <stdint.h>

#if defined(AXHEL_HAS_SVE) || defined(AXHEL_HAS_SVE2)

namespace unipi {
namespace axhel {

    void EltwiseAddModSVE(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod){
        
        //TODO: check

        const uint64_t lanes = svcntd();

        #pragma unroll 4
        for(uint64_t i=0; i<n; i+=lanes){
            svbool_t loop_pg = svwhilelt_b64(i, n);

            svuint64_t op1_seg = svld1_u64(loop_pg, op1 + i);
            svuint64_t op2_seg = svld1_u64(loop_pg, op2 + i);
            svuint64_t res_seg = svadd_u64_m(loop_pg, op1_seg, op2_seg);

            svbool_t ge_cond = svcmpge_n_u64(loop_pg, res_seg, mod);
            res_seg = svsub_n_u64_m(ge_cond, res_seg, mod);

            svst1_u64(loop_pg, res + i, res_seg);           
        }
    }

    void EltwiseAddModSVE(uint64_t* res, const uint64_t* op1, const uint64_t op2, uint64_t n, uint64_t mod){

        const uint64_t lanes = svcntd();
        const svuint64_t op2_vec = svdup_n_u64(op2);

        #pragma unroll 4
        for (uint64_t i = 0; i < n; i += lanes) {
            svbool_t loop_pg = svwhilelt_b64(i, n);

            svuint64_t op1_seg = svld1_u64(loop_pg, op1 + i);
            svuint64_t res_seg = svadd_u64_m(loop_pg, op1_seg, op2_vec);

            svbool_t ge_cond = svcmpge_n_u64(loop_pg, res_seg, mod);
            res_seg = svsub_n_u64_m(ge_cond, res_seg, mod);

            svst1_u64(loop_pg, res + i, res_seg);
        }
    }

}
}

#endif
