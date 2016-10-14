

#ifndef _THREE_STRIPE_ANCESTOR_H_
#define _THREE_STRIPE_ANCESTOR_H_


#include <ea/digital_evolution.h>
#include <ea/metadata.h>


namespace ealib {
    /*! Generates a self-replicating ancestor that performs
     */
    struct three_stripe_ancestor {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] = ea.isa()["h_alloc"];
            repr[1] = ea.isa()["push"];
            repr[2] = ea.isa()["nop_a"];
            repr[3] = ea.isa()["h_search"];
            repr[4] = ea.isa()["nop_c"];
            repr[5] = ea.isa()["mov_head"];
            repr[6] = ea.isa()["rx_msg"];
            repr[7] = ea.isa()["is_neighbor_matrix"];
            repr[8] = ea.isa()["h_alloc"];
            repr[9] = ea.isa()["h_alloc"];
            repr[10] = ea.isa()["donate_res_to_group"];
            repr[11] = ea.isa()["if_less"];
            repr[12] = ea.isa()["nop_x"];
            repr[13] = ea.isa()["nop_x"];
            repr[14] = ea.isa()["if_not_equal"];
            repr[15] = ea.isa()["tx_msg_check_task_matrix"];
            repr[16] = ea.isa()["h_divide_soft_parent_reset"];
            repr[17] = ea.isa()["rotate_cw"];
            repr[18] = ea.isa()["if_label"];
            repr[19] = ea.isa()["jump_head"];
            repr[20] = ea.isa()["rotate_ccw"];
            repr[21] = ea.isa()["is_neighbor_matrix"];
            repr[22] = ea.isa()["rx_msg"];
            repr[23] = ea.isa()["h_alloc"];
            repr[24] = ea.isa()["if_not_equal"];
            repr[25] = ea.isa()["fixed_input"];
            repr[26] = ea.isa()["push"];
            repr[27] = ea.isa()["nop_c"];
            repr[28] = ea.isa()["pop"];
            repr[29] = ea.isa()["nand"];
            repr[30] = ea.isa()["output"];
            repr[31] = ea.isa()["bc_msg_matrix"];
            repr[32] = ea.isa()["inc"];
            repr[33] = ea.isa()["push"];
            repr[34] = ea.isa()["tx_msg_check_task_matrix"];
            repr[35] = ea.isa()["h_search"];
            repr[36] = ea.isa()["h_alloc"];
            repr[37] = ea.isa()["h_search"];
            repr[38] = ea.isa()["jump_head"];
            repr[39] = ea.isa()["donate_res_to_group"];
            repr[40] = ea.isa()["h_divide_soft_parent_reset"];
            repr[41] = ea.isa()["dec"];
            repr[42] = ea.isa()["nop_x"];
            repr[43] = ea.isa()["pop"];
            repr[44] = ea.isa()["h_divide_soft_parent_reset"];
            repr[45] = ea.isa()["nand"];
            repr[46] = ea.isa()["nand"];
            repr[47] = ea.isa()["bc_msg_matrix"];
            repr[48] = ea.isa()["tx_msg_check_task_matrix"];
            repr[49] = ea.isa()["h_search"];
            repr[50] = ea.isa()["fixed_input"];
            repr[51] = ea.isa()["nop_c"];
            repr[52] = ea.isa()["fixed_input"];
            repr[53] = ea.isa()["nand"];
            repr[54] = ea.isa()["output"];
            repr[55] = ea.isa()["dec"];
            repr[56] = ea.isa()["nop_x"];
            repr[57] = ea.isa()["dec"];
            repr[58] = ea.isa()["rotate_cw"];
            repr[59] = ea.isa()["h_alloc"];
            repr[60] = ea.isa()["inc"];
            repr[61] = ea.isa()["nand"];
            repr[62] = ea.isa()["is_neighbor_matrix"];
            repr[63] = ea.isa()["rotate_cw"];
            repr[64] = ea.isa()["bc_msg_matrix"];
            repr[65] = ea.isa()["push"];
            repr[66] = ea.isa()["nand"];
            repr[67] = ea.isa()["if_less"];
            repr[68] = ea.isa()["if_not_equal"];
            repr[69] = ea.isa()["inc"];
            repr[70] = ea.isa()["rotate_cw"];
            repr[71] = ea.isa()["output"];
            repr[72] = ea.isa()["pop"];
            repr[73] = ea.isa()["push"];
            repr[74] = ea.isa()["nop_x"];
            repr[75] = ea.isa()["inc"];
            repr[76] = ea.isa()["fixed_input"];
            repr[77] = ea.isa()["nand"];
            repr[78] = ea.isa()["nand"];
            repr[79] = ea.isa()["output"];
            repr[80] = ea.isa()["nop_c"];
            repr[81] = ea.isa()["rotate"];
            repr[82] = ea.isa()["h_divide_soft_parent_reset"];
            repr[83] = ea.isa()["output"];
            repr[84] = ea.isa()["h_divide_soft_parent_reset"];
            repr[85] = ea.isa()["donate_res_to_group"];
            repr[86] = ea.isa()["nop_x"];
            repr[87] = ea.isa()["is_neighbor_matrix"];
            repr[88] = ea.isa()["rotate"];
            repr[89] = ea.isa()["tx_msg_check_task_matrix"];
            repr[90] = ea.isa()["is_neighbor_matrix"];
            repr[91] = ea.isa()["h_search"];
            repr[92] = ea.isa()["h_copy"];
            repr[93] = ea.isa()["h_copy"];
            repr[94] = ea.isa()["nop_a"];
            repr[95] = ea.isa()["if_label"];
            repr[96] = ea.isa()["h_divide_soft_parent_reset"];
            repr[97] = ea.isa()["mov_head"];
            repr[98] = ea.isa()["nop_a"];
            repr[99] = ea.isa()["nop_b"];
            
            return repr;
        }
    };
    
}

#endif
