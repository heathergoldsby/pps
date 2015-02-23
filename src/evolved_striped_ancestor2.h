/* evolved_striped_ancestor.h
 *
 * This file is part of EALib.
 *
 * Copyright 2014 David B. Knoester, Heather J. Goldsby.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef ps_evolved_striped_ancestor2_h
#define ps_evolved_striped_ancestor2_h

namespace ealib {
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct evolved_striped_ancestor2 {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] = ea.isa()["h_alloc"];
            repr[1] = ea.isa()["nop_c"];
            repr[2] = ea.isa()["nop_a"];
            repr[3] = ea.isa()["h_search"];
            repr[4] = ea.isa()["nop_c"];
            repr[5] = ea.isa()["mov_head"];
            repr[6] = ea.isa()["nop_x"];
            repr[7] = ea.isa()["h_search"];
            repr[8] = ea.isa()["pop"];
            repr[9] = ea.isa()["nop_x"];
            repr[10] = ea.isa()["deploy_propagule"]; // deploy propagule!
            repr[11] = ea.isa()["pop"];
            repr[12] = ea.isa()["nop_x"];
            repr[13] = ea.isa()["nop_x"];
            repr[14] = ea.isa()["h_alloc"];
            repr[15] = ea.isa()["rotate"];
            repr[16] = ea.isa()["nop_x"];
            repr[17] = ea.isa()["rotate_cw"];
            repr[18] = ea.isa()["inc"];
            repr[19] = ea.isa()["output"];
            repr[20] = ea.isa()["if_equal"];
            repr[21] = ea.isa()["nop_x"];
            repr[22] = ea.isa()["if_label"];
            repr[23] = ea.isa()["rx_msg"];
            repr[24] = ea.isa()["donate_res_to_group"];
            repr[25] = ea.isa()["fixed_input"];
            repr[26] = ea.isa()["push"];
            repr[27] = ea.isa()["nop_c"];
            repr[28] = ea.isa()["pop"];
            repr[29] = ea.isa()["nand"];
            repr[30] = ea.isa()["output"];
            repr[31] = ea.isa()["pop"];
            repr[32] = ea.isa()["inc"];
            repr[33] = ea.isa()["donate_res_to_group"];
            repr[34] = ea.isa()["rotate"];
            repr[35] = ea.isa()["if_equal"];
            repr[36] = ea.isa()["nop_x"];
            repr[37] = ea.isa()["if_label"];
            repr[38] = ea.isa()["inc"];
            repr[39] = ea.isa()["is_neighbor"];
            repr[40] = ea.isa()["h_alloc"];
            repr[41] = ea.isa()["nop_x"];
            repr[42] = ea.isa()["nand"];
            repr[43] = ea.isa()["if_label"];
            repr[44] = ea.isa()["donate_res_to_group"];
            repr[45] = ea.isa()["if_label"];
            repr[46] = ea.isa()["is_neighbor"];
            repr[47] = ea.isa()["h_copy"];
            repr[48] = ea.isa()["rx_msg"];
            repr[49] = ea.isa()["bc_msg"];
            repr[50] = ea.isa()["nop_x"];
            repr[51] = ea.isa()["nop_x"];
            repr[52] = ea.isa()["tx_msg"];
            repr[53] = ea.isa()["nop_x"];
            repr[54] = ea.isa()["if_label"];
            repr[55] = ea.isa()["tx_msg_check_task"];
            repr[56] = ea.isa()["tx_msg_check_task"];
            repr[57] = ea.isa()["if_label"];
            repr[58] = ea.isa()["h_copy"];
            repr[59] = ea.isa()["h_alloc"];
            repr[60] = ea.isa()["h_divide_soft_parent_reset"];
            repr[61] = ea.isa()["nop_b"];
            repr[62] = ea.isa()["nop_x"];
            repr[63] = ea.isa()["bc_msg"];
            repr[64] = ea.isa()["inc"];
            repr[65] = ea.isa()["nop_b"];
            repr[66] = ea.isa()["fixed_input"];
            repr[67] = ea.isa()["nop_a"];
            repr[68] = ea.isa()["if_label"];
            repr[69] = ea.isa()["rotate"];
            repr[70] = ea.isa()["nop_x"];
            repr[71] = ea.isa()["nand"];
            repr[72] = ea.isa()["rx_msg"];
            repr[73] = ea.isa()["is_neighbor"];
            repr[74] = ea.isa()["fixed_input"];
            repr[75] = ea.isa()["push"];
            repr[76] = ea.isa()["output"];
            repr[77] = ea.isa()["nand"];
            repr[78] = ea.isa()["output"];
            repr[79] = ea.isa()["nop_x"];
            repr[80] = ea.isa()["nop_x"];
            repr[81] = ea.isa()["is_neighbor"];
            repr[82] = ea.isa()["nop_c"];
            repr[83] = ea.isa()["push"];
            repr[84] = ea.isa()["jump_head"];
            repr[85] = ea.isa()["h_copy"];
            repr[86] = ea.isa()["if_not_equal"];
            repr[87] = ea.isa()["nop_x"];
            repr[88] = ea.isa()["rotate_ccw"];
            repr[89] = ea.isa()["nop_b"];
            repr[90] = ea.isa()["rotate_cw"];
            repr[91] = ea.isa()["h_search"];
            repr[92] = ea.isa()["h_copy"];
            repr[93] = ea.isa()["nop_c"];
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
