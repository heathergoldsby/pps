/* digital_evolution/ancestors/multibirth_not_nand_prop_ancestor.h
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
#ifndef _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_NAND_ANCESTOR2_H_
#define _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_NAND_ANCESTOR2_H_


#include <ea/digital_evolution.h>
#include <ea/metadata.h>


namespace ealib {
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_not_nand_prop_ancestor2 {
        template <typename EA>
        typename EA::genome_type operator()(EA& ea) {
            typename EA::genome_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            //            repr[40] = ea.isa()["create_propagule"]; // create propagule
            repr[49] = ea.isa()["if_prop_cell_absent"];
            repr[50] = ea.isa()["deploy_propagule"]; // deploy propagule!
            
            // nand
            repr[74] = ea.isa()["fixed_input"]; // input
            repr[75] = ea.isa()["nop_c"]; // nopc
            repr[76] = ea.isa()["fixed_input"]; // input
            repr[77] = ea.isa()["nand"]; // nand
            repr[78] = ea.isa()["output"]; //output
            repr[79] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
            repr[90] =  ea.isa()["rotate_cw"];
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            repr[96] =  ea.isa()["h_divide_soft_parent_reset"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
            
            return repr;
        }
    };
    
}

#endif
