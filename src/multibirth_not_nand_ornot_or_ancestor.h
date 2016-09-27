//
//  multibirth_not_nand_ornot_or.h
//  ps
//
//  Created by Heather Goldsby on 9/27/16.
//  Copyright © 2016 Michigan State University. All rights reserved.
//

#ifndef _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_NAND_ORNOT_OR_ANCESTOR_H_
#define _EA_DIGITAL_EVOLUTION_ANCESTORS_MULTI_BIRTH_SELFREP_NOT_NAND_ORNOT_OR_ANCESTOR_H_


#include <ea/digital_evolution.h>
#include <ea/metadata.h>


namespace ealib {
    /*! Generates a self-replicating ancestor that performs 
     */
    struct multibirth_selfrep_not_nand_ornot_or_ancestor {
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
            repr[14] = ea.isa()["fixed_input"]; // input
            repr[15] = ea.isa()["fixed_input"]; // input
            repr[16] = ea.isa()["push"]; // push
            repr[17] = ea.isa()["nop_c"]; // nopc
            repr[18] = ea.isa()["pop"]; // pop
            repr[19] = ea.isa()["nand"]; // nand
            repr[20] = ea.isa()["output"]; //output
            repr[21] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            // nand
            repr[35] = ea.isa()["fixed_input"]; // input
            repr[36] = ea.isa()["nop_c"]; // nopc
            repr[37] = ea.isa()["fixed_input"]; // input
            repr[38] = ea.isa()["nand"]; // nand
            repr[39] = ea.isa()["output"]; //output
            repr[40] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            // ornot
            repr[50] = ea.isa()["fixed_input"]; // input
            repr[51] = ea.isa()["nop_c"]; // nopc
            repr[52] = ea.isa()["fixed_input"]; // input
            repr[53] = ea.isa()["nand"]; // nand
            repr[54] = ea.isa()["nand"]; // nand
            repr[55] = ea.isa()["output"]; //output
            repr[56] = ea.isa()["donate_res_to_group"]; // donate_res_to_group

            
            // or
            repr[65] = ea.isa()["fixed_input"]; // input
            repr[66] = ea.isa()["push"]; // push
            repr[67] = ea.isa()["nop_c"]; // nopc
            repr[68] = ea.isa()["pop"]; // pop
            repr[69] = ea.isa()["nop_a"]; // nopa
            repr[70] = ea.isa()["nand"]; // nand
            repr[71] = ea.isa()["fixed_input"]; // input
            repr[72] = ea.isa()["push"]; // push
            repr[73] = ea.isa()["nop_c"]; // nopc
            repr[74] = ea.isa()["pop"]; // pop
            repr[75] = ea.isa()["nand"]; // nand
            repr[76] = ea.isa()["nop_c"]; // nopc
            repr[77] = ea.isa()["swap"]; // swap
            repr[78] = ea.isa()["nand"]; // nand
            repr[79] = ea.isa()["output"]; //output
            repr[80] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            
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
