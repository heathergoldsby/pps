//
//  evo_propagule_ins.h
//  ps
//
//  Created by Heather Goldsby on 4/22/16.
//  Copyright © 2016 Michigan State University. All rights reserved.
//
#ifndef _EALIFE_EVO_PROPAGULE_INS_H_
#define _EALIFE_EVO_PROPAGULE_INS_H_

#include "stripes.h"

using namespace ealib;


LIBEA_MD_DECL(IS_PROPAGULE, "ea.stripes.is_propagule", int); // 0 - no, 1 - prop
LIBEA_MD_DECL(DEPLOY_ONE, "ea.stripes.deploy_one", int);

DIGEVO_INSTRUCTION_DECL(deploy_propagule) {
    
    // If the cell only ever wants to deploy one propagule... use this
    if ((get<DEPLOY_ONE>(*p,0)) == 1) {
        if(get<PROP_COUNT>(ea,0) > 0) {
            return;
        }
    }
    
    if (get<IS_PROPAGULE>(*p,0) != 1) {
        put<IS_PROPAGULE>(1, *p);
        get<PROP_COUNT>(ea,0) += 1;
    }
}


DIGEVO_INSTRUCTION_DECL(deploy_one_propagule) {
    get<DEPLOY_ONE>(*p,0) = 1;
}

//! Execute the next instruction if the multicell does not have a propagule cell
DIGEVO_INSTRUCTION_DECL(if_prop_cell_absent) {
    if(get<PROP_COUNT>(ea,0)) {
        hw.advanceHead(Hardware::IP);
    }
}


//! Execute the next instruction if the multicell does not have a propagule cell
DIGEVO_INSTRUCTION_DECL(only_deploy_one) {
    if(get<PROP_COUNT>(ea,0)) {
        put<IS_PROPAGULE>(1, *p);
        get<PROP_COUNT>(ea,0) += 1;
    }
}

/* Questions related to fecundity: 
 - if you pick one... what happens if you make many offspring... 
 - what about if you pick a larger number? 
 - if you pick one... and then don't have any, then what?
 - what if you want propagule of size 1, but MANY of them.
 - relationship between size and number of offspring
 */


#endif /* evo_propagule_ins_h */
