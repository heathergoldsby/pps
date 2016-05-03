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
LIBEA_MD_DECL(PROPAGULE_SIZE, "ea.stripes.prop_size", int); // size of proposed propagule.
LIBEA_MD_DECL(START_PROPAGULE_SIZE, "ea.stripes.start_prop_size", int); // start size of proposed propagule.


DIGEVO_INSTRUCTION_DECL(deploy_propagule) {
    
    if (get<IS_PROPAGULE>(*p,0) != 1) {
        put<IS_PROPAGULE>(1, *p);
        get<PROP_COUNT>(ea,0) += 1;
    }
}


DIGEVO_INSTRUCTION_DECL(prop_size_1) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 1;
}
DIGEVO_INSTRUCTION_DECL(prop_size_2) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 2;
}
DIGEVO_INSTRUCTION_DECL(prop_size_3) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 3;
}
DIGEVO_INSTRUCTION_DECL(prop_size_4) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 4;
}
DIGEVO_INSTRUCTION_DECL(prop_size_5) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 5;
}
DIGEVO_INSTRUCTION_DECL(prop_size_6) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 6;
}
DIGEVO_INSTRUCTION_DECL(prop_size_7) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 7;
}
DIGEVO_INSTRUCTION_DECL(prop_size_8) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 8;
}
DIGEVO_INSTRUCTION_DECL(prop_size_9) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 9;
}
DIGEVO_INSTRUCTION_DECL(prop_size_10) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 10;
}
DIGEVO_INSTRUCTION_DECL(prop_size_11) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 11;
}
DIGEVO_INSTRUCTION_DECL(prop_size_12) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 12;
}
DIGEVO_INSTRUCTION_DECL(prop_size_13) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 13;
}
DIGEVO_INSTRUCTION_DECL(prop_size_14) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 14;
}
DIGEVO_INSTRUCTION_DECL(prop_size_15) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 15;
}
DIGEVO_INSTRUCTION_DECL(prop_size_16) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 16;
}
DIGEVO_INSTRUCTION_DECL(prop_size_17) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 17;
}
DIGEVO_INSTRUCTION_DECL(prop_size_18) {
    get<PROPAGULE_SIZE>(*p,get<START_PROPAGULE_SIZE>(*p)) = 18;
}

#endif /* evo_propagule_ins_h */
