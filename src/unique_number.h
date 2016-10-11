#ifndef _UNIQUE_NUMBER_H_
#define _UNIQUE_NUMBER_H_

#include "stripes.h"

using namespace ealib;



LIBEA_MD_DECL(OPINION, "ea.stripes.ea.stripes.opinion", int); // unique number

DIGEVO_INSTRUCTION_DECL(get_opinion) {
    int rbx = hw.modifyRegister();
    hw.setRegValue(rbx, get<OPINION>(*p,0));
}

DIGEVO_INSTRUCTION_DECL(set_opinion) {
    int rbx = hw.modifyRegister();
    get<OPINION>(*p, 0) = hw.getRegValue(rbx);
}


DIGEVO_INSTRUCTION_DECL(inc_opinion) {
    get<OPINION>(*p, 0) += 1;
}

DIGEVO_INSTRUCTION_DECL(dec_opinion) {
    get<OPINION>(*p, 0) -= 1;
}


//! fitness.
struct permute_unique_number : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_unique_number(EA& ea) {
        double tmp_fit = eval_unique_number(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_unique_number(sea));
        put<STRIPE_FIT>(f,sea);
        
        return f;
    }
};


template <typename EA>
double eval_unique_number(EA& ea) {
    
    double f = 0.0;
    std::set<int> s;

    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
            typename EA::environment_type::location_type* l = &ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }
            
            s.insert(get<OPINION>(*(l->inhabitant()), 0));

            
            
        }
    }
    
    
    
    f = s.size();
    
    put<STRIPE_FIT>(f,ea);
    return f;
    
    
}


#endif
