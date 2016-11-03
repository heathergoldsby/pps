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

//! Send a message with the cell's opinion to the currently-faced neighbor.
DIGEVO_INSTRUCTION_DECL(tx_op) {
    typename EA::environment_type::location_type& l=*ea.env().neighbor(p);
    if(l.occupied()) {
        l.inhabitant()->hw().deposit_message(get<OPINION>(*p,0), get<OPINION>(*p,0));
    }
}

//! Broadcast a message.
DIGEVO_INSTRUCTION_DECL(bc_op) {
    int rbx = hw.modifyRegister();
    int rcx = hw.nextRegister(rbx);
    
    typedef typename EA::environment_type::neighborhood_iterator neighborhood_iterator;
    std::pair<neighborhood_iterator,neighborhood_iterator> ni=ea.env().neighborhood(*p);
    for(; ni.first!=ni.second; ++ni.first) {
        typename EA::environment_type::location_type& l=*ni.first;
        if(l.occupied()) {
            l.inhabitant()->hw().deposit_message(get<OPINION>(*p,0), get<OPINION>(*p,0));
        }
    }
}



/*! epigenetic opinions
 */
template <typename EA>
struct epi_op_birth_event : birth_event<EA> {
    
    //! Constructor.
    epi_op_birth_event(EA& ea) : birth_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~epi_op_birth_event() {
    }
    
    /*! Called for every inheritance event.
     */
    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
                            typename EA::individual_type& parent, // individual parent
                            EA& ea) {
        //ea.env().face_org(parent, offspring);
        
        get<OPINION>(offspring, 0) = get<OPINION>(parent,0);
        
        
    }
};





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
            
            int op = get<OPINION>(*(l->inhabitant()), 0);
            
            if (op >= 0 && op <= get<POPULATION_SIZE>(ea)) {
                s.insert(op);
            }
            
        }
    }
    
    
    
    f = s.size();
    
    put<STRIPE_FIT>(f,ea);
    return f;
    
    
}


#endif
