

#ifndef _EALIFE_N_NEIGHBORS_H_
#define _EALIFE_N_NEIGHBORS_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <cmath>

#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/subpopulation_founder.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;

LIBEA_MD_DECL(N_NEIGHBORS, "ea.stripes.n_neighbors", int); // get desired number of neigbhors

//! Stripe fitness.
struct permute_n_neighbors : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_n_neighbors(EA& ea) {
        double tmp_fit = eval_n_neighbors(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_n_neighbors(sea));
        put<STRIPE_FIT>(f,sea);
        
        return f;
    }
};


template <typename EA>
double eval_n_neighbors(EA& ea) {
    
    double f = 0.0;
    int n_neighbors = get<N_NEIGHBORS>(ea);

    
    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
            typename EA::environment_type::location_type* l = &ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }

            int neighbors = 0;
            int ny = y - 1;
            if (ny < 0) {
                ny = get<SPATIAL_Y>(ea) - 1;
            }
            int sy = y + 1;
            if (sy > get<SPATIAL_Y>(ea)) {
                sy = 0;
            }
            int ex = x + 1;
            if (sy > get<SPATIAL_X>(ea)) {
                ex = 0;
            }
            int wx = x - 1;
            if (wx < 0) {
                wx = get<SPATIAL_X>(ea) - 1;
            }
            
            typename EA::environment_type::location_type* n;

            // north row

            n = &ea.env().location(ex, ny);
            if (n->occupied()) { ++neighbors; }
            
            
                
            n = &ea.env().location(x, ny);
            if (n->occupied()) { ++neighbors; }
            
            n = &ea.env().location(wx, ny);
            if (n->occupied()) { ++neighbors; }
            
            // middle row
            n = &ea.env().location(ex, y);
            if (n->occupied()) { ++neighbors; }
            
            n = &ea.env().location(wx, y);
            if (n->occupied()) { ++neighbors; }
            
            // south row
            n = &ea.env().location(ex, sy);
            if (n->occupied()) { ++neighbors; }
            
            n = &ea.env().location(x, sy);
            if (n->occupied()) { ++neighbors; }
            
            n = &ea.env().location(wx, sy);
            if (n->occupied()) { ++neighbors; }
            
            if (n_neighbors == neighbors) {
                f+=1.0;
            } else {
            
                f += (1.0 / (fabs(n_neighbors - neighbors) * n_neighbors));
            }
            
        }
    }
    
    
    

    
    put<STRIPE_FIT>(f,ea);
    return f;
    
    
}

#endif
