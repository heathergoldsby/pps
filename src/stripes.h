//
//  stripes.h
//  ealife
//
//  Created by Heather Goldsby on 11/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//


#ifndef _EALIFE_STRIPES_H_
#define _EALIFE_STRIPES_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


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


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int); // count the number of organisms that have the right color stripe

LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);

LIBEA_MD_DECL(IS_PROPAGULE, "ea.stripes.is_propagule", int); // 0 - no, 1 - cost paid, 2 - transfer


LIBEA_MD_DECL(FIT_MAX, "ea.stripes.fit_max", int);
LIBEA_MD_DECL(FIT_MIN, "ea.stripes.fit_min", int);
LIBEA_MD_DECL(FIT_GAMMA, "ea.stripes.fit_gamma", double);
LIBEA_MD_DECL(RES_UPDATE, "ea.stripes.res_update", int);
LIBEA_MD_DECL(PROP_SIZE, "ea.stripes.propagule_size", double);
LIBEA_MD_DECL(PROP_SIZE_OPTION, "ea.stripes.propagule_size_option", int);
LIBEA_MD_DECL(PROP_SIZE_BOUND, "ea.stripes.propagule_size_bound", int);

LIBEA_MD_DECL(REP_COUNT, "ea.stripes.rep_count", int); // count the number of times a multicell has replicated. 


LIBEA_MD_DECL(PROP_COUNT, "ea.stripes.prop_count", int);

//! create_propagule - creates the propagule cell


//! deploys the propagule to the holding tank
DIGEVO_INSTRUCTION_DECL(deploy_propagule) {
    
    if (get<IS_PROPAGULE>(*p, 0) == 1) {
        put<IS_PROPAGULE>(2, *p);
        get<PROP_COUNT>(ea,0) += 1;
    }
}


//! Execute the next instruction if the multicell does not have a propagule cell
DIGEVO_INSTRUCTION_DECL(if_prop_cell_absent) {
    if(get<PROP_COUNT>(ea,0)) {
        hw.advanceHead(Hardware::IP);
    }
}




//! Stripe fitness.
struct permute_stripes : public fitness_function<unary_fitness<double>, nonstationaryS > {
    template <typename EA>
    double eval_permute_stripes(EA& ea) {
        double tmp_fit = eval_two_stripes(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_stripes(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};

//! Stripe fitness.
struct permute_three_stripes : public fitness_function<unary_fitness<double>, nonstationaryS > {
    template <typename EA>
    double eval_permute_three_stripes(EA& ea) {
        double tmp_fit = eval_three_stripes(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_three_stripes(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};





template <typename EA>
std::string get_last_task(int x, int y, EA& ea) {
    typename EA::environment_type::location_type* l = &ea.env().location(x,y);
    std::string lt = "";
    if (l->occupied()) {
        lt = get<LAST_TASK>(*(l->inhabitant()),"");
    }
    return lt;
}

template <typename EA>
double eval_three_stripes(EA& ea) {
    
    double num_correct = 0;
    
    
    //    accumulator_set<double, stats<tag::mean, tag::max> > sfit;
    int num_neighbors = 4;
    std::vector<double> not_fit(num_neighbors);
    std::vector<double> nand_fit(num_neighbors);
    std::vector<double> ornot_fit(num_neighbors);
    std::vector<double> total_fit(num_neighbors);
    
    int max_x = get<SPATIAL_X>(ea);
    int max_y = get<SPATIAL_Y>(ea);
    
    for (int x=0; x < max_x; ++x) {
        for (int y=0; y< max_y; ++y){
            typename EA::environment_type::location_type* l = &ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
            if (lt == "" ) { continue; }
            
            std::vector<double> temp_t_fit(num_neighbors);
            
            // Get the relevant neighbors...
            // We only check NW (0) , N (1) , NE (2) , E (3) (the rest are covered as the grid moves)
            int n = y - 1;
            if (n < 0) {
                n += max_y;
            }
            int w = x - 1;
            if (w < 0) { w += max_x; }
            int e = x + 1;
            if (e >= max_x) {
                e -= max_x;
            }
            
            
            
            // NW
            if (lt == get_last_task(w, n, ea)) { ++temp_t_fit[0]; }
            
            // N
            if (lt == get_last_task(x, n, ea)) { ++temp_t_fit[1]; }
            
            // NE
            if (lt == get_last_task(e, n, ea)) { ++temp_t_fit[2]; }
            
            // E
            if (lt == get_last_task(e, y, ea)) { ++temp_t_fit[3]; }
            
            if (lt == "not") {
                not_fit[0] += temp_t_fit[0];
                not_fit[1] += temp_t_fit[1];
                not_fit[2] += temp_t_fit[2];
                not_fit[3] += temp_t_fit[3];
            } else if (lt == "nand") {
                nand_fit[0] += temp_t_fit[0];
                nand_fit[1] += temp_t_fit[1];
                nand_fit[2] += temp_t_fit[2];
                nand_fit[3] += temp_t_fit[3];
            } else if (lt == "ornot") {
                ornot_fit[0] += temp_t_fit[0];
                ornot_fit[1] += temp_t_fit[1];
                ornot_fit[2] += temp_t_fit[2];
                ornot_fit[3] += temp_t_fit[3];
            }
        }
    }
    
    total_fit[0] = (not_fit[0] + 1) * (nand_fit[0] + 1) * (ornot_fit[0] + 1);
    total_fit[1] = (not_fit[1] + 1) * (nand_fit[1] + 1) * (ornot_fit[1] + 1);
    total_fit[2] = (not_fit[2] + 1) * (nand_fit[2] + 1) * (ornot_fit[2] + 1);
    total_fit[3] = (not_fit[3] + 1) * (nand_fit[3] + 1) * (ornot_fit[3] + 1);
    
    double min_fit = 1;
    double max_fit = pow((get<POPULATION_SIZE>(ea) / 3), 3);
    
    double tmp_fit = std::max(total_fit[0], total_fit[1]);
    tmp_fit = std::max(tmp_fit, total_fit[2]);
    tmp_fit = std::max(tmp_fit, total_fit[3]);
    
    if (tmp_fit < min_fit) {
        tmp_fit = min_fit;
    }
    
    
    put<STRIPE_FIT>(tmp_fit, ea);
    return (tmp_fit);
    //    // rescale fitness!
    //    double rescaled_fit = (get<FIT_MAX>(ea) - get<FIT_MIN>(ea)) * pow (((tmp_fit - min_fit) / (max_fit - min_fit)), (get<FIT_GAMMA>(ea))) + get<FIT_MIN>(ea);
    //
    //
    //    put<STRIPE_FIT>(rescaled_fit,ea);
    
}




template <typename EA>
double eval_two_stripes(EA& ea) {
    // vert stripes
    double one_fit_not = 0;
    double one_fit_nand = 0;
    double two_fit_not = 0;
    double two_fit_nand = 0;
    // horizontal stripes
    double three_fit_not = 0;
    double three_fit_nand = 0;
    double four_fit_not = 0;
    double four_fit_nand = 0;
    // diagonal stripes
    double five_fit_not = 0;
    double five_fit_nand = 0;
    double six_fit_not = 0;
    double six_fit_nand = 0;
    
    
    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
            typename EA::environment_type::location_type* l = &ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
            // Vertical stripes!
            if((y % 2) == 0) {
                if (lt == "nand") { ++one_fit_nand; }
                if (lt == "not") { ++two_fit_not; }
            } else {
                if(lt == "not") { ++one_fit_not; }
                if (lt == "nand") { ++two_fit_nand; }
            }
            
            // Horizontal stripes
            if ((x % 2) == 0) {
                if (lt == "nand") { ++three_fit_nand; }
                if (lt == "not") { ++four_fit_not; }
            } else {
                if(lt == "not") { ++three_fit_not; }
                if (lt == "nand") { ++four_fit_nand; }
            }
            
            // Diagonal stripes
            if ((x % 2) == (y % 2)) {
                if(lt == "not") { ++five_fit_not; }
                if (lt == "nand") { ++six_fit_nand; }
            } else {
                if(lt == "nand") { ++five_fit_nand; }
                if (lt == "not") { ++six_fit_not; }
            }
        }
    }
    

    
    double tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
    double tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
    double tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
    double tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
    double tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
    double tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
    double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
    tmp_fit = std::max(tmp_fit, tmp_three_fit);
    tmp_fit = std::max(tmp_fit, tmp_four_fit);
    tmp_fit = std::max(tmp_fit, tmp_five_fit);
    tmp_fit = std::max(tmp_fit, tmp_six_fit);
    
    put<STRIPE_FIT>(tmp_fit,ea);
    return tmp_fit;


}


#endif
