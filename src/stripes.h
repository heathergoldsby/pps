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


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int); // count the number of organisms that have the right color stripe

LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);
LIBEA_MD_DECL(PROPAGULE_IDENTICAL, "ea.stripes.prop_identical", int); // 0 NO - different, 1 YES - identical
LIBEA_MD_DECL(PROPAGULE_COPY, "ea.stripes.prop_copy", int); // 0 NO - move cells, 1 YES - copy cells
LIBEA_MD_DECL(PARENT_PATTERN, "ea.stripes.parent_pattern", std::string); // strip representation of parent pattern
LIBEA_MD_DECL(PATTERN, "ea.stripes.pattern", std::string); // strip representation of parent pattern
LIBEA_MD_DECL(PARENT_WEIGHT, "ea.stripes.parent_weight", double); // 0 - no  parent weight; 1 only sim to parent


//LIBEA_MD_DECL(IS_PROPAGULE, "ea.stripes.is_propagule", int); // 0 - no, 1 - cost paid, 2 - transfer


//LIBEA_MD_DECL(FIT_MAX, "ea.stripes.fit_max", int);
//LIBEA_MD_DECL(FIT_MIN, "ea.stripes.fit_min", int);
//LIBEA_MD_DECL(FIT_GAMMA, "ea.stripes.fit_gamma", double);
//LIBEA_MD_DECL(RES_UPDATE, "ea.stripes.res_update", int);
//LIBEA_MD_DECL(PROP_SIZE, "ea.stripes.propagule_size", double);
//LIBEA_MD_DECL(PROP_SIZE_OPTION, "ea.stripes.propagule_size_option", int);
//LIBEA_MD_DECL(PROP_SIZE_BOUND, "ea.stripes.propagule_size_bound", int);

LIBEA_MD_DECL(REP_COUNT, "ea.stripes.rep_count", int); // count the number of times a multicell has replicated. 


LIBEA_MD_DECL(PROP_COUNT, "ea.stripes.prop_count", int);
LIBEA_MD_DECL(AGING, "ea.stripes.aging", int); // 0 off; 1 on
LIBEA_MD_DECL(AGE, "ea.stripes.age", double);
LIBEA_MD_DECL(AGE_AMOUNT, "ea.stripes.age_amount", double);

LIBEA_MD_DECL(TRANS, "ea.stripes.trans", int); // 0 off; 1 on
LIBEA_MD_DECL(BASE, "ea.stripes.base", double);


//
//
///*!  transmit parent pattern to offspring.
// */
//template <typename EA>
//struct transmit_pattern_event : birth_event<EA> {
//    
//    //! Constructor.
//    transmit_pattern_event(EA& ea) : birth_event<EA>(ea) {
//    }
//    
//    //! Destructor.
//    virtual ~transmit_pattern_event() {
//    }
//    
//    /*! Called for every inheritance event.
//     */
//    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
//                            typename EA::individual_type& parent, // individual parent
//                            EA& ea) {
//        
//        
//        
//        // sweep through parent. create pattern string. place in offspring
//        std::string pattern = get<PATTERN>(parent, "");
//        put<PARENT_PATTERN>(pattern, offspring);
//        
//    }
//};




//! Stripe fitness.
struct permute_stripes : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_stripes(EA& ea) {
        double tmp_fit = eval_two_stripes(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_stripes(sea));
        double adj_f =f;
        put<STRIPE_FIT>(f,sea);
        
        // add in aging effects
        if (get<AGING>(sea,0) == 1) {
            double age_factor = get<AGE_AMOUNT>(sea);

            double age = get<AGE>(sea, 0);
            adj_f *= (1.0 - (age * age_factor) );
        }
        
        
        // add in fitness transformation
        if (get<TRANS>(sea,0) == 1) {
            double base = get<BASE>(sea);
            adj_f = pow(base, adj_f);
        }
        
        return adj_f;
    }
};

//! Stripe fitness.
struct permute_three_stripes : public fitness_function<unary_fitness<double> > {
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



//! Square fitness
struct permute_square : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_square(EA& ea) {
        double tmp_fit = eval_square(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_square(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};


//!
struct solid_control : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double evaluate_solid_control(EA& ea) {
        double tmp_fit = eval_solid_control(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(evaluate_solid_control(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};


//! two_color_control
struct two_color_control : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_two_color_control(EA& ea) {
        double tmp_fit = eval_two_color(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_two_color_control(sea));
        put<STRIPE_FIT>(f,sea);
    
        return f;
    }
};

//! Bullseye
struct permute_bullseye : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_bullseye(EA& ea) {
        double tmp_fit = eval_bullseye(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_bullseye(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};

//! Bullseye
struct permute_bullseye_fixed : public fitness_function<unary_fitness<double> > {
    template <typename EA>
    double eval_permute_bullseye_fixed(EA& ea) {
        double tmp_fit = eval_bullseye_fixed(ea);
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_bullseye_fixed(sea));
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
double eval_solid_control(EA& ea) {
    
    double num_correct = 0;
    
    
    
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
            
                        // NW
            if (lt == "NOT" ) { ++num_correct; }
        }
    }
        
    num_correct = num_correct * num_correct;
        
    
    put<STRIPE_FIT>(num_correct, ea);

    return num_correct;
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



template <typename EA>
double eval_two_color(EA& ea) {
    double num_not = 0;
    double num_nand = 0;
    
    
    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
            typename EA::environment_type::location_type* l = &ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
            
            if (lt == "nand") { ++num_nand; }
            if (lt == "not") { ++num_not; }
        }
    }
    
    
    
    
    double fit = num_not * num_nand;
    
    put<STRIPE_FIT>(fit,ea);
    return fit;
    
    
}



template <typename EA>
double eval_square(EA& ea) {
    
    double num_correct = 0;
        
    int num_square = 3;
    std::vector<double> not_fit(num_square);
    std::vector<double> nand_fit(num_square);
    std::vector<double> ornot_fit(num_square);
    std::vector<double> total_fit(6);
    
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
            
            
            if (x == 0 || x == 5 || y == 0 || y == 5) {
                if (lt == "not") {
                    not_fit[0] += 1;
                } else if (lt == "nand") {
                    nand_fit[0] += 1;
                } else if (lt == "ornot") {
                    ornot_fit[0] +=1;
                }
            } else if (x == 1 || x == 4 || y == 1 || y == 4) {
                if (lt == "not") {
                    not_fit[1] += 1;
                } else if (lt == "nand") {
                    nand_fit[1] += 1;
                } else if (lt == "ornot") {
                    ornot_fit[1] +=1;
                }
            } else {
                if (lt == "not") {
                    not_fit[2] += 1;
                } else if (lt == "nand") {
                    nand_fit[2] += 1;
                } else if (lt == "ornot") {
                    ornot_fit[2] +=1;
                }
            }
            
        }
    }
    
    total_fit[0] = (not_fit[0] + 1) * (nand_fit[1] + 1) * (ornot_fit[2] + 1);
    total_fit[1] = (not_fit[0] + 1) * (nand_fit[2] + 1) * (ornot_fit[1] + 1);
    total_fit[2] = (not_fit[1] + 1) * (nand_fit[0] + 1) * (ornot_fit[2] + 1);
    total_fit[3] = (not_fit[1] + 1) * (nand_fit[2] + 1) * (ornot_fit[0] + 1);
    total_fit[4] = (not_fit[2] + 1) * (nand_fit[0] + 1) * (ornot_fit[1] + 1);
    total_fit[5] = (not_fit[2] + 1) * (nand_fit[1] + 1) * (ornot_fit[0] + 1);

    
    double tmp_fit = std::max(total_fit[0], total_fit[1]);
    tmp_fit = std::max(tmp_fit, total_fit[2]);
    tmp_fit = std::max(tmp_fit, total_fit[3]);
    tmp_fit = std::max(tmp_fit, total_fit[4]);
    tmp_fit = std::max(tmp_fit, total_fit[5]);
    

    
    put<STRIPE_FIT>(tmp_fit, ea);
    return (tmp_fit);
    
}


template <typename EA>
double eval_bullseye(EA& ea) {
    
    int best_fit = 0;
    
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
            
            if (lt == "or") {
                // check the grid.
                int offset_x = 3 - x;
                int offset_y = 3 - y;
                std::vector<double> outer_ring(3);
                std::vector<double> middle_ring(3);
                std::vector<double> inner_ring(3);
                
                for (int check_x=0; check_x < max_x; ++check_x) {
                    for (int check_y=0; check_y< max_y; ++check_y){
                
                        int act_x = check_x + offset_x;
                        if (act_x >= max_x) {
                            act_x -= max_x;
                        }
                        if (act_x < 0) {
                            act_x += max_x;
                        }
                        
                        int act_y = check_y + offset_y;
                        if (act_y >= max_y) {
                            act_y -= max_y;
                        }
                        if (act_y < 0) {
                            act_y += max_y;
                        }
                        
                        typename EA::environment_type::location_type* l2 = &ea.env().location(act_x,act_y);
                        if (!l2->occupied()) {
                            continue;
                        }

                        std::string act_lt = get<LAST_TASK>(*(l2->inhabitant()),"");
                        if (act_lt == "") {
                            continue;
                        }
                        
                        
                        // the outer ring.
                        if (check_x == 0 || check_x == 6 || check_y == 0 || check_y == 6) {
                            if(act_lt == "not" ) {
                                outer_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                outer_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                outer_ring[2] +=1;
                            }
                            
                        } else if (check_x == 1 || check_x == 5 || check_y == 1 || check_y == 5) { // second ring
                            if(act_lt == "not" ) {
                                middle_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                middle_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                middle_ring[2] +=1;
                            }
                        } else if (check_x == 2 || check_x == 4 || check_y == 2 || check_y == 4) {
                            if(act_lt == "not" ) {
                                inner_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                inner_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                inner_ring[2] +=1;
                            }
                        }
                        
                    }
                }
                
                // compute fit of this particular position...
                std::vector<double> total_fit(6);
    
                total_fit[0] = (outer_ring[0] + 1) * (middle_ring[1] + 1) * (inner_ring[2] + 1);
                total_fit[1] = (outer_ring[0] + 1) * (middle_ring[2] + 1) * (inner_ring[1] + 1);
                total_fit[2] = (outer_ring[1] + 1) * (middle_ring[0] + 1) * (inner_ring[2] + 1);
                total_fit[3] = (outer_ring[1] + 1) * (middle_ring[2] + 1) * (inner_ring[0] + 1);
                total_fit[4] = (outer_ring[2] + 1) * (middle_ring[0] + 1) * (inner_ring[1] + 1);
                total_fit[5] = (outer_ring[2] + 1) * (middle_ring[1] + 1) * (inner_ring[0] + 1);
                
                
                double tmp_fit = std::max(total_fit[0], total_fit[1]);
                tmp_fit = std::max(tmp_fit, total_fit[2]);
                tmp_fit = std::max(tmp_fit, total_fit[3]);
                tmp_fit = std::max(tmp_fit, total_fit[4]);
                tmp_fit = std::max(tmp_fit, total_fit[5]);
                
                if (tmp_fit > best_fit) {
                    best_fit = tmp_fit;
                }
            }
        }
    }
    
    
    
    
    put<STRIPE_FIT>(best_fit, ea);
    return (best_fit);

    
}




template <typename EA>
double eval_bullseye_fixed(EA& ea) {
    
    int best_fit = 0;
    
    int max_x = get<SPATIAL_X>(ea);
    int max_y = get<SPATIAL_Y>(ea);
    
    std::string pattern;
    
    int x = 0;
    int y = 0;
    typename EA::environment_type::location_type* l = &ea.env().location(x,y);

    
    std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
        
                // check the grid.
                int offset_x = 3 - x;
                int offset_y = 3 - y;
                std::vector<double> outer_ring(3);
                std::vector<double> middle_ring(3);
                std::vector<double> inner_ring(3);
                
                for (int check_x=0; check_x < max_x; ++check_x) {
                    for (int check_y=0; check_y< max_y; ++check_y){
                        
                        
                        int act_x = check_x + offset_x;
                        if (act_x >= max_x) {
                            act_x -= max_x;
                        }
                        if (act_x < 0) {
                            act_x += max_x;
                        }
                        
                        int act_y = check_y + offset_y;
                        if (act_y >= max_y) {
                            act_y -= max_y;
                        }
                        if (act_y < 0) {
                            act_y += max_y;
                        }
                        
                        typename EA::environment_type::location_type* l2 = &ea.env().location(act_x,act_y);
                        if (!l2->occupied()) {
                            // add on to pattern.
                            pattern += "9";
                            continue;
                        }
                        
                        
                        std::string act_lt = get<LAST_TASK>(*(l2->inhabitant()),"");
                        if (act_lt == "") {
                            pattern += "0";
                            continue;
                        }
                        
                        if(lt == "not") {
                            pattern += "1";
                        }
                        if (lt == "nand") {
                            pattern += "2";
                        }
                        if (lt == "or") {
                            pattern += "3";
                        }
                        if (lt == "ornot") {
                            pattern += "4";
                        }
                       
                        
                        // the outer ring.
                        if (check_x == 0 || check_x == 6 || check_y == 0 || check_y == 6) {
                            if(act_lt == "not" ) {
                                outer_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                outer_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                outer_ring[2] +=1;
                            }
                            
                        } else if (check_x == 1 || check_x == 5 || check_y == 1 || check_y == 5) { // second ring
                            if(act_lt == "not" ) {
                                middle_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                middle_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                middle_ring[2] +=1;
                            }
                        } else if (check_x == 2 || check_x == 4 || check_y == 2 || check_y == 4) {
                            if(act_lt == "not" ) {
                                inner_ring[0] +=1;
                            } else if (act_lt == "nand") {
                                inner_ring[1] +=1;
                            } else if (act_lt == "ornot") {
                                inner_ring[2] +=1;
                            }
                        }
                        
                    }
                }
        
                // compute fit of this particular position...
                std::vector<double> total_fit(6);
                
                total_fit[0] = (outer_ring[0] + 1) * (middle_ring[1] + 1) * (inner_ring[2] + 1);
                total_fit[1] = (outer_ring[0] + 1) * (middle_ring[2] + 1) * (inner_ring[1] + 1);
                total_fit[2] = (outer_ring[1] + 1) * (middle_ring[0] + 1) * (inner_ring[2] + 1);
                total_fit[3] = (outer_ring[1] + 1) * (middle_ring[2] + 1) * (inner_ring[0] + 1);
                total_fit[4] = (outer_ring[2] + 1) * (middle_ring[0] + 1) * (inner_ring[1] + 1);
                total_fit[5] = (outer_ring[2] + 1) * (middle_ring[1] + 1) * (inner_ring[0] + 1);
                
                
                double tmp_fit = std::max(total_fit[0], total_fit[1]);
                tmp_fit = std::max(tmp_fit, total_fit[2]);
                tmp_fit = std::max(tmp_fit, total_fit[3]);
                tmp_fit = std::max(tmp_fit, total_fit[4]);
                tmp_fit = std::max(tmp_fit, total_fit[5]);
                
                if (tmp_fit > best_fit) {
                    best_fit = tmp_fit;
                }
    
    // add on if the middle of the bullseye is correct
    if (lt == "or")
    {
        best_fit *= 2;
    }
    
    
    double parent_weight = get<PARENT_WEIGHT>(ea,0);
    double sim_val = 0;
    
    if (parent_weight) {
        std::string p_string = get<PARENT_PATTERN>(ea,"");
        // compute similarity
        
        if (p_string != "") {
        
            for (int i = 0; i < p_string.length(); ++i) {
                if (p_string[i] == pattern[i]) {
                    ++sim_val;
                }
            }
        
            sim_val /= p_string.length();
        
            best_fit = best_fit + (sim_val * parent_weight);
        }

    }
    
    
    put<STRIPE_FIT>(best_fit, ea);
    put<PATTERN>(pattern, ea);
    return (best_fit);
    
    
}





#endif
