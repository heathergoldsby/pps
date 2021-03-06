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
LIBEA_MD_DECL(STRIPE_FIT_PT, "ea.stripes.fit.pt", int);
LIBEA_MD_DECL(ANCESTOR, "ea.stripes.ancestor", int);
LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);

LIBEA_MD_DECL(IS_PROPAGULE, "ea.stripes.is_propagule", int); // 0 - no, 1 - cost paid, 2 - transfer
LIBEA_MD_DECL(PROPAGULE_COST, "ea.stripes.propagule_cost", int);
LIBEA_MD_DECL(PROPAGULE_BASE_COST, "ea.stripes.propagule_base_cost", int);
LIBEA_MD_DECL(PROPAGULE_PER_CELL_COST, "ea.stripes.propagule_per_cell_cost", int);
LIBEA_MD_DECL(PROPAGULE_FAIL_PROB, "ea.stripes.propagule_fail_prob", double);

LIBEA_MD_DECL(GERM_STATUS, "ea.stripes.germ_status", int);


LIBEA_MD_DECL(STRIPE_FIT_FUNC, "ea.stripes.fit_func", int); // 0 = sum; 1 = prod
LIBEA_MD_DECL(MULTICELL_REP_TIME, "ea.stripes.mc_rep_time", int);
LIBEA_MD_DECL(MC_RESOURCE_UNITS, "ea.stripes.mc_res_units", int);
LIBEA_MD_DECL(FIT_MAX, "ea.stripes.fit_max", int);
LIBEA_MD_DECL(FIT_MIN, "ea.stripes.fit_min", int);
LIBEA_MD_DECL(FIT_GAMMA, "ea.stripes.fit_gamma", double);
LIBEA_MD_DECL(RES_UPDATE, "ea.stripes.res_update", int);
LIBEA_MD_DECL(PROP_SIZE, "ea.stripes.propagule_size", double);
LIBEA_MD_DECL(PROP_SIZE_OPTION, "ea.stripes.propagule_size_option", int);
LIBEA_MD_DECL(PROP_SIZE_BOUND, "ea.stripes.propagule_size_bound", int);

LIBEA_MD_DECL(DEATH_PROB, "ea.stripes.death_prob", double);
LIBEA_MD_DECL(SL_PERIOD, "ea.stripes.swap_location_period", int);
LIBEA_MD_DECL(NUM_SWAPS, "ea.stripes.num_swaps", int);

LIBEA_MD_DECL(PROP_COUNT, "ea.stripes.prop_count", int);
LIBEA_MD_DECL(JUV_PERIOD, "ea.stripes.juv_period", int);
LIBEA_MD_DECL(DEPLOY_ONE, "ea.stripes.deploy_one", int);
LIBEA_MD_DECL(CELL_COUNT_AT_SEQUESTER, "ea.stripes.cell_count_at_sequester", int);


LIBEA_MD_DECL(PARENT, "ea.stripes.parent", int); // 0 - no offspring, ++ offspring.


////! create_propagule - creates the propagule cell
//DIGEVO_INSTRUCTION_DECL(create_propagule) {
//    p->hw().add_cost(get<PROPAGULE_COST>(ea));
//    
//    if (ea.rng().p() > get<PROPAGULE_FAIL_PROB>(ea)) {
//        put<IS_PROPAGULE>(1, *p);
//    }
//    
//}
//
//! deploys the propagule to the holding tank
DIGEVO_INSTRUCTION_DECL(deploy_propagule) {
    
    // If the cell only ever wants to deploy one propagule... use this
//    if (get<DEPLOY_ONE>(ea,0) == 1) {
    if (get<DEPLOY_ONE>(*p,0) == 1) {
        if(get<PROP_COUNT>(ea,0) > 0) {
            return;
        }
    }
    
    if (get<IS_PROPAGULE>(*p,0) != 2) {
        put<IS_PROPAGULE>(2, *p);
        get<PROP_COUNT>(ea,0) += 1;
    }
}



DIGEVO_INSTRUCTION_DECL(deploy_one_propagule) {
//    get<DEPLOY_ONE>(ea,0) = 1;
    get<DEPLOY_ONE>(*p,0) = 1;
}

//! Increment the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(inc_propagule_size){
    get<PROP_SIZE>(*p, get<NUM_PROPAGULE_CELL>(ea))++;
}

//! Decrement the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(dec_propagule_size){
    double p1 = get<PROP_SIZE>(*p, get<NUM_PROPAGULE_CELL>(ea));
    p1 -= 1.0;
}

//! Get the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(get_propagule_size){
    hw.setRegValue(hw.modifyRegister(), get<PROP_SIZE>(*p, get<NUM_PROPAGULE_CELL>(ea)));
}

DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
}


//! Execute the next instruction if the multicell does not have a propagule cell
DIGEVO_INSTRUCTION_DECL(if_prop_cell_absent) {
    if(get<PROP_COUNT>(ea,0)) {
        hw.advanceHead(Hardware::IP);
    }
}



template <typename EA>
struct ps_birth_event : birth_event<EA> {
    
    //! Constructor.
    ps_birth_event(EA& ea) : birth_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~ps_birth_event() {
    }
    
    /*! Called for every inheritance event.
     */
    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
                            typename EA::individual_type& parent, // individual parent
                            EA& ea) {
        
        get<PARENT>(parent,0) ++;
        
    }
};




//! Stripe fitness.
struct permute_stripes : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    int eval_permute_stripes(EA& ea) {
        
//        
//        // sterile
//        if(get<PROP_COUNT>(ea,0) == 0) {
//            return 1.0;
//        }

        
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
                typename EA::location_type& l = ea.env().location(x,y);
                if (!l.occupied()) {
                    continue;
                }
                
                std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
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
        
        //return tmp_fit;
        return tmp_one_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_stripes(sea));
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

//! Stripe fitness.
struct permute_three_stripes : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    double eval_permute_three_stripes(EA& ea) {
        
//        // sterile
//        if(get<PROP_COUNT>(ea,0) == 0) {
//            return 1.0;
//        }
        
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
        
        //double f = pow(1.05, tmp_fit);
        double f = tmp_fit;
        
        
        return f;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_three_stripes(sea));
        put<STRIPE_FIT>(f,sea);
        return f;
    }
};


        
        

template <typename EA>
void eval_permute_three_stripes(EA& ea) {
    
    
//    // sterile
//    if(get<PROP_COUNT>(ea,0) == 0) {
//        put<STRIPE_FIT>(0.0, ea);
//        return;
//    }
    
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
    
    //double f = pow(1.05, tmp_fit);
    double f = tmp_fit;
    
    put<STRIPE_FIT>(f, ea);
    
    

    
}



//! Kill a subset of the population....
template <typename EA>
struct swap_locations : end_of_update_event<EA> {
    //! Constructor.
    swap_locations(EA& ea) : end_of_update_event<EA>(ea), _df("swap.dat") {}
    
    
    //! Destructor.
    virtual ~swap_locations() {
    }
    
    virtual void operator()(EA& ea) {
        
        int sl = get<SL_PERIOD>(ea,1);
        if ((ea.current_update() % sl) == 1) {
            
            // go through each multicell...
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                for (int j=0; j<get<NUM_SWAPS>(ea,0); ++j) {
                    
                    int s = get<POPULATION_SIZE>(*i);
                    std::size_t pos1 = i->rng()(s);
                    std::size_t pos2 = i->rng()(s);
            
                    i->env().swap_locations(pos1, pos2);
                }

            
            }
        
        }
    }
    
    datafile _df;
    
};



//! Kill a subset of the population....
template <typename EA>
struct random_death : end_of_update_event<EA> {
    //! Constructor.
    random_death(EA& ea) : end_of_update_event<EA>(ea), _df("random_death.dat") {
        _df.add_field("update")
        .add_field("pop_size")
        .add_field("num_killed");
    }
    
    
    //! Destructor.
    virtual ~random_death() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
                // go through each multicell...
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    for(typename EA::subpopulation_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        // get prob...
                        if ((*j)->alive()) {
                            if (ea.rng().p() < get<DEATH_PROB>(ea,0)) {
                                (*j)->alive() = false;
                                i->events().death(**j,*i);
                                
                            }
                        }
                    }
                    
                }
    
        }
    
    datafile _df;

};



//! Performs multicell replication using germ lines. The number of propagule cells determines the propagule size.
// One of these cells is selected, mutated, and then used to create the appropriate number of cells. Thus, the starting multicell
// offspring is clonal.
template <typename MEA>
struct stripes_ep3 : end_of_update_event<MEA> {
    //! Constructor.
    stripes_ep3(MEA& mea) : end_of_update_event<MEA>(mea), _df("stripes_replication_evo_plane.dat") {
        _df.add_field("update")
        .add_field("mean_rep_time")
        .add_field("last_fitness")
        .add_field("last_fitness_pt")
        .add_field("multicell_size")
        .add_field("prop_size")
        .add_field("replication_count");
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~stripes_ep3() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(MEA& mea) {
        
        int ru = get<RES_UPDATE>(mea,1);
        if ((mea.current_update() % ru) == 0) {
            
            // See if any subpops have exceeded the resource threshold
            typename MEA::population_type offspring;
            for(typename MEA::iterator i=mea.begin(); i!=mea.end(); ++i) {
                
                // Calc fitness for each subpop
                eval_permute_stripes(*i);
                
                // copy the stripe fit to the accumulator and also the subpop
                double sf =get<STRIPE_FIT>(*i);
                put<STRIPE_FIT>(sf, *i);
                get<MC_RESOURCE_UNITS>(*i,0) += sf;
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                
                // figure out which individuals from the parent comprise the propagule:
                typedef typename MEA::subpopulation_type::population_type propagule_type;
                propagule_type propagule;
                
                
                int num_prop = get<PROP_COUNT>(*i,0);
                
                int prop_total_cost = get<PROPAGULE_BASE_COST>(*i) + (num_prop * get<PROPAGULE_PER_CELL_COST>(*i));
                
                
                int num_founder = i->traits().founder()->population().size();
                if (((get<MC_RESOURCE_UNITS>(*i) > prop_total_cost) && (num_prop > 0)) && ((i->population().size() > num_founder ) || (num_founder == get<POPULATION_SIZE>(*i)))){
                    
                    
                    // get a new subpopulation:
                    typename MEA::individual_ptr_type p = mea.make_individual();
                    p->initialize(mea.md());
                    p->reset_rng(mea.rng().seed());
                    
                    
                    
                    // create a mutator
                    configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
                    
                    // shuffle the population
                    std::random_shuffle(i->population().begin(), i->population().end(), i->rng());
                    
                    
                    // find the first germ cell...
                    typename MEA::subpopulation_type::individual_ptr_type germ;

                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if (get<IS_PROPAGULE>(**j, 0) == 2) {
                            typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                            (*j)->genome().begin()+(*j)->hw().original_size());
                            germ = p->make_individual(r);
                            inherits_from(**j, *germ, *p);
                            mutate(*germ,m,*p);
                            // double check breaks out of loop
                            break;
                        }
                    }
                    
                    int s = get<POPULATION_SIZE>(mea);
                    std::vector<int> open_pos (s);
                    for( int n = 0 ; n < s ; ++n ) {
                        open_pos[ n ] = n;
                    }
                    
                    
                    std::size_t t = p->rng()(open_pos.size());
                    std::size_t pos = open_pos[t];
                    open_pos.erase(open_pos.begin() + t);
                    
                    p->insert_at(p->end(), germ, p->env().location(pos).position());
                    
                    for (int k=1; k<num_prop; ++k) {
                        // check if this is a valid way to copy an individual
                        typename MEA::subpopulation_type::individual_ptr_type o(germ);
                        
                        t = p->rng()(open_pos.size());
                        pos = open_pos[t];
                        open_pos.erase(open_pos.begin() + t);
                        p->insert_at(p->end(), o, p->env().location(pos).position());
                        
                    }
                    
                    
                    offspring.insert(offspring.end(),p);
                    
                    
                    // track last fitness AND time to rep
                    multicell_rep.push_back(get<MULTICELL_REP_TIME>(*i));
                    multicell_last_fitness.push_back(get<STRIPE_FIT>(*i));
                    multicell_size.push_back(i->size());
                    multicell_last_fitness_pt.push_back(get<STRIPE_FIT_PT>(*i));
                    multicell_prop_size.push_back(num_prop);
                    ++num_rep;
                    
                    // reset parent multicell
                    i->resources().reset();
                    put<MC_RESOURCE_UNITS>(0,*i);
                    put<MULTICELL_REP_TIME>(0,*i);
                    
                    // i == parent individual;
                    typename MEA::population_type parent_pop, offspring_pop;
                    parent_pop.push_back(*i.base());
                    offspring_pop.push_back(p);
                    inherits(parent_pop, offspring_pop, mea);
                    
                    if (multicell_rep.size() > 100) {
                        multicell_rep.pop_front();
                        multicell_last_fitness.pop_front();
                        multicell_last_fitness_pt.pop_front();
                        multicell_size.pop_front();
                        multicell_prop_size.pop_front();
                    }
                    
                }
            }
            
            
            // select surviving parent groups
            if (offspring.size() > 0) {
                int n = get<METAPOPULATION_SIZE>(mea) - offspring.size();
                
                typename MEA::population_type survivors;
                select_n<selection::random< > >(mea.population(), survivors, n, mea);
                
                // add the offspring to the list of survivors:
                survivors.insert(survivors.end(), offspring.begin(), offspring.end());
                
                // and swap 'em in for the current population:
                std::swap(mea.population(), survivors);
            }
        }
        
        
        if ((mea.current_update() % 100) == 0) {
            if (multicell_rep.size() > 0) {
                _df.write(mea.current_update())
                .write(std::accumulate(multicell_rep.begin(), multicell_rep.end(), 0.0)/multicell_rep.size())
                .write(std::accumulate(multicell_last_fitness.begin(), multicell_last_fitness.end(), 0.0)/multicell_last_fitness.size())
                .write(std::accumulate(multicell_last_fitness_pt.begin(), multicell_last_fitness_pt.end(), 0.0)/multicell_last_fitness_pt.size())
                .write(std::accumulate(multicell_size.begin(), multicell_size.end(), 0.0)/multicell_size.size())
                .write(std::accumulate(multicell_prop_size.begin(), multicell_prop_size.end(), 0.0)/multicell_prop_size.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
            } else {
                _df.write(mea.current_update())
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(num_rep)
                .endl();
            }
        }
        
    }
    
    datafile _df;
    std::deque<double> multicell_rep;
    std::deque<double> multicell_last_fitness;
    std::deque<double> multicell_last_fitness_pt;
    std::deque<double> multicell_size;
    std::deque<double> multicell_prop_size;
    
    int num_rep;
    
    
};


//template <typename EA>
//int int_get_prop_size(EA& ea) {
//    int ps = 0;
//    accumulator_set<double, stats<tag::median> > p_size;
//    
//    switch (get<PROP_SIZE_OPTION>(ea,0)) {
//        case 0: { // config
//            ps = get<PROP_SIZE_OPTION>(ea,0);
//            break;
//        }
//        case 1:{ // vote
//            for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
//                for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
//                    typename EA::environment_type::location_type* l = &ea.env().location(x,y);
//
//                    if (!l->occupied()) {
//                        continue;
//                    }
//                    
//                    p_size(get<PROP_SIZE>(*(l->inhabitant()),get<NUM_PROPAGULE_CELL>(ea)));
//                    
//                }
//            }
//            
//            ps = median(p_size);
//            if (ps > get<POPULATION_SIZE>(ea)) { ps = get<POPULATION_SIZE>(ea); }
//            break;
//        }
//        case 2: {// num germ
//            for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
//                for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
//                    typename EA::environment_type::location_type* l = &ea.env().location(x,y);
//                    if (!l->occupied()) {
//                        continue;
//                    }
//                    
//                    if (get<GERM_STATUS>(*(l->inhabitant()),1)) {
//                        ++ps;
//                    }
//                }
//            }
//            
//            break;
//        }
//    }
//    return ps;
//    
//}

template <typename EA>
void eval_permute_stripes(EA& ea) {
    
//    // sterile
//    if(get<PROP_COUNT>(ea,0) == 0) {
//            put<STRIPE_FIT>(0.0, ea);
//            return;
//    }
    
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
    
    // fix to use fit func for prod or sum
    double tmp_one_fit;
    double tmp_two_fit;
    double tmp_three_fit;
    double tmp_four_fit;
    double tmp_five_fit;
    double tmp_six_fit;
    double min_fit;
    double max_fit;
    
    switch (get<STRIPE_FIT_FUNC>(ea,1)) {
        case 0: {
            tmp_one_fit = one_fit_not +  one_fit_nand;
            tmp_two_fit = two_fit_not + two_fit_nand;
            tmp_three_fit = three_fit_not + three_fit_nand;
            tmp_four_fit = four_fit_not + four_fit_nand;
            tmp_five_fit = five_fit_not + five_fit_nand;
            tmp_six_fit = six_fit_not + six_fit_nand;
            min_fit = get<POPULATION_SIZE>(ea) / 2;
            max_fit = get<POPULATION_SIZE>(ea);
            break;
        }
        case 1:
        default:
            tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
            tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
            tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
            tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
            tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
            tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
            min_fit = 0;
            max_fit = pow(((get<POPULATION_SIZE>(ea) / 2) + 1), 2);
            break;
    }
    
    
    double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
    tmp_fit = std::max(tmp_fit, tmp_three_fit);
    tmp_fit = std::max(tmp_fit, tmp_four_fit);
    tmp_fit = std::max(tmp_fit, tmp_five_fit);
    tmp_fit = std::max(tmp_fit, tmp_six_fit);
    
    if (tmp_fit < min_fit) {
        tmp_fit = min_fit;
    }
    
    
    put<STRIPE_FIT_PT>(tmp_fit, ea);
    put<STRIPE_FIT>(tmp_fit,ea);

//    // rescale fitness!
//    double rescaled_fit = (get<FIT_MAX>(ea) - get<FIT_MIN>(ea)) * pow (((tmp_fit - min_fit) / (max_fit - min_fit)), (get<FIT_GAMMA>(ea))) + get<FIT_MIN>(ea);
//    
//    
//    put<STRIPE_FIT>(rescaled_fit,ea);

}

/*! If an organism is a propagule, it updates the number of propagule cells accordingly.
 */
template <typename EA>
struct prop_death_event : death_event<EA> {
    
    //! Constructor.
    prop_death_event(EA& ea) : death_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~prop_death_event() {
    }
    
    /*! Called for every inheritance event. We are using the orientation of the first parent...
     */
    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
                            EA& ea) {
        if (get<IS_PROPAGULE>(offspring, 0) == 2) {
            get<PROP_COUNT>(ea,0) -= 1;
        }
        
    }
};

#endif
