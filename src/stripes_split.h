

#ifndef _EALIFE_STRIPES_SPLIT_H_
#define _EALIFE_STRIPES_SPLIT_H_

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


#include "stripes.h"

using namespace ealib;




//! Performs multicell replication using germ lines. The number of propagule cells determines the propagule size.
// One of these cells is selected, mutated, and then used to create the appropriate number of cells. Thus, the starting multicell
// offspring is clonal.
template <typename MEA>
struct stripes_split : end_of_update_event<MEA> {
    //! Constructor.
    stripes_split(MEA& mea) : end_of_update_event<MEA>(mea), _df("stripes_replication_evo_plane.dat") {
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
    virtual ~stripes_split() {
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
                
                
                double num_prop = ceil(get<PROP_COUNT>(*i,0) / 2.0);
                
                int prop_total_cost = get<PROPAGULE_BASE_COST>(*i) + (num_prop * get<PROPAGULE_PER_CELL_COST>(*i));
                
                
                int num_founder = i->traits().founder()->population().size();
                if (((get<MC_RESOURCE_UNITS>(*i) > prop_total_cost) && (num_prop > 0)) && ((i->population().size() > num_founder ) || (num_founder == get<POPULATION_SIZE>(*i)))){
                    
                    
                    // get a new subpopulation:
                    typename MEA::individual_ptr_type p = mea.make_individual();
                    p->initialize(mea.md());
                    p->reset_rng(mea.rng().seed());
                    
                    
                    // shuffle the population
                    std::random_shuffle(i->population().begin(), i->population().end(), i->rng());
                    
                    int num_moved = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((get<IS_PROPAGULE>(**j, 0) == 2) && (*j)->alive()) {
                            typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                            (*j)->genome().begin()+(*j)->hw().original_size());
                            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                            
                            inherits_from(**j, *q, *p);
                            
                            
                            const position_type pos = (*j)->position();
                            
                            
                            //iterator insert_at(iterator i, individual_ptr_type x, const position_type& pos) {
                            p->insert_at(p->end(), q, pos);
                            
                            (*j)->alive() = false;
                            i->events().death(**j,*i);
                            
                            
                            ++num_moved;
                            
                        }
                        if (num_moved == num_prop) {
                            break;
                        }
                    }
                    
                    
                    
                    
                    
                    offspring.insert(offspring.end(),p);
                    
                    int alive_count = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            alive_count++;
                        }

                    }
                    
                    // track last fitness AND time to rep
                    multicell_rep.push_back(get<MULTICELL_REP_TIME>(*i));
                    multicell_last_fitness.push_back(get<STRIPE_FIT>(*i));
                    multicell_size.push_back(alive_count);
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


#endif
