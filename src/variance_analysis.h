#ifndef _PS_VAR_ANALYSIS_H_
#define _PS_VAR_ANALYSIS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
//#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "stripes.h"
namespace ealib {
    namespace analysis {
        
        LIBEA_ANALYSIS_TOOL(variance_analysis) {
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            int max_fit = 0;
            typename EA::individual_type best;
            accumulator_set<double, stats<tag::min, tag::mean, tag::max, tag::variance> > fit;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                // not preserving location
                typename EA::individual_ptr_type p = ea.make_individual();
                p->initialize(ea.md());
                //p->reset_rn(ea.rng().seed());
                
                //typename EA::individual_type tmp(*i->traits().founder());
                typedef typename EA::subpopulation_type::population_type propagule_type;
                for (typename propagule_type::iterator j=(*i->traits().founder()).population().begin(); j!= (*i->traits().founder()).population().end(); ++j) {
                    int pos = 0;
                    typename EA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                   (*j)->genome().begin()+(*j)->hw().original_size());
                    typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                    inherits_from(**j, *q, *p);
                    p->insert_at(p->end(),q, p->env().location((*j)->position()).position());
                }
                
                for (int j=0; j<=update_max; ++j) {
                    p->update();
                }
                recalculate_fitness(*p, ea);
                double sf =get<STRIPE_FIT>(*p);
                
                if (sf > max_fit) {
                    max_fit = sf;
                    best = *i;
                }
            }
            
            
            datafile df("dominant_fitness_variance.dat");
            for (int q=0; q < 1000; q++) {
                
                int pop_size = get<POPULATION_SIZE>(ea);
                
                std::vector<int> used_pos;
                std::vector<int> used_pos_with_avail_neighbors;
                
                //typename EA::individual_type best_founder = *best.traits().founder();
                // not preserving location
                typename EA::individual_ptr_type p = ea.make_individual();
                p->initialize(ea.md());
                
                // use different seeds
                p->rng().reset(rand());
                ea.rng().resent(rand()); 
                
                int pos = ea.rng().uniform_integer(0,pop_size);
                used_pos.push_back(pos);
                used_pos_with_avail_neighbors.push_back(pos);
                
                
                typedef typename EA::subpopulation_type::population_type propagule_type;
                for (typename propagule_type::iterator j=(*best.traits().founder()).population().begin(); j!= (*best.traits().founder()).population().end(); ++j) {
                    int pos = 0;
                    typename EA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                   (*j)->genome().begin()+(*j)->hw().original_size());
                    typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                    
                    inherits_from(**j, *q, *p);
                    
                    
                    int max_x =get<SPATIAL_X>(ea);
                    
                    bool not_placed = true;
                    int place = -1;
                    
                    while (not_placed) {
                        int n = ea.rng().uniform_integer(0,used_pos_with_avail_neighbors.size());
                        int neighbor = used_pos_with_avail_neighbors[n];
                        int dir_try = 0;
                        int dir = ea.rng().uniform_integer(0,3);
                        
                        
                        
                        while (dir_try <= 4) {
                            
                            dir = (dir + 1);
                            if (dir > 3) { dir = 0; }
                            dir_try++;
                            
                            // N
                            if (dir == 0) {
                                if ((neighbor - max_x) < 0) { continue; }
                                place = neighbor - max_x;
                            } else if (dir == 1) { // E
                                if ((neighbor % max_x) == (max_x -1)) { continue; }
                                place = neighbor + 1;
                            } else if (dir == 2) { // S
                                if ((neighbor + max_x) >= pop_size) { continue; }
                                place = neighbor + max_x;
                            } else if (dir == 3) { // W
                                if ((neighbor % max_x) == 0) { continue; }
                                place = neighbor - 1;
                            }
                            
                            // already used
                            if (std::find(used_pos.begin(), used_pos.end(), place) != used_pos.end()) {
                                place = -1;
                                continue;
                            } else {
                                break;
                            }
                            
                        }
                        
                        // not placed, then try again
                        if (place == -1) {
                            used_pos_with_avail_neighbors.erase(used_pos_with_avail_neighbors.begin()+n);
                            continue;
                        }
                        
                        
                        
                        p->insert_at(p->end(),q, p->env().location(place).position());
                        used_pos.push_back(place);
                        used_pos_with_avail_neighbors.push_back(place);
                        not_placed = false;
                    }
                    
                }
                
                
                for (int j=0; j<=update_max; ++j) {
                    p->update();
                }
                recalculate_fitness(*p, ea);
                
                double sf =get<STRIPE_FIT>(*p);
                fit(sf);
                
                
                //df.endl();
                
            }
            
            
            df.write(mean(fit))
            .write(max(fit))
            .write(variance(fit))
            .write(sqrt(variance(fit)));
            df.endl();
        }
    }
}
#endif
