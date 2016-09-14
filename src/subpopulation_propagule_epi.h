/* subpopulation_propagule_split.h
 *
 * This file is part of EALib.
 *
 * Copyright 2014 David B. Knoester, Heather J. Goldsby.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _SUBPOPULATION_PROPAGULE_EPI_H_
#define _SUBPOPULATION_PROPAGULE_EPI_H_

#include <algorithm>
#include <utility>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <ea/metadata.h>
#include "stripes.h"


using namespace ealib;


/*! This recombination operator is metapopulation-specific; it copies a
 propagule from a parent subpopulation, mutates it, and adds it to an
 offspring subpopulation.
 
 Cells are flagged as propagules using the create_propagule and deploy_propagule methods.
 */
struct subpopulation_propagule_epi {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        
        int p_size = get<NUM_PROPAGULE_CELL>(*parents[0]);
        
        if ((p_size == 0) || (p_size >= parents[0]->size())) {
            return;
        }
        
        
        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        
        int num_moved = 0;
        int s = get<POPULATION_SIZE>(mea);
        std::vector<int> open_pos (s);
        for( int n = 0 ; n < s ; ++n ) {
            open_pos[ n ] = n;
        }
        
        
        // first grab the propagule cells...
        for(typename propagule_type::iterator j=parents[0]->population().begin(); j!=parents[0]->population().end(); ++j) {
            if ((*j)->alive()) {
                
                typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                (*j)->genome().begin()+(*j)->hw().original_size());

                

                // If the cells are identical, we make copies of the first one rather than moving them all.
                if (get<PROPAGULE_IDENTICAL>(*parents[0],0)){
                    // for the first propagule cell make however many copies.
                    if (num_moved == 0) {
                        for(int i = 0; i < p_size; ++i) {
                            std::size_t t = p->rng()(open_pos.size());
                            std::size_t pos = open_pos[t];
                            open_pos.erase(open_pos.begin() + t);
                
                            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                        
                            inherits_from(**j, *q, *p);
                            p->insert_at(p->end(), q, p->env().location(pos).position());
                        }
                    } 
                } else {
                    std::size_t t = p->rng()(open_pos.size());
                    std::size_t pos = open_pos[t];
                    open_pos.erase(open_pos.begin() + t);
                    
                    typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                    
                    inherits_from(**j, *q, *p);
                    p->insert_at(p->end(), q, p->env().location(pos).position());

                }
                
                
                // If the propagule cell is moved. (propagule_copy == 0)
                if (get<PROPAGULE_COPY>(*parents[0],0) == 0) {
                    (*j)->alive() = false;
                    parents[0]->events().death(**j,*(parents[0]));
                }
                
                
                ++num_moved;
            }
            if (num_moved == p_size) {
                break;
            }
        }
        
        // set propagule size for
        
        int new_p_size = p_size;
        if (mea.rng().p(get<GERM_MUTATION_PER_SITE_P>(mea))) {
            new_p_size = mea.rng()(19);
        }
        
        put<NUM_PROPAGULE_CELL>(new_p_size, *p);
        get<REP_COUNT>(*(parents[0]),0) ++;
        
        
        // update parent prop count
        offspring.insert(offspring.end(),p);
    }
};

#endif
