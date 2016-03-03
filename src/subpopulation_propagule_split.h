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
#ifndef _SUBPOPULATION_PROPAGULE_SPLIT_H_
#define _SUBPOPULATION_PROPAGULE_SPLIT_H_

#include <algorithm>
#include <utility>

#include <ea/metadata.h>
#include "stripes.h"

using namespace ealib;


/*! This recombination operator is metapopulation-specific; it copies a
 propagule from a parent subpopulation, mutates it, and adds it to an
 offspring subpopulation.
 
 Cells are flagged as propagules using the create_propagule and deploy_propagule methods.
 */
struct subpopulation_propagule_split {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
    
        double num_prop = ceil(get<PROP_COUNT>(*(parents[0]),0) / 2.0);
        
        if (num_prop == 0) return;
        

        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        // shuffle the population
        std::random_shuffle(parents[0]->population().begin(), parents[0]->population().end(), parents[0]->rng());
                
        int num_moved = 0;
        int s = get<POPULATION_SIZE>(mea);
        std::vector<int> open_pos (s);
        for( int n = 0 ; n < s ; ++n ) {
            open_pos[ n ] = n;
        }
        
        for(typename propagule_type::iterator j=parents[0]->population().begin(); j!=parents[0]->population().end(); ++j) {
            if ((get<IS_PROPAGULE>(**j, 0) == 2) && (*j)->alive()) {
                typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                (*j)->genome().begin()+(*j)->hw().original_size());
                typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                
                inherits_from(**j, *q, *p);

                
                std::size_t t = p->rng()(open_pos.size());
                std::size_t pos = open_pos[t];
                open_pos.erase(open_pos.begin() + t);
                
                p->insert_at(p->end(), q, p->env().location(pos).position());
                
                (*j)->alive() = false;
                parents[0]->events().death(**j,*(parents[0]));
                
                
                ++num_moved;
                get<PROP_COUNT>*(parents[0]))--; 
            }
            if (num_moved == num_prop) {
                break;
            }
        }
        
        // update parent prop count
        offspring.insert(offspring.end(),p);
    }
};

#endif
