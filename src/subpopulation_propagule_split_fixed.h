/* subpopulation_propagule_split_fixed.h
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
#ifndef _SUBPOPULATION_PROPAGULE_SPLIT_FIXED_H_
#define _SUBPOPULATION_PROPAGULE_SPLIT_FIXED_H_

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
struct subpopulation_propagule_split_fixed {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        // the number of parents selected is the propagule size or 1, if
        // the propagule's composition is clonal.
        std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
        assert(prop_size > 0);
        if (prop_size > (parents[0]->size() / 2.0)) {
            return;
        }
        
//        if (prop_size > get<PROP_COUNT>(*(parents[0]),0)) {
//            prop_size = get<PROP_COUNT>(*(parents[0]),0);
//        }
//        
//        
        
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
            if ((*j)->alive()) {

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
                get<PROP_COUNT>(*(parents[0]))--;
            }
            if (num_moved == prop_size) {
                break;
            }
        }
        
        // update parent prop count
        offspring.insert(offspring.end(),p);
    }

};

struct subpopulation_propagule_split_fixed_partial_pick {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        // the number of parents selected is the propagule size or 1, if
        // the propagule's composition is clonal.
        std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
        assert(prop_size > 0);
        if (prop_size > parents[0]->size()) {
            prop_size = parents[0]->size();
        }
        
        
        double num_prop = get<PROP_COUNT>(*(parents[0]),0);
                               
        
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
        
        // pick the propagules first.
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
                get<PROP_COUNT>(*(parents[0]))--;
            }
            
            // If we moved the total number of the propagule size OR we have moved all of the propagule cells...
            if ((num_moved == prop_size) || (num_moved == num_prop)){
                break;
            }
        }
        
        // if we need more propagule cells...
        if (num_moved < prop_size) {
                                   
            for(typename propagule_type::iterator j=parents[0]->population().begin(); j!=parents[0]->population().end(); ++j) {
                // we are specifically NOT picking propagule cells (since they were already selected)
                if ((get<IS_PROPAGULE>(**j, 0) != 2) && (*j)->alive()) {
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
                }
                
                // If we moved the total number of the propagule size OR we have moved all of the propagule cells...
                if (num_moved == prop_size) {
                    break;
                }
            }
        }
                               
        // update parent prop count
        offspring.insert(offspring.end(),p);
    }
};



struct subpopulation_propagule_split_fixed_pick {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        // the number of parents selected is the propagule size or 1, if
        // the propagule's composition is clonal.
        std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
        assert(prop_size > 0);
        if (prop_size > parents[0]->size()) {
            prop_size = parents[0]->size();
        }
        
        
        double num_prop = ceil(get<PROP_COUNT>(*(parents[0]),0) / 2.0);

        if ((prop_size > num_prop) && (num_prop > 0)) {
            prop_size = num_prop;
        }
        
        
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
                get<PROP_COUNT>(*(parents[0]))--;
            }
            if (num_moved == prop_size) {
                break;
            }
        }
        
        // update parent prop count
        offspring.insert(offspring.end(),p);
    }
};


#endif
