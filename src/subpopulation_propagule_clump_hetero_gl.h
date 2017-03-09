/*
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
#ifndef _SUBPOPULATION_PROPAGULE_CLUMP_HETERO_GL_H_
#define _SUBPOPULATION_PROPAGULE_CLUMP_HETERO_GL_H_

#include <algorithm>
#include <utility>
#include <iterator>


#include <ea/metadata.h>



/*! An organism rotates to face its parent....
 */
template <typename EA>
struct gl_birth_event : birth_event<EA> {
    
    //! Constructor.
    gl_birth_event(EA& ea) : birth_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~gl_birth_event() {
    }
    
    /*! Called for every inheritance event. We are using the orientation of the first parent...
     */
    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
                            typename EA::individual_type& parent, // individual parent
                            EA& ea) {
        put<DIVISION_ROUND>((get<DIVISION_ROUND>(parent,0) + 1), offspring);
        
    }
};





/*! This recombination operator is metapopulation-specific; it copies a
 propagule from a parent subpopulation, mutates it, and adds it to an
 offspring subpopulation.
 */
struct subpopulation_propagule_clump_hetero_gl {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        // the number of parents selected is the propagule size or 1, if
        // the propagule's composition is clonal.
        std::size_t prop_size = get<NUM_PROPAGULE_GERM>(mea,1);
        
        // sterility clause
        if (prop_size == parents[0]->size()) {
            return;
        }
        
        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        
        
        // figure out which individuals from the parent comprise the propagule:
        // grab the first N cells....
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        propagule_type propagule;
        /*
        mea.rng().sample_without_replacement(parents[0]->population().begin(),
                                             parents[0]->population().end(),
                                             std::back_inserter(propagule),
                                             prop_size);*/
        int round_of_picks = 0;
        while (propagule.size() < prop_size) {
            for(typename propagule_type::iterator i=parents[0]->population().begin(); i!=parents[0]->population().end(); ++i) {

            // go through cells... add those in this division round...
                int dr = get<DIVISION_ROUND>(*i,0);
                if (dr == round_of_picks) {
                    propagule->insert(*i);
                }
                if (propagule.size() == prop_size) {
                    break;
                }
            }
            round_of_picks ++;
        }
        
        
        int pop_size = get<POPULATION_SIZE>(mea);
        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
        
        std::vector<int> used_pos;
        std::vector<int> used_pos_with_avail_neighbors;
        
        
        int max_x =get<SPATIAL_X>(mea);
        
        // now add a new individual built from each of the propagules to the
        // subpopulation:
        for(typename propagule_type::iterator i=propagule.begin(); i!=propagule.end(); ++i) {
            // grab the original part of the propagule's genome; note that it could have been
            // changed (implicit-like mutations):
            typename MEA::subpopulation_type::genome_type r((*i)->genome().begin(),
                                                            (*i)->genome().begin()+(*i)->hw().original_size());
            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
            
            inherits_from(**i, *q, *p);
            
            mutate(*q,m,*p);
            
            put<DIVISION_ROUND>(0,*q);
            
            if (used_pos.size() == 0) {
                int pos = mea.rng().uniform_integer(0,pop_size);
                used_pos.push_back(pos);
                used_pos_with_avail_neighbors.push_back(pos);
                p->insert_at(p->end(),q, p->env().location(pos).position());
            } else {
                
                //
                //            int max_x =get<SPATIAL_X>(mea);
                //
                
                typename MEA::subpopulation_type::individual_ptr_type o(q);
                
                bool not_placed = true;
                int place = -1;
                
                while (not_placed) {
                    int n = mea.rng().uniform_integer(0,used_pos_with_avail_neighbors.size());
                    int neighbor = used_pos_with_avail_neighbors[n];
                    int dir_try = 0;
                    int dir = mea.rng().uniform_integer(0,3);
                    
                    
                    
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
                            if ((neighbor + max_x) >= get<POPULATION_SIZE>(mea)) { continue; }
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
                
                //}
            }
        }
        
        
        
        offspring.insert(offspring.end(),p);
    }
};

#endif
