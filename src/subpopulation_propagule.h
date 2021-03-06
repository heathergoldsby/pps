/* subpopulation_propagule.h
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
#ifndef _SUBPOPULATION_PROPAGULE_H_
#define _SUBPOPULATION_PROPAGULE_H_

#include <algorithm>
#include <utility>
#include <iterator>


#include <ea/metadata.h>

/*! This recombination operator is metapopulation-specific; it copies a
 propagule from a parent subpopulation, mutates it, and adds it to an
 offspring subpopulation.
 */
struct subpopulation_propagule {
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
        
        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        // figure out which individuals from the parent comprise the propagule:
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        propagule_type propagule;
        mea.rng().sample_without_replacement(parents[0]->population().begin(),
                                             parents[0]->population().end(),
                                             std::back_inserter(propagule),
                                             prop_size);
        
        
        int s = get<POPULATION_SIZE>(mea);
        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));
        
        std::vector<int> open_pos (s);
        for( int n = 0 ; n < s ; ++n ) {
            open_pos[ n ] = n;
        }
        

        // now add a new individual built from each of the propagules to the
        // subpopulation:
        for(typename propagule_type::iterator i=propagule.begin(); i!=propagule.end(); ++i) {
            // grab the original part of the propagule's genome; note that it could have been
            // changed (implicit-like mutations):
            typename MEA::subpopulation_type::genome_type r((*i)->genome().begin(),
                                                            (*i)->genome().begin()+(*i)->hw().original_size());
            typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
            
            inherits_from(**i, *q, *p);
            
            // mutate
            mutate(*q,m,*p);
            
            
            std::size_t t = p->rng()(open_pos.size());
            std::size_t pos = open_pos[t];
            open_pos.erase(open_pos.begin() + t);
            
            p->insert_at(p->end(),q, p->env().location(pos).position());

            
            for (int k=1; k<get<NUM_PROPAGULE_CELL>(mea); ++k) {
                // check if this is a valid way to copy an individual
                typename MEA::subpopulation_type::individual_ptr_type o(q);

                t = p->rng()(open_pos.size());
                pos = open_pos[t];
                open_pos.erase(open_pos.begin() + t);
                p->insert_at(p->end(), o, p->env().location(pos).position());

             }

        }
        
        
        
            for (int x=0; x < get<SPATIAL_X>(mea); ++x) {
                for (int y=0; y<get<SPATIAL_Y>(mea); ++y){
                    
                    //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
                    typename MEA::subpopulation_type::environment_type::location_type l = (p->env().location(x,y));
                    
                    
//                    if (l.occupied()) {
//                        std::cout << x << " " << y << std::endl;
//                    }
                    
                }
            }
//        std::cout << std::endl << std::endl;

        
        
        offspring.insert(offspring.end(),p);
    }
};

#endif
