/* subpopulation_propagule_plane.h
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
#ifndef _SUBPOPULATION_PROPAGULE_PLANE_H_
#define _SUBPOPULATION_PROPAGULE_PLANE_H_

#include <algorithm>
#include <utility>

#include <ea/metadata.h>

/*! This recombination operator is metapopulation-specific; it copies a
 propagule from a parent subpopulation, mutates it, and adds it to an
 offspring subpopulation. 
 
 Cells are flagged as propagules using the create_propagule and deploy_propagule methods.
 */
struct subpopulation_propagule_plane {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {

        
        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        // figure out which individuals from the parent comprise the propagule:
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        propagule_type propagule;

        for(typename propagule_type::iterator i=parents[0]->population().begin(); i!=parents[0]->population().end(); ++i) {
            if (get<IS_PROPAGULE>(**i, 0) == 2) {
                propagule.push_back(*i);
            }
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
            mutate(*q,*p);
            
            p->insert(p->end(),q);
        }
        
        offspring.insert(offspring.end(),p);
    }
};

#endif
