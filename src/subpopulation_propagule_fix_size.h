
#ifndef _SUBPOPULATION_PROPAGULE_FIX_SIZE_H_
#define _SUBPOPULATION_PROPAGULE_FIX_SIZE_H_

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
struct subpopulation_propagule_fix_size {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        // the number of parents selected is the propagule size or 1, if
        // the propagule's composition is clonal.
        std::size_t num_prop = get<NUM_PROPAGULE_GERM>(mea,1);
        assert(num_prop > 0);
        if (num_prop > parents[0]->size()) {
            num_prop = parents[0]->size();
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
            if  ((*j)->alive() &&
                 (get<IS_PROPAGULE>(**j,0) == 2) &&
                 (get<PARENT>(**j,0) > 0)) {
                typename MEA::subpopulation_type::genome_type r((*j)->genome().begin(),
                             
                                                                (*j)->genome().begin()+(*j)->hw().original_size());
                typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                
                inherits_from(**j, *q, *p);
                
                
                std::size_t t = p->rng()(open_pos.size());
                std::size_t pos = open_pos[t];
                open_pos.erase(open_pos.begin() + t);
                
                
                //iterator insert_at(iterator i, individual_ptr_type x, const position_type& pos) {
                p->insert_at(p->end(), q, p->env().location(pos).position());
                
                
                // Rotate the org to a random direction...
                int n = p->rng()(8);
                position_type& p1 = q->position();
                for (int q=0; q<=n; ++q){
                    int blah =0;
                    p1.rotate_cw();
                    
                }
                
                
                
                (*j)->alive() = false;
                parents[0]->events().death(**j,*(parents[0]));
                
                
                ++num_moved;
                
            }
            if (num_moved == num_prop) {
                break;
            }
        }
        
        offspring.insert(offspring.end(),p);
    }
};

#endif



