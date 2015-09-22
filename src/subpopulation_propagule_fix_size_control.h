
#ifndef _SUBPOPULATION_PROPAGULE_FIX_SIZE_CONTROL_H_
#define _SUBPOPULATION_PROPAGULE_FIX_SIZE_CONTROL_H_

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
struct subpopulation_propagule_fix_size_control {
    std::size_t capacity() const { return 1; }
    
    template <typename Population, typename MEA>
    void operator()(Population& parents, Population& offspring, MEA& mea) {
        
        
        // get a new subpopulation:
        typename MEA::individual_ptr_type p = mea.make_individual();
        p->initialize(mea.md());
        p->reset_rng(mea.rng().seed());
        
        
        typedef typename MEA::subpopulation_type::population_type propagule_type;
        

        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(mea));

        /*
        // if there are no founders, work around it and add the first member of the pop.
        // This SHOULD just happen during the first replication cycle
        if (parents[0]->traits().founder()->population().size() == 0) {
            int x = 5;
        } else {
            int i = 1+1;
        }*/
        
        // shuffle the founder population
        //std::random_shuffle(parents[0]->traits().founder()->population().begin(), parents[0]->traits().founder()->population().end(), parents[0]->rng());
        
        
        
        int s = get<POPULATION_SIZE>(mea);
        std::vector<int> open_pos (s);
        for( int n = 0 ; n < s ; ++n ) {
            open_pos[ n ] = n;
        }
        
        // first founder is the propagule.
        typename MEA::subpopulation_type::genome_type r(parents[0]->population()[0]->genome().begin(), parents[0]->population()[0]->genome().begin()+parents[0]->population()[0]->hw().original_size());
        
        typename MEA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
        
        // mutate
        mutate(*q,m,*p);
        
        inherits_from(*(parents[0]->population()[0]), *q, *p);
        
        for (int i=0; i< get<NUM_PROPAGULE_CELL>(mea,1);  ++i) {
            typename MEA::subpopulation_type::individual_ptr_type o(q);

            std::size_t t = p->rng()(open_pos.size());
            std::size_t pos = open_pos[t];
            open_pos.erase(open_pos.begin() + t);
            
            p->insert_at(p->end(), o, p->env().location(pos).position());
            
            // Rotate the org to a random direction...
            int n = p->rng()(8);
            position_type& p1 = o->position();
            for (int q=0; q<=n; ++q){
                p1.rotate_cw();
                
            }
        }
        
        offspring.insert(offspring.end(),p);

    }
};

#endif



