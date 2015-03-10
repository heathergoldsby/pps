
#ifndef _PS_MOVIE_EVO_PLANE_H_
#define _PS_MOVIE_EVO_PLANE_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>


namespace ealib {
    namespace analysis {
        
        /* movie for evo plane */
        LIBEA_ANALYSIS_TOOL(movie_evo_plane) {
            
            int count_updates = 0;
            
            while (count_updates < 1000) {
                count_updates ++;
                ea.update();
            }
            
            
            
            count_updates = 0;
            
            while (count_updates < 10000) {
                count_updates ++;
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    i->update();
                    // Calc fitness for each subpop
                    eval_permute_stripes(*i);
                    
                    // copy the stripe fit to the accumulator and also the subpop
                    double sf =get<STRIPE_FIT>(*i);
                    put<STRIPE_FIT>(sf, *i);
                    get<MC_RESOURCE_UNITS>(*i,0) += sf;
                    
                    // track time since group rep
                    get<MULTICELL_REP_TIME>(*i,0) +=1;
                    
                    // figure out which individuals from the parent comprise the propagule:
                    typedef typename EA::subpopulation_type::population_type propagule_type;
                    propagule_type propagule;
                    
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if (get<IS_PROPAGULE>(**j, 0) == 2) {
                            propagule.push_back(*j);
                        }
                    }
                    
                    int prop_total_cost = get<PROPAGULE_BASE_COST>(*i) + propagule.size() * get<PROPAGULE_PER_CELL_COST>(*i);
                    
                    if ((get<MC_RESOURCE_UNITS>(*i) > prop_total_cost) && (propagule.size() > 0)){
                        
                        datafile df("movie.dat");
                        df.write(get<SPATIAL_X>(ea));
                        df.write(get<SPATIAL_Y>(ea));
                        df.endl();
                        
                        typename EA::individual_type best_founder(*i->traits().founder());
                        
                        int num_updates = get<MULTICELL_REP_TIME>(*i);
//                        for (int j=0; j<=num_updates; ++j) {
                        int update_count =0;
                        while (get<MC_RESOURCE_UNITS>(best_founder,0) < prop_total_cost) {
                            best_founder.update();
                            df.write(update_count);
                            ++update_count;
                            // grab info based on location...
                            for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                                for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                                    
                                    typename EA::individual_type::environment_type::location_type* l = &best_founder.env().location(x,y);
                                    
                                    
                                    if (l->occupied()) {
                                        std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
                                        
                                        if(lt == "not") {
                                            df.write("1");
                                        }
                                        if (lt == "nand") {
                                            df.write("2");
                                        }
                                        if (lt == "") {
                                            df.write("0");
                                        }
                                        
                                        
                                    } else {
                                        df.write("-1");
                                    }
                                    
                                }
                            }
                            df.endl();
                            
                        }
                        
                        df.endl();
                        return;
                        
                        
                    }
                }
            }
        }
    }
}
#endif

