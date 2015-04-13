#ifndef _PS_MOVIE_H_
#define _PS_MOVIE_H_

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
        
        /* Run through each multicell... 
         - reset it
         - run for competition period
         - save if its fitness is the best!
         - write the best to a file.
         */
        LIBEA_ANALYSIS_TOOL(movie_for_competitions) {
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            int max_fit = 0;
            typename EA::individual_type best;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                typename EA::individual_type tmp(*i->traits().founder());
                
                for (int j=0; j<=update_max; ++j) {
                    tmp.update();
                }
                
                // Run for prescribed number of updates...
                
                
                // Calc fitness for each subpop
                eval_permute_stripes(tmp);
                
                // copy the stripe fit to the accumulator and also the subpop
                double sf =get<STRIPE_FIT>(tmp);
                
                if (sf > max_fit) {
                    max_fit = sf;
                    best = *i;
                }
            }
            

            datafile df("movie.dat");
            df.write(get<SPATIAL_X>(ea));
            df.write(get<SPATIAL_Y>(ea));
            df.endl();

            typename EA::individual_type best_founder(*best.traits().founder());

            for (int j=0; j<=update_max; ++j) {
                best_founder.update();
                df.write(j);
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                        
                        //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
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
            
        }


        

        
        
    }
}
#endif
