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
                eval_two_stripes(tmp);
                
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
                            if (lt == "or") {
                                df.write("3");
                            }
                            if (lt == "ornot") {
                                df.write("4");
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


        
        LIBEA_ANALYSIS_TOOL(movie_bullseye) {
            
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
                eval_bullseye(tmp);
                
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
                            if (lt == "or") {
                                df.write("3");
                            }
                            if (lt == "ornot") {
                                df.write("4");
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

        
        LIBEA_ANALYSIS_TOOL(movie_n_neighbors) {
            
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
                eval_n_neighbors(tmp);
                
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
                            if (lt == "or") {
                                df.write("3");
                            }
                            if (lt == "ornot") {
                                df.write("4");
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

    
//        /*! lod_movie
//         */
//        LIBEA_ANALYSIS_TOOL(movie) {
//
//            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
//            
//            typename line_of_descent<EA>::reverse_iterator i=lod.rbegin(); ++i;
//            datafile df("movie.dat");
//            typename EA::individual_ptr_type control_sp = ea.make_individual();
//            control_sp->ea().rng().reset(get<RNG_SEED>((i->ea())));
//            
//        
//            typename EA::individual_type::ea_type::individual_type g (i->ea().founder());
//            typename EA::individual_type::ea_type::individual_ptr_type o = i->ea().copy_individual(g);
//            o->hw().initialize();
//
//            
//            control_sp->ea().insert(control_sp->ea().end(), o);
//            
//
//            
//            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
//
//            df.write(get<SPATIAL_X>(ea));
//            df.write(get<SPATIAL_Y>(ea));
//            df.endl();
//            for (int i=0; i<=update_max; ++i) {
//                control_sp->ea().update();
//                df.write(i);
//                // grab info based on location...
//                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
//                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
//                        typename EA::individual_type::ea_type::environment_type::location_ptr_type l = control_sp->ea().env().location(x,y);
//                        if (l->occupied()) {
//                            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
//
//                            if(lt == "not") {
//                                df.write("1");
//                            }
//                            if (lt == "nand") {
//                                df.write("2");
//                            }
//                            if (lt == "") {
//                                df.write("0");
//                            }
//                            
//            
//                        } else {
//                            df.write("-1");
//                        }
//                    
//                    }
//                }
//                df.endl();
//
//            }
//            
//            df.endl();
//            
//             }
    }
}
#endif
