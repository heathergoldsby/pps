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


#include "unique_number.h"
/*
 recalculate_fitness(*i, ea);
 fs(static_cast<int>(ealib::fitness(*i,ea)));
 */

namespace ealib {
    namespace analysis {
        
        LIBEA_ANALYSIS_TOOL(movie_for_competitions_spatial) {
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            int max_fit = 0;
            typename EA::individual_type best;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                // not preserving location
                typename EA::individual_ptr_type p = ea.make_individual();
                p->initialize(ea.md());
                //p->reset_rn(ea.rng().seed());
                
                //typename EA::individual_type tmp(*i->traits().founder());
                typedef typename EA::subpopulation_type::population_type propagule_type;
                for (typename propagule_type::iterator j=(*i->traits().founder()).population().begin(); j!= (*i->traits().founder()).population().end(); ++j) {
                    int pos = 0;
                    typename EA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                                    (*j)->genome().begin()+(*j)->hw().original_size());
                    typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);

                    inherits_from(**j, *q, *p);

                    p->insert_at(p->end(),q, p->env().location((*j)->position()).position());

                }
                
                
                for (int j=0; j<=update_max; ++j) {
                    p->update();
                }
                recalculate_fitness(*p, ea);
                
                
                
                
                // Calc fitness for each subpop
                //                eval_two_stripes(tmp);
                
                // copy the stripe fit to the accumulator and also the subpop
                double sf =get<STRIPE_FIT>(*p);
                
                if (sf > max_fit) {
                    max_fit = sf;
                    best = *i;
                }
            }
            
            
            datafile df("movie.dat");
            df.write(get<SPATIAL_X>(ea));
            df.write(get<SPATIAL_Y>(ea));
            df.endl();
            
            int pop_size = get<POPULATION_SIZE>(ea);
            
            std::vector<int> used_pos;
            std::vector<int> used_pos_with_avail_neighbors;
            
            //typename EA::individual_type best_founder = *best.traits().founder();
            // not preserving location
            typename EA::individual_ptr_type p = ea.make_individual();
            p->initialize(ea.md());
            //p->reset_rn(ea.rng().seed());
            
            
            int pos = ea.rng().uniform_integer(0,pop_size);
            used_pos.push_back(pos);
            used_pos_with_avail_neighbors.push_back(pos);
            
            
            //typename EA::individual_type tmp(*i->traits().founder());
            typedef typename EA::subpopulation_type::population_type propagule_type;
            for (typename propagule_type::iterator j=(*best.traits().founder()).population().begin(); j!= (*best.traits().founder()).population().end(); ++j) {
                int pos = 0;
                typename EA::subpopulation_type::genome_type r((*j)->genome().begin(),
                                                               (*j)->genome().begin()+(*j)->hw().original_size());
                typename EA::subpopulation_type::individual_ptr_type q = p->make_individual(r);
                
                inherits_from(**j, *q, *p);
                
  
                int max_x =get<SPATIAL_X>(ea);
                

                    
                    bool not_placed = true;
                    int place = -1;
                    
                    while (not_placed) {
                        int n = ea.rng().uniform_integer(0,used_pos_with_avail_neighbors.size());
                        int neighbor = used_pos_with_avail_neighbors[n];
                        int dir_try = 0;
                        int dir = ea.rng().uniform_integer(0,3);
                        
                        
                        
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
                                if ((neighbor + max_x) >= pop_size) { continue; }
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
                    
                }
                
    

            
            for (int j=0; j<=update_max; ++j) {
                p->update();
                df.write(j);
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                        
                        //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
                        typename EA::individual_type::environment_type::location_type* l = &p->env().location(x,y);
                        
                        
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
            
            //df.endl();
            
        }
        

        
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
                
                // not preserving location
                typename EA::individual_type tmp(*i->traits().founder());
                
                for (int j=0; j<=update_max; ++j) {
                    tmp.update();
                }
                recalculate_fitness(tmp, ea);
                
                
                
                
                // Calc fitness for each subpop
//                eval_two_stripes(tmp);
                
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

//            typename EA::individual_type best_founder(*best.traits().founder());
            typename EA::individual_type best_founder = *best.traits().founder();

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

        
                
        LIBEA_ANALYSIS_TOOL(movie_unique_number) {
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            int max_fit = 0;
            typename EA::individual_type best;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                typename EA::individual_type tmp(*i->traits().founder());
                
                for (int j=0; j<=update_max; ++j) {
                    tmp.update();
                }
                recalculate_fitness(tmp, ea);
                
                
                
                
                // Calc fitness for each subpop
                //                eval_two_stripes(tmp);
                
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
                            int op = get<OPINION>(*(l->inhabitant()),0);
                            df.write(op);
                            
                        } else {
                            df.write(-1);
                        }
                        
                    }
                }
                df.endl();
                
            }
            
            df.endl();
            
        }
        

        LIBEA_ANALYSIS_TOOL(print_dom) {
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            int max_fit = 0;
            typename EA::individual_type best;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                typename EA::individual_type tmp(*i->traits().founder());
                
                for (int j=0; j<=update_max; ++j) {
                    tmp.update();
                }
                recalculate_fitness(tmp, ea);
                
                
                
                
                // Calc fitness for each subpop
                //                eval_two_stripes(tmp);
                
                // copy the stripe fit to the accumulator and also the subpop
                double sf =get<STRIPE_FIT>(tmp);
                
                if (sf > max_fit) {
                    max_fit = sf;
                    best = *i;
                }
            }
            
            
            datafile df("dominant.dat");
            
            // (*j)->genome().begin()
            
            typename EA::subpopulation_type::genome_type r((best[0])->genome().begin(),
                                                            (best[0])->genome().begin()+(best[0])->hw().original_size());
//            
//            // iterate through a genome?
//            for (k = r.begin(); k < r.end(); ++k) {
//                
//            }
//            
            
            //(*j)[0]->genome();
            
//            df.write(get<SPATIAL_X>(ea));
//            df.write(get<SPATIAL_Y>(ea));
//            df.endl();
            
            //typename EA::individual_type best_founder(*best.traits().founder());
            
//            for (int j=0; j<=update_max; ++j) {
//                best_founder.update();
//                df.write(j);
//                // grab info based on location...
//                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
//                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
//                        
//                        //typename EA::individual_type::environment_type::location_type* l = &best_founder.traits().founder()->env().location(x,y);
//                        typename EA::individual_type::environment_type::location_type* l = &best_founder.env().location(x,y);
//                        
//                        
//                        if (l->occupied()) {
//                            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
//                            
//                            if(lt == "not") {
//                                df.write("1");
//                            }
//                            if (lt == "nand") {
//                                df.write("2");
//                            }
//                            if (lt == "or") {
//                                df.write("3");
//                            }
//                            if (lt == "ornot") {
//                                df.write("4");
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
