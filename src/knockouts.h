#ifndef _PS_KO_H_
#define _PS_KO_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
//#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "stripes.h"


namespace ealib {
    namespace analysis {
        
        /* Run through each multicell...
         - reset it
         - run for competition period
         - analyze knockouts
         */
        LIBEA_ANALYSIS_TOOL(knockouts_for_competition) {
            
            using namespace boost::accumulators;
            accumulator_set<double, stats<tag::mean> > c;
            accumulator_set<double, stats<tag::mean> > r;
            accumulator_set<double, stats<tag::mean> > o;
            accumulator_set<double, stats<tag::mean> > a;
            accumulator_set<double, stats<tag::mean> > n;
            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                typename EA::individual_type control(*i->traits().founder());
                
                typename EA::individual_type ko_rx(*i->traits().founder());
                knockout<instructions::rx_msg,instructions::nop_x>(ko_rx);
                
                typename EA::individual_type ko_one_prop(*i->traits().founder());
                knockout<deploy_one_propagule,instructions::nop_x>(ko_one_prop);

                
                typename EA::individual_type ko_prop_absent(*i->traits().founder());
                knockout<if_prop_cell_absent,instructions::nop_x>(ko_prop_absent);

                typename EA::individual_type ko_is_neighbor(*i->traits().founder());
                knockout<instructions::is_neighbor,instructions::nop_x>(ko_is_neighbor);

                
                for (int j=0; j<=update_max; ++j) {
                    control.update();
                    ko_rx.update();
                    ko_one_prop.update();
                    ko_prop_absent.update();
                    ko_is_neighbor.update();
                }
                
                // Run for prescribed number of updates...
                
                
                // Calc fitness for each subpop
                eval_permute_stripes(control);
                eval_permute_stripes(ko_rx);
                eval_permute_stripes(ko_one_prop);
                eval_permute_stripes(ko_prop_absent);
                eval_permute_stripes(ko_is_neighbor);
                
                c(get<STRIPE_FIT>(control,0));
                r(get<STRIPE_FIT>(ko_rx,0));
                o(get<STRIPE_FIT>(ko_one_prop,0));
                a(get<STRIPE_FIT>(ko_prop_absent,0));
                n(get<STRIPE_FIT>(ko_is_neighbor,0));
                
            }
            
            datafile df("knockouts.dat");
            df.add_field("no_ko")
            .add_field("rx_ko")
            .add_field("one_prop_ko")
            .add_field("prop_absent_ko")
            .add_field("is_neighbor_ko");
                
            df.write(mean(c));
            df.write(mean(r));
            df.write(mean(o));
            df.write(mean(a));
            df.write(mean(n));
                
            df.endl();
//            
        
        }
        
        
        
        
        
        
    }
}
#endif
