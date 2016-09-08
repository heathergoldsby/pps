/* propagule.h
 *
 * This file is part of EALib.
 *
 * Copyright 2014 David B. Knoester.
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
#ifndef _EA_DATAFILES_PROPAGULE_H_
#define _EA_DATAFILES_PROPAGULE_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>


#include <ea/datafile.h>
#include <ea/traits.h>


namespace ealib {
    namespace datafiles {
        
        /*! Datafile for mean generation and min, mean, and max fitness.
         */
        template <typename EA>
        struct propagule_dat : record_statistics_event<EA> {
            propagule_dat(EA& ea) : record_statistics_event<EA>(ea), _df("propagule.dat") {
                _df.add_field("update")
                .add_field("mean_num_prop_cells")
                .add_field("mean_num_prop_cells_act")
                .add_field("num_zero_prop_cells")
                .add_field("mean_multicell_size")
                .add_field("mc_rep_count_over_0")
                .add_field("mean_size_mc_rep_count_over_0") // this one is the mean number of births for those over 0.
                .add_field("mean_prop_size")
                .add_field("min_prop_size")
                .add_field("max_prop_size")
                .add_field("mean_stripes_fit")
                .add_field("max_stripes_fit");

            }
            
            virtual ~propagule_dat() {
            }
            
            virtual void operator()(EA& ea) {
                using namespace boost::accumulators;
                accumulator_set<double, stats<tag::mean> > num_prop_cells;
                accumulator_set<double, stats<tag::mean> > num_prop_cells_act;

                accumulator_set<double, stats<tag::mean, tag::max> > fit;
                accumulator_set<double, stats<tag::mean, tag::min, tag::max> > prop_size;

                accumulator_set<double, stats<tag::mean> > multicell_size;
                accumulator_set<double, stats<tag::mean> > size_mc_rep_count_over_0;

                int rep_count_over_0 = 0;
                int num_zero_prop = 0;
                
                typedef typename EA::subpopulation_type::population_type propagule_type;
                
                int prop_count;
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    accumulator_set<double, stats<tag::mean> > subpop_prop_size;

                    prop_count = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            if (get<IS_PROPAGULE>(**j,0) == 1) {
                                prop_count++;
                                
                            }
                            if (exists<PROPAGULE_SIZE>(**j)) {
                                subpop_prop_size(get<PROPAGULE_SIZE>(**j));
                            }

                        }
                        
                    }
                    
                    if (prop_count == 1) {
                        num_zero_prop++;
                    }
                    
                    
                    
                    num_prop_cells(get<PROP_COUNT>(*i,0));
                    num_prop_cells_act(prop_count);
                    multicell_size(i->size());
                    prop_size(floor(mean(subpop_prop_size)));
                    fit(get<STRIPE_FIT>(*i));
                    
                    if ((get<REP_COUNT>(*i,0)) > 0) {
                        ++rep_count_over_0;
                        size_mc_rep_count_over_0(get<REP_COUNT>(*i));
                    }
                }
                
                _df.write(ea.current_update())
                .write(mean(num_prop_cells))
                .write(mean(num_prop_cells_act))
                .write(num_zero_prop)
                .write(mean(multicell_size))
                .write(rep_count_over_0)
                .write(mean(size_mc_rep_count_over_0))
                .write(mean(prop_size))
                .write(min(prop_size))
                .write(max(prop_size))
                .write(mean(fit))
                .write(max(fit))
                .endl();
            }
            
            datafile _df;
        };
        
    } // datafiles
} // ealib

#endif
