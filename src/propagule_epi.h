/* propagule_epi.h
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
#ifndef _EA_DATAFILES_PROPAGULE_EPI_H_
#define _EA_DATAFILES_PROPAGULE_EPI_H_

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
        struct propagule_epi_dat : record_statistics_event<EA> {
            propagule_epi_dat(EA& ea) : record_statistics_event<EA>(ea), _df("propagule_epi.dat") {
                _df.add_field("update")
                .add_field("mean_multicell_size")
                .add_field("mean_prop_size")
                .add_field("min_prop_size")
                .add_field("max_prop_size")
                .add_field("mean_stripes_fit")
                .add_field("max_stripes_fit");
                
            }
            
            virtual ~propagule_epi_dat() {
            }
            
            virtual void operator()(EA& ea) {
                using namespace boost::accumulators;
                
                accumulator_set<double, stats<tag::mean, tag::max> > fit;
                accumulator_set<double, stats<tag::mean, tag::min, tag::max> > prop_size;
                
                accumulator_set<double, stats<tag::mean> > multicell_size;
                

                typedef typename EA::subpopulation_type::population_type propagule_type;
                
                int prop_count;
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    
                    
                    prop_size(get<NUM_PROPAGULE_CELL>(*i,0));
                    multicell_size(i->size());
                    fit(get<STRIPE_FIT>(*i));
                }
                
                _df.write(ea.current_update())
                .write(mean(multicell_size))
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
