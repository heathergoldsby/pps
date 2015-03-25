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
                .add_field("mean_propagule_size")
                .add_field("mean_multicell_size");
//                .add_field("mean_fitness")
//                .add_field("max_fitness");
            }
            
            virtual ~propagule_dat() {
            }
            
            virtual void operator()(EA& ea) {
                using namespace boost::accumulators;
                accumulator_set<double, stats<tag::mean> > propagule_size;
                accumulator_set<double, stats<tag::mean> > multicell_size;
                
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    double num_prop = ceil(get<PROP_COUNT>(*i,0) / 2.0);

                    propagule_size(num_prop);
                    multicell_size(i->size());
                }
                
                _df.write(ea.current_update())
                .write(mean(propagule_size))
                .write(mean(multicell_size))
                .endl();
            }
            
            datafile _df;
        };
        
    } // datafiles
} // ealib

#endif
