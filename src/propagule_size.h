/* propagule_size.h
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
#ifndef _EA_DATAFILES_PROPAGULE_SIZE_H_
#define _EA_DATAFILES_PROPAGULE_SIZE_H_

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
        struct propagule_size_dat : record_statistics_event<EA> {
            propagule_size_dat(EA& ea) : record_statistics_event<EA>(ea), _df("propagule_size.dat") {
                _df.add_field("update")
                .add_field("psNA")
                .add_field("ps1")
                .add_field("ps2")
                .add_field("ps3")
                .add_field("ps4")
                .add_field("ps5")
                .add_field("ps6")
                .add_field("ps7")
                .add_field("ps9")
                .add_field("ps10")
                .add_field("ps11")
                .add_field("ps12")
                .add_field("ps13")
                .add_field("ps14")
                .add_field("ps15")
                .add_field("ps16")
                .add_field("ps17")
                .add_field("ps18");
                
            }
            
            virtual ~propagule_size_dat() {
            }
            
            virtual void operator()(EA& ea) {
                using namespace boost::accumulators;
                
                int psSizes[19] = {0};
                
                typedef typename EA::subpopulation_type::population_type propagule_type;
                
                int prop_count;
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    accumulator_set<double, stats<tag::mean> > subpop_prop_size;
                    
                    prop_count = 0;
                    for(typename propagule_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                        if ((*j)->alive()) {
                            
                            if (exists<PROPAGULE_SIZE>(**j)) {
                                int ps = get<PROPAGULE_SIZE>(**j);
                                psSizes[ps] += 1;
                            } else {
                                psSizes[0] += 1;
                            }
                            
                        }
                        
                    }
                    
                }
                
                _df.write(ea.current_update())
                .write(psSizes[0])
                .write(psSizes[1])
                .write(psSizes[2])
                .write(psSizes[3])
                .write(psSizes[4])
                .write(psSizes[5])
                .write(psSizes[6])
                .write(psSizes[7])
                .write(psSizes[8])
                .write(psSizes[9])
                .write(psSizes[10])
                .write(psSizes[11])
                .write(psSizes[12])
                .write(psSizes[13])
                .write(psSizes[14])
                .write(psSizes[15])
                .write(psSizes[16])
                .write(psSizes[17])
                .write(psSizes[18])
                .endl();
            }
            
            datafile _df;
        };
        
    } // datafiles
} // ealib

#endif
