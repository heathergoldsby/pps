#include <ea/digital_evolution.h>
#include <ea/cmdline_interface.h>
#include <ea/subpopulation_founder.h>
#include <ea/line_of_descent.h>
#include <ea/analysis/archive.h>
#include <ea/generational_models/periodic_competition.h>
//#include <ea/generational_models/moran_process.h>
#include <ea/selection/rank.h>
#include <ea/selection/proportionate.h>

#include <ea/datafiles/fitness.h>
#include <ea/digital_evolution/extra_instruction_sets/matrix.h>

#include "evo_propagule_ins.h"
#include "multibirth_not_nand_prop_ancestor.h"

#include "subpopulation_propagule_split.h"

#include "stripes.h"
#include "meta_moran_process.h"
#include "propagule.h"
#include "movie.h"
//#include "knockouts.h"




using namespace ealib;



//! Configuration object for an EA.
struct lifecycle : public default_lifecycle {
    
    //! Called as the final step of EA construction (must not depend on configuration parameters)
    template <typename EA>
    void after_initialization(EA& ea) {
        if(ea.isa().size()) {
            return;
        }
        
        using namespace instructions;
        append_isa<nop_a>(0,ea);
        append_isa<nop_b>(0,ea);
        append_isa<nop_c>(0,ea);
        append_isa<nop_x>(ea);
        append_isa<mov_head>(ea);
        append_isa<if_label>(ea);
        append_isa<h_search>(ea);
        append_isa<nand>(ea);
        append_isa<push>(ea);
        append_isa<pop>(ea);
        append_isa<swap>(ea);
        append_isa<inc>(ea);
        append_isa<dec>(ea);
        append_isa<tx_msg_check_task_matrix>(ea);
        append_isa<tx_msg_matrix>(ea);
        append_isa<rx_msg>(ea);
        append_isa<bc_msg_matrix>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea);
        append_isa<h_alloc>(ea);
        append_isa<h_copy>(ea);
        append_isa<h_divide_soft_parent_reset>(ea);
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<is_neighbor_matrix>(ea);

        append_isa<deploy_propagule>(ea);
        append_isa<prop_size_1>(ea);
        append_isa<prop_size_2>(ea);
        append_isa<prop_size_3>(ea);
        append_isa<prop_size_4>(ea);
        append_isa<prop_size_5>(ea);
        append_isa<prop_size_6>(ea);
        append_isa<prop_size_7>(ea);
        append_isa<prop_size_8>(ea);
        append_isa<prop_size_9>(ea);
        append_isa<prop_size_10>(ea);
        append_isa<prop_size_11>(ea);
        append_isa<prop_size_12>(ea);
        append_isa<prop_size_13>(ea);
        append_isa<prop_size_14>(ea);
        append_isa<prop_size_15>(ea);
        append_isa<prop_size_16>(ea);
        append_isa<prop_size_17>(ea);
        append_isa<prop_size_18>(ea);
        
        
        
        add_event<task_resource_consumption>(ea);
        add_event<task_switching_cost>(ea);
        //        add_event<prop_death_event>(ea);
        
        add_event<ts_birth_event>(ea);
        
        typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
        typedef typename EA::resource_ptr_type resource_ptr_type;
        
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        
        resource_ptr_type resA = make_resource("resA", ea);
        resource_ptr_type resB = make_resource("resB", ea);
        
        task_not->consumes(resA);
        task_nand->consumes(resB);
        
        
    }
    
};

template <typename T>
struct subpop_trait : subpopulation_founder_trait<T>, fitness_trait<T> {
    typedef subpopulation_founder_trait<T> parent1_type;
    typedef fitness_trait<T> parent2_type;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp("subpopulation_founder_trait", boost::serialization::base_object<parent1_type>(*this));
        ar & boost::serialization::make_nvp("fitness_trait", boost::serialization::base_object<parent2_type>(*this));
    }
};




typedef digital_evolution
< lifecycle
, recombination::asexual
, round_robin
, multibirth_not_nand_prop_ancestor
, empty_facing_neighbor_matrix
, dont_stop
, generate_single_ancestor
> sea_type;

// original , generational_models::periodic_competition < generational_models::meta_moran_process< selection::random< >, selection::rank< > >, generational_models::isolated_subpopulations > // generational_models::moran_process< >, isolated_subpopulations


// flipped to avoid age: // , generational_models::periodic_competition < generational_models::meta_moran_process< selection::tournament< >, selection::random < > >, generational_models::isolated_subpopulations >
//, generational_models::periodic_competition < generational_models::meta_moran_process< selection::random< >, selection::proportionate< > >, generational_models::isolated_subpopulations > // generational_models::moran_process< >, isolated_subpopulations


typedef metapopulation
< sea_type
, permute_stripes
, mutation::operators::no_mutation
, subpopulation_propagule_split
, generational_models::periodic_competition < generational_models::meta_moran_process< selection::random< >, selection::proportionate< > >, generational_models::isolated_subpopulations > // generational_models::moran_process< >, isolated_subpopulations
, ancestors::default_subpopulation
, dont_stop
, fill_metapopulation
, default_lifecycle
, subpop_trait
> mea_type;




/*!
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        add_option<SPATIAL_X>(this);
        add_option<SPATIAL_Y>(this);
        add_option<METAPOPULATION_SIZE>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<SCHEDULER_RESOURCE_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<MORAN_REPLACEMENT_RATE_P>(this);
        
        add_option<ANALYSIS_INPUT>(this);

        
        // ts specific options
        add_option<TASK_SWITCHING_COST>(this);

        
        add_option<TOURNAMENT_SELECTION_K>(this);
        add_option<TOURNAMENT_SELECTION_N>(this);
        
        add_option<METAPOP_COMPETITION_PERIOD>(this);
        add_option<PROPAGULE_SIZE>(this);

        
        
    }
    
    virtual void gather_tools() {
        add_tool<ealib::analysis::movie_for_competitions>(this);

        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<subpopulation_founder_event>(ea);
        add_event<datafiles::fitness_dat>(ea);
        add_event<datafiles::propagule_dat>(ea);

        
    }
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);