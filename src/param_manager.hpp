#pragma once
#include "adaptors/adaptor.hpp"
#include "global_typedefs.hpp"
#include <array>
#include <list>
#include <string>
#include <unordered_map>

class Adaptor;

//Based on http://stackoverflow.com/questions/1008019/c-singleton-design-pattern
class ParamManager
{
private:
    //
public:
    static std::string runName;
    static std::string filePath;

    static bool verbose;
    static bool unique_initial_strains;

    static unsigned int run_time;
    static unsigned int output_interval;
    static unsigned int burn_in_period;
    static unsigned int num_hosts;
    static unsigned int initial_num_mosquitoes;
    static unsigned int initial_antigen_diversity;
    static unsigned int initial_num_strains;
    static unsigned int initial_num_mosquito_infections;
    static unsigned int num_phenotypes;
    static unsigned int num_genotype_only_bits;
    static unsigned int genotypic_space_size;
    static unsigned int repertoire_size;
    static unsigned int mean_mosquito_life_expectancy; //In days, Bellan2010.
    static unsigned int mosquito_eip; //Extrinsic inoculation period in days, Deitz1974.

    static float bite_rate;
    static float intergenic_recombination_p;
    static float intragenic_recombination_p;
    static float recombination_scale;
    static float infection_duration_scale;
    static float infectivity_scale;
    static float cross_immunity;
    static float immunityScale; //Linear scaling of immunity.

    ////Output management
    static bool output_antigen_frequency; //Outputs the frequency with which antigens are present in the parasite population.
    static bool output_host_susceptibility; //Outputs a number (ranging between 0 and 1) indicating the mean susceptibility of the host popualtion to currently circulating parasite population.
    static bool output_strain_structure; //Output a list of all strain vector frequencies each output interval (uses multiple files).

    ////Dynamic support parameters.
    //Dynamic mosquito population (MosquitoPopulationAdaptor).
    static bool dyn_num_mosquitoes; //Used to know whether or not to output timeseries of number of mosquitoes for example.
    static unsigned int max_num_mosquitoes; //Used to reserve space in vectors.
    //Dynamic bite rate (BiteRateAdaptor).
    static bool dyn_bite_rate;
    static bool dyn_intragenic_recombination_p;

    static unsigned int output_size_needed; //Calculated as a derived parameter based on any output_interval adaptors

    static std::array<float, 2> recombination_cumu_p;

    static BITE_FREQUENCY_TABLE cumulativeBiteFrequencyDistribution;
    static std::list<float> immunityMask;

    static std::list<Adaptor*> adaptors;

    static std::string run_name() { return runName; }
    static std::string file_path() { return filePath; }

    //Functions
    static void set_param(const std::string name, const std::string value);

    static void add_adaptor(Adaptor* adaptor);
    static bool is_compatable_adaptor(const Adaptor& adaptor);

    static void update_adaptors(const unsigned int t);

    static bool recalculate_derived_parameters();
    static bool recalculate_recombination_distributions();
    static void recalculate_cumulative_bite_frequency_distribution();
    static const BITE_FREQUENCY_TABLE& get_cumulative_bite_frequency_distribution() { return cumulativeBiteFrequencyDistribution; }
    static void recalculate_output_array_size_needed();
    static void recalculate_immunity_mask();
    static const std::list<float>& get_immunity_mask() { return immunityMask; }

    ParamManager(ParamManager const&) = delete; //disable copy construction
    void operator=(ParamManager const&) = delete; //disable copy assignment

    ~ParamManager();
};
