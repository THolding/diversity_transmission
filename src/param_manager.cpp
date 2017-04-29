#include "param_manager.hpp"
#include "mosquito.hpp"
#include "utilities.hpp"
#include "adaptors/output_interval_adaptor.hpp"
#include "diversity_monitor.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

#include <iostream>
#include <omp.h>

std::string ParamManager::runName = "default";
std::string ParamManager::filePath = "";

bool ParamManager::verbose = false;
bool ParamManager::unique_initial_strains = false;

unsigned int ParamManager::run_time = 10000;//50000;
unsigned int ParamManager::output_interval = 250;//250;
unsigned int ParamManager::burn_in_period = 3000;//3000;
unsigned int ParamManager::num_hosts = 8000;
unsigned int ParamManager::initial_num_mosquitoes = 8000;
unsigned int ParamManager::initial_antigen_diversity = 2500;
unsigned int ParamManager::initial_num_strains = 40;
unsigned int ParamManager::initial_num_mosquito_infections = 500;
unsigned int ParamManager::num_phenotypes = 50000;
unsigned int ParamManager::num_genotype_only_bits = 7;
unsigned int ParamManager::genotypic_space_size = std::numeric_limits<unsigned int>::max();
unsigned int ParamManager::repertoire_size = 60;
unsigned int ParamManager::mean_mosquito_life_expectancy = 32; //In days, Bellan2010.
unsigned int ParamManager::mosquito_eip = 10; //Extrinsic inoculation period in days, Deitz1974.

float ParamManager::bite_rate = 0.12f;
float ParamManager::intergenic_recombination_p = 0.01f;
float ParamManager::intragenic_recombination_p = 0.002f;
float ParamManager::recombination_scale = 50.0f;
float ParamManager::infection_duration_scale = 2.0f;
float ParamManager::infectivity_scale = 0.5f;
float ParamManager::cross_immunity = 0.0f;
float ParamManager::immunityScale = 1.0f;

////Output management
bool ParamManager::output_antigen_frequency = false; //Outputs the frequency with which antigens are present in the parasite population.
bool ParamManager::output_host_susceptibility = false; //Outputs a number (ranging between 0 and 1) indicating the mean susceptibility of the host popualtion to currently circulating parasite population.
bool ParamManager::output_strain_structure = false; //Output a list of all strain vector frequencies each output interval (uses multiple files).

////Dynamic support parameters.
//Dynamic mosquito population (MosquitoPopulationAdaptor).
bool ParamManager::dyn_num_mosquitoes = false; //Used to know whether or not to output timeseries of number of mosquitoes for example.
unsigned int ParamManager::max_num_mosquitoes = initial_num_mosquitoes; //Used to reserve space in vectors.
//Dynamic bite rate (BiteRateAdaptor).
bool ParamManager::dyn_bite_rate = false;
bool ParamManager::dyn_intragenic_recombination_p = false;

std::list<float> ParamManager::immunityMask;
BITE_FREQUENCY_TABLE ParamManager::cumulativeBiteFrequencyDistribution;
unsigned int ParamManager::output_size_needed = 0;
std::array<float, 2> ParamManager::recombination_cumu_p;
std::list<Adaptor*> ParamManager::adaptors;

bool ParamManager::recalculate_derived_parameters()
{
    if (!recalculate_recombination_distributions())
        throw std::runtime_error("ParamManager::recalculate_derived_parameters failed to calculate recombination distributions");

    //if (!next_function())
        //throw std::runtime_error("");

    recalculate_cumulative_bite_frequency_distribution();
    recalculate_output_array_size_needed();
    recalculate_immunity_mask();

    //if (paramsBool["output_parasite_adaptedness"])
    //    paramsBool["output_antigen_frequency"] = true;

    DiversityMonitor::reset();

    return true;
}

bool ParamManager::recalculate_recombination_distributions()
{
     recombination_cumu_p = { intergenic_recombination_p, intergenic_recombination_p+intragenic_recombination_p};
     return true;
}

void ParamManager::recalculate_cumulative_bite_frequency_distribution()
{
    const float biteRate = ParamManager::bite_rate;
    BITE_FREQUENCY_TABLE pdfPoisson;
    for (unsigned int k=0; k<pdfPoisson.size(); k++)
        pdfPoisson[k] = ((float)std::pow(biteRate, k) * std::exp(-biteRate)) / (float) utilities::factorial(k);

    //Calculate cumulative of pdfPoisson.
    cumulativeBiteFrequencyDistribution[0] = pdfPoisson[0];
    for (unsigned int i=1; i<pdfPoisson.size(); i++)
        cumulativeBiteFrequencyDistribution[i] = cumulativeBiteFrequencyDistribution[i-1]+pdfPoisson[i];
}

void ParamManager::recalculate_immunity_mask()
{
    immunityMask.clear();

    //Generate one half of a normal distribution with mean of 0 and variance 'cross_immunity'.

    immunityMask.push_back(1.0); //peak.

    if (ParamManager::cross_immunity == 0.0f) //If cross-immunity is turned off, just return the peak.
        return;
    else
    {
        //Calculate one tail.
        unsigned int x = 1;
        float curVal = 1.0;
        float mu = 0.0;
        float var = ParamManager::cross_immunity; //Variance
        float amp = 1.0;

        while (curVal >= 0.01) //Keep going until the distribution is very low.
        {
            curVal = amp * std::exp( -(std::pow(x-mu,2) / (2*std::pow(var,2))) );
            immunityMask.push_back(curVal);
            ++x;
        }

        //Copy right hand tail to the left.
        unsigned int tailSize = immunityMask.size()-1;
        auto itr = immunityMask.begin();
        itr++;
        for (unsigned int i=0; i<tailSize; i++)
        {
            immunityMask.insert(immunityMask.begin(), *itr);
            itr++;
        }
    }
}

void ParamManager::recalculate_output_array_size_needed()
{
    unsigned int sizeNeeded = (run_time / output_interval) + 1;

    for (const Adaptor* adaptor : adaptors)
    {
        //For each OutputIntervalAdaptor adjust the baseline size needed by the appropriate amount.
        if (adaptor->get_adaptor_name() == "OutputIntervalAdaptor")
        {
            //std::cout << "before: " << output_size_needed << "\n";
            unsigned int adaptorSizeNeeded = (adaptor->get_stop_t() - adaptor->get_start_t()) / static_cast<const OutputIntervalAdaptor*>(adaptor)->get_target_output_interval();
            unsigned int baselineSizeNeeded = (adaptor->get_stop_t() - adaptor->get_start_t()) / output_interval;
            int change = adaptorSizeNeeded - baselineSizeNeeded;
            sizeNeeded += change;
            //std::cout << "before: " << sizeNeeded << "\n";
        }
    }

    output_size_needed = sizeNeeded;
}

void ParamManager::set_param(const std::string name, const std::string value)
{
    if (name == "run_name")
        runName = value;
    else if (name == "file_path")
        filePath = value;

    else if (name == "verbose")
        verbose = (value == "true" || value == "1" || value == "True" || value == "TRUE");
    else if (name == "unique_initial_strains")
        unique_initial_strains = (value == "true" || value == "1" || value == "True" || value == "TRUE");

    else if (name == "run_time")
        run_time = std::stoi(value);
    else if (name == "output_interval")
        output_interval = std::stoi(value);
    else if (name == "burn_in_period")
        burn_in_period = std::stoi(value);
    else if (name == "num_hosts")
        num_hosts = std::stoi(value);
    else if (name == "initial_num_mosquitoes")
        initial_num_mosquitoes = std::stoi(value);
    else if (name == "initial_antigen_diversity")
        initial_antigen_diversity = std::stoi(value);
    else if (name == "initial_num_strains")
        initial_num_strains = std::stoi(value);
    else if (name == "initial_num_mosquito_infections")
        initial_num_mosquito_infections = std::stoi(value);
    else if (name == "num_phenotypes")
        num_phenotypes = std::stoi(value);
    else if (name == "num_genotype_only_bits")
        num_genotype_only_bits = std::stoi(value);
    else if (name == "genotypic_space_size")
        genotypic_space_size = std::stoi(value);
    else if (name == "repertoire_size")
        repertoire_size = std::stoi(value);
    else if (name == "mean_mosquito_life_expectancy")
        mean_mosquito_life_expectancy = std::stoi(value);
    else if (name == "mosquito_eip")
        mosquito_eip = std::stoi(value);

    else if (name == "bite_rate")
        bite_rate = std::stof(value);
    else if (name == "intergenic_recombination_p")
        intergenic_recombination_p = std::stof(value);
    else if (name == "intragenic_recombination_p")
        intragenic_recombination_p = std::stof(value);
    else if (name == "recombination_scale")
        recombination_scale = std::stof(value);
    else if (name == "infection_duration_scale")
        infection_duration_scale = std::stof(value);
    else if (name == "infectivity_scale")
        infectivity_scale = std::stof(value);
    else if (name == "cross_immunity")
        cross_immunity = std::stof(value);
    else if (name == "immunity_scale")
        immunityScale = std::stof(value);

    else if (name == "output_antigen_frequency")
        output_antigen_frequency = (value == "true" || value == "1" || value == "True" || value == "TRUE");
    else if (name == "output_host_susceptibility")
        output_host_susceptibility = (value == "true" || value == "1" || value == "True" || value == "TRUE");
    else if (name == "output_strain_structure")
        output_strain_structure = (value == "true" || value == "1" || value == "True" || value == "TRUE");

    else if (name == "dyn_num_mosquitoes")
        dyn_num_mosquitoes = (value == "true" || value == "1" || value == "True" || value == "TRUE");
    else if (name == "max_num_mosquitoes")
        max_num_mosquitoes = std::stoi(value);
    else if (name == "dyn_bite_rate")
        dyn_bite_rate = (value == "true" || value == "1" || value == "True" || value == "TRUE");
    else if (name == "dyn_intragenic_recombination_p")
        dyn_intragenic_recombination_p = (value == "true" || value == "1" || value == "True" || value == "TRUE");

    else
        throw std::runtime_error("ParamManager::set_param cannot create a new parameter with name '" + name + "'. You can only set values for pre-existing parameters.");
}

void ParamManager::add_adaptor(Adaptor* adaptor)
{
    if (is_compatable_adaptor(*adaptor))
        adaptors.push_back(adaptor);
    else
        throw std::runtime_error("ParamManager::add_adaptor: Cannot add adaptor (adaptor name: '"+adaptor->get_adaptor_name()+"') which starts at t="+std::to_string(adaptor->get_start_t())+" and ends at t="+std::to_string(adaptor->get_stop_t())+" because it is not compatable with adaptors already added.");
}

//Returns true if the proposed adaptor is compatable with current list of adaptors.
bool ParamManager::is_compatable_adaptor(const Adaptor& proposed)
{
    //Check start and stop times make sense.
    if (proposed.get_stop_t() < proposed.get_start_t())
        throw std::runtime_error("ParamManager::is_compatable_adaptor: Cannot add adaptor with an earlier stop time than start time!");

    //iterate through each adaptor and see if they match parameter and overlap start/stop interval. Return false if so.
    for (const Adaptor* adaptor : adaptors) {
        if (adaptor->get_adaptor_name() == proposed.get_adaptor_name())
        {
            //Overlapping periods of effect are not allowed for the same parameter.
            //If the start overlaps another's period.
            if (proposed.get_start_t() >= adaptor->get_start_t() && proposed.get_start_t() < adaptor->get_stop_t())
                return false;
            //If the stop overlaps another's period
            else if (proposed.get_stop_t() > adaptor->get_start_t() && proposed.get_stop_t() <= adaptor->get_stop_t())
                return false;
            //If proposed range encompasses another range
            else if (proposed.get_start_t() <= adaptor->get_start_t() && proposed.get_stop_t() >= adaptor->get_stop_t())
                return false;
        }
    }
    return true;
}

void ParamManager::update_adaptors(const unsigned int t)
{
    for (Adaptor* adaptor : adaptors)
        adaptor->update(t);
}


ParamManager::~ParamManager()
{
    //Delete all our adaptors. This leaves our list of adaptors full of invalid pointers - so we must then clear it!
    for (Adaptor* adaptor : adaptors)
        delete adaptor;
    adaptors.clear();
}
