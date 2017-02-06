#include "param_manager.hpp"
#include "mosquito.hpp"
#include "utilities.hpp"
#include "adaptors/output_interval_adaptor.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

#include <iostream>

//Default parameters setup in private constructor.
ParamManager::ParamManager()
{
    std::cout << "Initialised ParamManager singleton instance.\n";

    paramsBool["verbose"] = false;
    paramsInt["run_time"] = 50000;//50000;
    paramsInt["output_interval"] = 250;//250;
    paramsInt["burn_in_period"] = 3000;//3000;

    paramsInt["num_hosts"] = 8000;
    paramsInt["initial_num_mosquitoes"] = 8000;

    paramsInt["initial_antigen_diversity"] = 2500;
    paramsInt["initial_num_strains"] = 40;
    paramsInt["initial_num_mosquito_infections"] = 500;
    paramsInt["num_phenotypes"] = 50000;
    paramsInt["num_genotype_only_bits"] = 7;
    paramsInt["genotypic_space_size"] = std::numeric_limits<unsigned int>::max();
    paramsInt["repertoire_size"] = 60;
    paramsInt["mean_mosquito_life_expectancy"] = 32; //In days, Bellan2010.
    paramsInt["mosquito_eip"] = 10; //Extrinsic inoculation period in days, Deitz1974.
    paramsFloat["bite_rate"] = 0.12f;
    paramsFloat["intergenic_recombination_p"] = 0.01f;
    paramsFloat["intragenic_recombination_p"] = 0.002f;
    paramsFloat["recombination_scale"] = 50.0f;
    paramsFloat["infection_duration_scale"] = 2.0f;
    paramsFloat["infectivity_scale"] = 0.5f;
    paramsFloat["cross_immunity"] = 0.0f;

    ////Dynamic support parameters.
    //Dynamic mosquito population (MosquitoPopulationAdaptor).
    paramsBool["dyn_num_mosquitoes"] = false; //Used to know whether or not to output timeseries of number of mosquitoes for example.
    paramsInt["max_num_mosquitoes"] = paramsInt["initial_num_mosquitoes"]; //Used to reserve space in vectors.
    //Dynamic bite rate (BiteRateAdaptor).
    paramsBool["dyn_bite_rate"] = false;

    //paramsBool["dyn_num_hosts"] = false;
    //paramsInt["num_hosts_max"] = 0; //Used to reserve space in vectors.

    recalculate_derived_parameters();
}

bool ParamManager::recalculate_derived_parameters()
{
    if (!recalculate_recombination_distributions())
        throw std::runtime_error("ParamManager::recalculate_derived_parameters failed to calculate recombination distributions");

    //if (!next_function())
        //throw std::runtime_error("");

    recalculate_cumulative_bite_frequency_distribution();
    recalculate_output_array_size_needed();
    recalculate_immunity_mask();

    return true;
}

bool ParamManager::recalculate_recombination_distributions()
{
     recombination_cumu_p = { paramsFloat["intergenic_recombination_p"], paramsFloat["intergenic_recombination_p"]+paramsFloat["intragenic_recombination_p"]};
     return true;
}

void ParamManager::recalculate_cumulative_bite_frequency_distribution()
{
    const float biteRate = paramsFloat["bite_rate"];
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

    if (paramsFloat["cross_immunity"] == 0) //If cross-immunity is turned off, just return the peak.
        return;
    else
    {
        //Calculate one tail.
        unsigned int x = 1;
        float curVal = 1.0;
        float mu = 0.0;
        float var = paramsFloat["cross_immunity"]; //Variance
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
    unsigned int sizeNeeded = (paramsInt["run_time"] / paramsInt["output_interval"]) + 1;

    for (const Adaptor* adaptor : adaptors)
    {
        //For each OutputIntervalAdaptor adjust the baseline size needed by the appropriate amount.
        if (adaptor->get_adaptor_name() == "OutputIntervalAdaptor")
        {
            std::cout << "before: " << paramsInt["output_size_needed"] << "\n";
            unsigned int adaptorSizeNeeded = (adaptor->get_stop_t() - adaptor->get_start_t()) / static_cast<const OutputIntervalAdaptor*>(adaptor)->get_target_output_interval();
            unsigned int baselineSizeNeeded = (adaptor->get_stop_t() - adaptor->get_start_t()) / paramsInt["output_interval"];
            int change = adaptorSizeNeeded - baselineSizeNeeded;
            sizeNeeded += change;
            std::cout << "before: " << sizeNeeded << "\n";
        }
    }

    paramsInt["output_size_needed"] = sizeNeeded;
}

//Getters
float ParamManager::get_float(std::string paramName) const
{
    if (paramsFloat.count(paramName) == 1)
        return paramsFloat.at(paramName);
    else
        throw std::runtime_error("ParamManager::get_float Attempted to access unknown parameter with key: "+paramName);
}

unsigned int ParamManager::get_int(std::string paramName) const
{
    if (paramsInt.count(paramName) == 1)
        return paramsInt.at(paramName);
    else
        throw std::runtime_error("ParamManager::get_int Attempted to access unknown parameter with key: "+paramName);
}

bool ParamManager::get_bool(std::string paramName) const
{
    if (paramsBool.count(paramName) == 1)
        return paramsBool.at(paramName);
    else
        throw std::runtime_error("ParamManager::get_bool Attempted to access unknown parameter with key: "+paramName);
}

void ParamManager::set_param(const std::string name, const std::string value)
{
    if (paramsInt.count(name) == 1)
        set_int(name, std::stoi(value));
    else if (paramsFloat.count(name) == 1)
        set_float(name, std::stof(value));
    else if (paramsBool.count(name) == 1) {
        bool boolVal = (value == "true" || value == "1" || value == "True" || value == "TRUE");
        set_bool(name, boolVal);
        }
    else if (name == "run_name")
        set_run_name(value);
    else if (name == "file_path")
        set_file_path(value);
    else
    {
        throw std::runtime_error("ParamManager::set_param cannot create a new parameter with name " + name + ". You can only set values for pre-existing parameters.");
    }
}


//Setters
void ParamManager::set_float(const std::string paramName, const float param)
{
    if (paramsFloat.count(paramName) == 1)
        paramsFloat[paramName] = param;
    else
        throw std::runtime_error("ParamManager::set_float Attempted to set unknown parameter with key: "+paramName);
}

void ParamManager::set_int(const std::string paramName, const int param)
{
    if (paramsInt.count(paramName) == 1)
        paramsInt[paramName] = param;
    else
        throw std::runtime_error("ParamManager::set_int Attempted to set unknown parameter with key: "+paramName);
}

void ParamManager::set_bool(const std::string paramName, const bool param)
{
    if (paramsBool.count(paramName) == 1)
        paramsBool[paramName] = param;
    else
        throw std::runtime_error("ParamManager::set_bool Attempted to set unknown parameter with key: "+paramName);
}

void ParamManager::add_adaptor(Adaptor* adaptor)
{
    if (is_compatable_adaptor(*adaptor))
        adaptors.push_back(adaptor);
    else
        throw std::runtime_error("ParamManager::add_adaptor: Cannot add adaptor (adaptor name: '"+adaptor->get_adaptor_name()+"') which starts at t="+std::to_string(adaptor->get_start_t())+" and ends at t="+std::to_string(adaptor->get_stop_t())+" because it is not compatable with adaptors already added.");
}

//Returns true if the proposed adaptor is compatable with current list of adaptors.
bool ParamManager::is_compatable_adaptor(const Adaptor& proposed) const
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
