#include "infection.hpp"
#include "param_manager.hpp"
#include "utilities.hpp"
#include "diversity_monitor.hpp"

#include <iostream>

void Infection::reset()
{
    if (infected) //register loss of antigen abundance
    {
        DiversityMonitor::register_lost_strain(strain);
    }

    infected = false;
    durationRemaining = 0;
    infectivity = 0.0f;
}

float infectivity_kernal(const Strain& strain, const ImmuneState& immuneState)
{
    return 1.0*ParamManager::infectivity_scale;
}


unsigned short duration_kernal(const Strain& strain, const ImmuneState& immuneState)
{
    float duration = 0.0;
    for (const Antigen antigen : strain)
    {
        duration += ParamManager::infection_duration_scale * (1.0-immuneState.at(get_phenotype_id(antigen)));
        //immuneState[get_phenotype_id(antigen)] = 1.0; //Temp - stops multiple expression
        //std::cout << "\t\tDurationCalc: " << duration << "\timmuneStata[x]: " << immuneState[get_phenotype_id(antigen)] << "\n";
    }

    return int(duration);
}


void exposure_kernal(const Strain& strain, ImmuneState& immuneState)
{
    const std::list<float>& immunityMask = ParamManager::get_immunity_mask();
    //std::cout << immunityMask.size();
    unsigned int tailSize = (immunityMask.size()-1)/2;

    for (const Antigen parasiteAntigen : strain)
    {
        unsigned int targetAntigenID = get_phenotype_id(parasiteAntigen);
        auto itr = immunityMask.begin();
        unsigned int curAntigen = utilities::wrap((int)targetAntigenID-tailSize, 0, ParamManager::num_phenotypes);
        while (itr != immunityMask.end())
        {
            immuneState[curAntigen] += ParamManager::immunityScale*(*itr);
            if (immuneState[curAntigen] > 1.0)
                immuneState[curAntigen] = 1.0;

            //increment counters
            curAntigen = utilities::wrap(curAntigen+1, 0, ParamManager::num_phenotypes);
            itr++;
        }
    }
}
