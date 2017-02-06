#include "infection.hpp"
#include "param_manager.hpp"
#include "utilities.hpp"

#include <iostream>

void Infection::reset()
{
    infected = false;
    durationRemaining = 0;
    infectivity = 0.0f;
}

float infectivity_kernal(const Strain& strain, const ImmuneState& immuneState)
{
    return 1.0*ParamManager::instance().get_float("infectivity_scale");
}


unsigned short duration_kernal(const Strain& strain, const ImmuneState& immuneState)
{
    float duration = 0.0;
    for (const Antigen antigen : strain)
    {
        duration += ParamManager::instance().get_float("infection_duration_scale") * (1.0-immuneState.at(get_phenotype_id(antigen)));
        //immuneState[get_phenotype_id(antigen)] = 1.0; //Temp - stops multiple expression
        //std::cout << "\t\tDurationCalc: " << duration << "\timmuneStata[x]: " << immuneState[get_phenotype_id(antigen)] << "\n";
    }

    return int(duration);
}


void exposure_kernal(const Strain& strain, ImmuneState& immuneState)
{
    const std::list<float>& immunityMask = ParamManager::instance().get_immunity_mask();
    unsigned int tailSize = (immunityMask.size()-1)/2;

    for (const Antigen parasiteAntigen : strain)
    {
        unsigned int targetAntigenID = get_phenotype_id(parasiteAntigen);
        auto itr = immunityMask.begin();
        unsigned int curAntigen = utilities::wrap((int)targetAntigenID-tailSize, 0, ParamManager::instance().get_int("num_phenotypes"));
        while (itr != immunityMask.end())
        {
            immuneState[curAntigen] += (*itr);
            if (immuneState[curAntigen] > 1.0)
                immuneState[curAntigen] = 1.0;

            //increment counters
            curAntigen = utilities::wrap(curAntigen+1, 0, ParamManager::instance().get_int("num_phenotypes"));
            itr++;
        }
    }
}
