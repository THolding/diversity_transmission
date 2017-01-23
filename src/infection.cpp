#include "infection.hpp"
#include "param_manager.hpp"

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


unsigned short duration_kernal(const Strain& strain, ImmuneState& immuneState)
{
    float duration = 0.0;
    for (const Antigen antigen : strain)
    {
        duration += ParamManager::instance().get_float("infection_duration_scale") * (1.0-immuneState[get_phenotype_id(antigen)]);
        immuneState[get_phenotype_id(antigen)] = 1.0; //Temp - stops multiple expression
        //std::cout << "\t\tDurationCalc: " << duration << "\timmuneStata[x]: " << immuneState[get_phenotype_id(antigen)] << "\n";
    }

    return int(duration);
}


/*void exposure_kernal(const Strain& strain, ImmuneState& immuneState)
{
    for (const Antigen antigen : strain)
    {
        immuneState[get_phenotype_id(antigen)] = 1.0;
        //if (immuneState[get_phenotype_id(antigen)] < 1.0)
        //    immuneState[get_phenotype_id(antigen)] += 0.9;

        //if (immuneState[get_phenotype_id(antigen)] > 1.0)
        //    immuneState[get_phenotype_id(antigen)] = 1.0;
    }
}*/
