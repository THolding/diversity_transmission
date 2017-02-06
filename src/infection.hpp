#pragma once
#include "strain.hpp"

class Infection
{
public:
    bool infected = false;
    unsigned int durationRemaining = 0; //Either latent period (in mosquitoes) or infection duration (in hosts)
    Strain strain;
    float infectivity = 0.0f;

    void reset();
    std::string to_string() const;
};

//Temporary - need to be replaced with function pointers and modularised.
float infectivity_kernal(const Strain& strain, const ImmuneState& immuneState);
unsigned short duration_kernal(const Strain& strain, const ImmuneState& immuneState);
void exposure_kernal(const Strain& strain, ImmuneState& immuneState);
