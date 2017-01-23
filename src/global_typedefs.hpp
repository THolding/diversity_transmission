#pragma once
#include <array>

//Type definitions.
typedef std::array<float, 500> PTABLE; //Holds e.g. daily probability of death or survival by age (days for mosquitoes, years for hosts.
typedef unsigned int Antigen;
#define Strain std::vector<Antigen>
#define ImmuneState std::vector<float>
#define BITE_FREQUENCY_TABLE std::array<float, 10>
