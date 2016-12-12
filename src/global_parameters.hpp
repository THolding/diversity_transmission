#pragma once
#include <array>
#include <limits>

extern std::string RUN_NAME;
extern std::string FILE_PATH;
extern bool VERBOSE;

extern unsigned int RUN_TIME;
extern unsigned int OUTPUT_INTERVAL;
extern unsigned int BURN_IN_PERIOD; //Recombination only allowed after this period.

extern unsigned int NUM_HOSTS;
extern unsigned int NUM_MOSQUITOES;
extern unsigned int INITIAL_ANTIGEN_DIVERSITY; //Number of antigen variants sampled to provide initial antigen diversity.
extern unsigned int INITIAL_NUM_STRAINS; //Number of strains which are constructed from the initial antigen pool at the start of simulation. Used to initialise initial infections.
extern unsigned int INITIAL_NUM_MOSQUITO_INFECTIONS;

extern unsigned int NUM_PHENOTYPES; //Number of antigenically distinct phenotypes for antigens.
extern unsigned int NUM_GENOTYPE_ONLY_BITS; //The first x bits encode genotype only variation - no phenotypic effects.
extern unsigned int GENOTYPIC_SPACE_SIZE; //Total genotypic space across all phenotypes.

extern unsigned int REPERTOIRE_SIZE;

extern unsigned int MEAN_MOSQUITO_LIFE_EXPECTANCY; //In days, Bellan2010.
extern unsigned int MOSQUITO_EIP; //Extrinsic inoculation period in days, Deitz1974.
extern float BITE_RATE; //Daily probability of mosquito biting.

extern float INTERGENIC_RECOMBINATION_P; //Per gene intergenic recombination probability (swap genes).
extern float INTRAGENIC_RECOMBINATION_P; //Per gene intragenic recombination probability (hybrid gene).
extern std::array<float, 2> RECOMBINATION_CUMU_P;
extern float RECOMBINATION_SCALE; //Scales the degree of genotypic change resulting from intragenic recombination

extern float INFECTION_DURATION_SCALE;
extern float INFECTIVITY_SCALE;

//Type definitions.
typedef std::array<float, 500> PTABLE; //Holds e.g. daily probability of death or survival by age (days for mosquitoes, years for hosts.
typedef unsigned int Antigen;
#define Strain std::vector<Antigen>
#define ImmuneState std::vector<float>
#define BITE_FREQUENCY_TABLE std::array<float, 10>
