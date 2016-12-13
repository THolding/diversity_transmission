#include "tmp_simulations.hpp"
#include "model_driver.hpp"

void run_dyn_mosquitoes()
{
    RUN_NAME = "dyn_mos";
    FILE_PATH = "";
    VERBOSE = false;
    OUTPUT_INTERVAL = 500;

    DYN_NUM_MOSQUITOES = true;
    NUM_MOSQUITOES_MAX = 15000;
    BURN_IN_PERIOD = 10000; //Recombination only allowed after this period.
    START_MOS_INC = 20000; //time point in days
    STOP_MOS_INC = 25000; //time point in days
    START_MOS_DCR = 125000; //time point in days
    STOP_MOS_DCR = 130000; //time point in days
    RUN_TIME = 250000;

    NUM_HOSTS = 7500;
    NUM_MOSQUITOES = 7500;
    INITIAL_ANTIGEN_DIVERSITY = 1000; //Number of antigen variants sampled to provide initial antigen diversity.
    INITIAL_NUM_STRAINS = 16; //Number of strains which are constructed from the initial antigen pool at the start of simulation. Used to initialise initial infections.
    INITIAL_NUM_MOSQUITO_INFECTIONS = 500;

    NUM_PHENOTYPES = 50000; //Number of antigenically distinct phenotypes for antigens.
    NUM_GENOTYPE_ONLY_BITS = 7; //The first x bits encode genotype only variation - no phenotypic effects.

    REPERTOIRE_SIZE = 60;

    MEAN_MOSQUITO_LIFE_EXPECTANCY = 32; //In days, Bellan2010.
    MOSQUITO_EIP = 10; //Extrinsic inoculation period in days, Deitz1974.
    BITE_RATE = 0.1f; //Daily probability of mosquito biting.

    INTERGENIC_RECOMBINATION_P = 0.01f; //Per gene intergenic recombination probability (swap genes).
    INTRAGENIC_RECOMBINATION_P = 0.002f; //Per gene intragenic recombination probability (hybrid gene).
    RECOMBINATION_SCALE = 600.0f; //Scales the degree of genotypic change resulting from intragenic recombination

    INFECTION_DURATION_SCALE = 1.0;
    INFECTIVITY_SCALE = 0.2;


    ModelDriver model;
    model.run_model();
}
