#include <iostream>
#include "global_parameters.hpp"
#include "model_driver.hpp"
#include "testing.hpp"
#include <fstream>

std::string RUN_NAME = "default";
std::string FILE_PATH = "";
bool VERBOSE = true;

unsigned int RUN_TIME = 5000;
unsigned int OUTPUT_INTERVAL = 500;
unsigned int BURN_IN_PERIOD = 2000; //Recombination only allowed after this period.

unsigned int NUM_HOSTS = 5000;
unsigned int NUM_MOSQUITOES = 10000;
unsigned int INITIAL_ANTIGEN_DIVERSITY = 2*60; //Number of antigen variants sampled to provide initial antigen diversity.
unsigned int INITIAL_NUM_STRAINS = 4; //Number of strains which are constructed from the initial antigen pool at the start of simulation. Used to initialise initial infections.
unsigned int INITIAL_NUM_MOSQUITO_INFECTIONS = 1000;

unsigned int NUM_PHENOTYPES = 50000; //Number of antigenically distinct phenotypes for antigens.
unsigned int NUM_GENOTYPE_ONLY_BITS = 7; //The first x bits encode genotype only variation - no phenotypic effects.
unsigned int GENOTYPIC_SPACE_SIZE = std::numeric_limits<unsigned int>::max(); //Total genotypic space across all phenotypes.

unsigned int REPERTOIRE_SIZE = 60;

unsigned int MEAN_MOSQUITO_LIFE_EXPECTANCY = 32; //In days, Bellan2010.
unsigned int MOSQUITO_EIP = 10; //Extrinsic inoculation period in days, Deitz1974.
float BITE_RATE = 0.2f; //Daily probability of mosquito biting.

float INTERGENIC_RECOMBINATION_P = 0.1f; //Per gene intergenic recombination probability (swap genes).
float INTRAGENIC_RECOMBINATION_P = 0.2f; //Per gene intragenic recombination probability (hybrid gene).
std::array<float, 2> RECOMBINATION_CUMU_P = {INTERGENIC_RECOMBINATION_P, INTERGENIC_RECOMBINATION_P+INTRAGENIC_RECOMBINATION_P};
float RECOMBINATION_SCALE = 600.0f; //Scales the degree of genotypic change resulting from intragenic recombination

float INFECTION_DURATION_SCALE = 1.0;
float INFECTIVITY_SCALE = 0.0001;//0.25;


bool set_globals_from_cmd(int argc, char* argv[]);
bool parse_cmd_argument(const std::string& token, const std::string& value);
void recalc_derived_parameters(); //Some parameters aren't set directly but depend on others. After setting parameters these must be recalculated.

int main(int argc, char* argv[])
{
    if (!set_globals_from_cmd(argc, argv))
        return -1;

    ModelDriver model;
    model.run_model();

    //testing::run_tests();


    return 0;
}

//http://stackoverflow.com/questions/5290089/how-to-convert-a-number-to-string-and-vice-versa-in-c
bool parse_cmd_argument(const std::string& token, const std::string& value)
{
    if (token == "run_name")
        RUN_NAME = value;
    else if (token == "file_path")
    {
        if (value.back() != '/') //File path must end with a '/'...
            FILE_PATH = value+"/";
        else
            FILE_PATH = value;
    }
    else if (token == "run_time")
        RUN_TIME = stoi(value);
    else if (token == "output_interval")
        OUTPUT_INTERVAL = stoi(value);
    else if (token == "burn_in_period")
        BURN_IN_PERIOD = stoi(value);
    else if (token == "num_hosts")
        NUM_HOSTS = stoi(value);
    else if (token == "num_mosquitoes")
        NUM_MOSQUITOES = stoi(value);
    else if (token == "initial_antigen_diversity")
        INITIAL_ANTIGEN_DIVERSITY = stoi(value);
    else if (token == "initial_num_strains")
        INITIAL_NUM_STRAINS = stoi(value);
    else if (token == "initial_num_mosquito_infections")
        INITIAL_NUM_MOSQUITO_INFECTIONS = stoi(value);
    else if (token == "num_phenotypes")
        NUM_PHENOTYPES = stoi(value);
    else if (token == "num_genotype_only_bits")
        NUM_GENOTYPE_ONLY_BITS = stoi(value);
    else if (token == "repertoire_size")
        REPERTOIRE_SIZE = stoi(value);
    else if (token == "mean_mosquito_life_expectancy")
        MEAN_MOSQUITO_LIFE_EXPECTANCY = stoi(value);
    else if (token == "mosquito_eip")
        MOSQUITO_EIP = stoi(value);
    else if (token == "bite_rate")
        BITE_RATE = stof(value);
    else if (token == "intergenic_recombination_p")
        INTERGENIC_RECOMBINATION_P = stof(value);
    else if (token == "intragenic_recombination_p")
        INTRAGENIC_RECOMBINATION_P = stof(value);
    else if (token == "recombination_scale")
        RECOMBINATION_SCALE = stof(value);
    else if (token == "infection_duration_scale")
        INFECTION_DURATION_SCALE = stof(value);
    else if (token == "infectivity_scale")
        INFECTIVITY_SCALE = stof(value);
    else if (token == "verbose")
    {
        if (value == "true")
            VERBOSE = true;
        else
            VERBOSE = false;
    }
    else
    {
        std::cout << "\n\nUnrecognised parameter token!\n";
        return false;
    }

    return true;
}

bool set_globals_from_cmd(int argc, char* argv[])
{
    //Check number of cmdline arguments
    if (argc == 1) //One argument is fine as this is the program being called.
        return true;
    else if ((argc-1) % 2 != 0) //Must be an odd number of arguments, as it is the first argument plus a number of parameter:value pairs.
    {
        std::cout << "\n\nMismatched number of command line arguments. Cannot parse token:value pairs.\n";
        std::cout << argc << " input arguments detected:\n";
        for (int i=0; i<argc; ++i)
            std::cout << "\t" << argv[i] << "\n";
        return false;
    }

    //Parse token:value pairs.
    for (int i=1; i<argc; i+=2)
    {
        std::string token(argv[i]);
        std::string value(argv[i+1]);
        if (parse_cmd_argument(token, value) == false)
        {
            std::cout << "unrecognised token: " << token << ", value: " << value << "\n";
            return false;
        }
    }
    recalc_derived_parameters(); //Some parameters are derived from others, now user parameters have changed these need calculating.

    //Write called arguments to file
    std::ofstream file;
    file.open(FILE_PATH+RUN_NAME+"_calling_arguments.txt", std::ofstream::out | std::ofstream::trunc);
    file << argv[0] << "\n\n";
    for (int i=1; i<argc; i+=2)
        file << argv[i] << " " << argv[i+1] << "\n";
    file.flush();
    file.close();

    return true;
}

//Some parameters aren't set directly but depend on others. After setting parameters these must be recalculated.
void recalc_derived_parameters()
{
    RECOMBINATION_CUMU_P = {INTERGENIC_RECOMBINATION_P, INTERGENIC_RECOMBINATION_P+INTRAGENIC_RECOMBINATION_P};
    GENOTYPIC_SPACE_SIZE = std::numeric_limits<unsigned int>::max(); //Total genotypic space across all phenotypes.
}
