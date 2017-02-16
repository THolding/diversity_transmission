#include <iostream>
#include "param_manager.hpp"
#include "model_driver.hpp"
#include "testing.hpp"
#include "utilities.hpp"
#include <fstream>
#include <sstream>
#include <limits>
#include "adaptors/mosquito_population_adaptor.hpp"
#include "adaptors/bite_rate_adaptor.hpp"
#include "adaptors/output_interval_adaptor.hpp"

#include "host.hpp"

void parse_parameters_from_cmd(int argc, char* argv[], ModelDriver& model);
void parse_adaptor_from_cmd(const std::string& adaptorType, const std::string& argList, ModelDriver& model);

int main(int argc, char* argv[])
{
    ParamManager::instance(); //Initialise the instance!
    ModelDriver model; //Empty project so we can setup references / pointers.
    parse_parameters_from_cmd(argc, argv, model); //Throws exception if fails.
    ParamManager::instance().recalculate_derived_parameters();


    //ParamManager::instance().set_int("run_time", 3000);
    //ParamManager::instance().set_int("burn_in_period", 1000);
    //ParamManager::instance().set_int("output_interval", 150);
    //ParamManager::instance().set_bool("output_antigen_frequency", true);
    //ParamManager::instance().recalculate_derived_parameters();

    //ParamManager::instance().set_float("cross_immunity", 3);
    ////ParamManager::instance().set_int("num_phenotypes", 300);
    //ParamManager::instance().recalculate_derived_parameters();
    //auto mask = ParamManager::instance().get_immunity_mask();
    //utilities::arrayToFile(mask, "immunity_mask.csv");

    //Host host;
    //std::vector<unsigned int> strain1 = {100*128};
    //host.infect(strain1);
    //utilities::arrayToFile(host.immuneState, "test_immune_state1.csv");

    //std::vector<unsigned int> strain2 = {150*128};
    //host.infect(strain2);
    //utilities::arrayToFile(host.immuneState, "test_immune_state2.csv");

    //for (float f : ParamManager::instance().get_immunity_mask())
    //    std::cout << f << "\n";

    model.run_model();

    return 0;
}

void parse_parameters_from_cmd(int argc, char* argv[], ModelDriver& model)
{
    //Check number of cmdline arguments
    if ((argc-1) % 2 != 0) //Must be an odd number of arguments, as it is the first argument plus a number of parameter:value pairs.
        throw std::runtime_error("main.cpp, set_globals_from_cmd(): Mismatched number of command line arguments. Cannot parse token:value pairs.");

    //Parse token:value pairs.
    for (int i=1; i<argc; i+=2)
    {
        std::string token(argv[i]);
        std::string value(argv[i+1]);

        if (token.find("adaptor") != std::string::npos) //If it's an adaptor definition parse it as an adaptor
            parse_adaptor_from_cmd(token, value, model);
        else //...otherwise assume it is setting a parameter
            ParamManager::instance().set_param(token, value); //Throwes exception if token doesn't name an existing parameter.
    }

    if (ParamManager::instance().recalculate_derived_parameters() == false) //Some parameters are derived from others, now user parameters have changed these need calculating.
        std::runtime_error("ERROR: Cannot set derived parameters - this indicates inconsistent parameter states e.g. stop time before start time in dynamic parameters.");

    //Write called arguments to file
    std::ofstream file;
    file.open(ParamManager::instance().file_path()+ParamManager::instance().run_name()+"_calling_arguments.txt", std::ofstream::out | std::ofstream::trunc);
    file << argv[0] << "\n\n";
    for (int i=1; i<argc; i+=2)
        file << argv[i] << " " << argv[i+1] << "\n";
    file.flush();
    file.close();
}

void parse_adaptor_from_cmd(const std::string &adaptorType, const std::string &argList, ModelDriver& model)
{
    const char delim = '+';

    std::stringstream ss;
    ss.str(argList);
    std::string curItem;

    //Parse start time.
    std::getline(ss, curItem, delim);
    //std::cout << "tStart:" << curItem << std::endl;
    unsigned int tStart = std::stoi(curItem);

    //Parse stop time.
    std::getline(ss, curItem, delim);
    //std::cout << "tStop:" << curItem << std::endl;
    unsigned int tStop = std::stoi(curItem);

    //Load target value into curItem but don't make any assumptions about the type until we know the adaptor type.
    std::getline(ss, curItem, delim);
    //std::cout << "target:" << curItem << std::endl;

    //Create the adaptor based on adaptorType string.
    if (adaptorType.find("mosquito_population_adaptor") != std::string::npos) {
        unsigned int targetPopulation = std::stoi(curItem);
        ParamManager::instance().add_adaptor(new MosquitoPopulationAdaptor(tStart, tStop, targetPopulation, model.get_mos_manager()));
        std::cout << "Added mosquito_population_adaptor.\n";
    }
    else if (adaptorType.find("bite_rate_adaptor") != std::string::npos) {
        float targetBiteRate = std::stof(curItem);
        ParamManager::instance().add_adaptor(new BiteRateAdaptor(tStart, tStop, targetBiteRate));
        std::cout << "Added bite_rate_adaptor.\n";
    }
    else if (adaptorType.find("output_interval_adaptor") != std::string::npos) {
        unsigned int targetOutputInterval = std::stoi(curItem);
        ParamManager::instance().add_adaptor(new OutputIntervalAdaptor(tStart, tStop, targetOutputInterval));
        std::cout << "Added output_interval_adaptor.\n";
    }
}
