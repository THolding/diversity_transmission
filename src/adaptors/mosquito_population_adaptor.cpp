#include "mosquito_population_adaptor.hpp"
#include "../mosquito_manager.hpp"
#include <cmath>
#include <iostream>

MosquitoPopulationAdaptor::MosquitoPopulationAdaptor(unsigned int tStart, unsigned int tStop, unsigned int targetPopulation, MosquitoManager* mManager)
    : tStart(tStart), tStop(tStop), adaptorName("MosquitoPopulationAdaptor"), targetPopulation(targetPopulation), mManager(mManager)
{
    ParamManager::instance().set_bool("dyn_num_mosquitoes", true);

    if (targetPopulation > ParamManager::instance().get_int("max_num_mosquitoes"))
        ParamManager::instance().set_int("max_num_mosquitoes", targetPopulation);
}


void MosquitoPopulationAdaptor::update(unsigned int time)
{
    //std::cout << "***Num mosquitoes (current): " << mManager->get_count() << "\n";
    if (time >= tStart && time < tStop)
    {
        int numLeftToChange = targetPopulation - mManager->get_count();
        unsigned int timeLeft = tStop - time;
        //std::cout << "\tnumLeftToChange: " << numLeftToChange << "\n";
        //std::cout << "\ttimeLeft: " << timeLeft << "\n";
        fractionalChange += (float)numLeftToChange / (float)timeLeft; //How much should the mosquito population be changed by?
        int integerToChange = 0;
        if (fractionalChange > 0)
            integerToChange = std::floor(fractionalChange);
        else if (fractionalChange < 0)
            integerToChange = std::ceil(fractionalChange);
        fractionalChange -= integerToChange; //We are about to modify mosquito population by this number so we don't need to track it anymore
        mManager->modify_population(integerToChange); //Modify mosquito population.
    }
}
