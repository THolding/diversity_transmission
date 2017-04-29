#include "mosquito_manager.hpp"
#include "utilities.hpp"

void MosquitoManager::initialise(std::vector<Mosquito>* mosquitoesArray)
{
    numMosquitoes = 0;
    mosquitoes = mosquitoesArray;
    activeMosquitoes.reserve(ParamManager::max_num_mosquitoes);
    inactiveMosquitoes.reserve(ParamManager::max_num_mosquitoes);
    for (unsigned int i=0; i<mosquitoes->size(); ++i)
    {
        if (mosquitoes->at(i).is_active())
        {
            activeMosquitoes.push_back(i);
            numMosquitoes++;
        }
        else
            inactiveMosquitoes.push_back(i);
    }
}

void MosquitoManager::remove_mosquito(unsigned int numToRemove)
{
    for (unsigned int i=0; i<numToRemove; ++i)
    {
        if (!activeMosquitoes.empty())
        {
            unsigned int iM = activeMosquitoes.back();
            mosquitoes->at(iM).kill();
            mosquitoes->at(iM).active = false;
            activeMosquitoes.pop_back();
            inactiveMosquitoes.push_back(iM);
            --numMosquitoes;
        }
    }
}
void MosquitoManager::add_mosquito(unsigned int numToAdd)
{
    for (unsigned int i=0; i<numToAdd; ++i)
    {
        if (inactiveMosquitoes.empty()) //then we need to create a new mosquitoe...
        {
            Mosquito newMosquito;
            newMosquito.kill();
            newMosquito.active = true;
            mosquitoes->push_back(newMosquito);
            activeMosquitoes.push_back(mosquitoes->size()-1);
            numMosquitoes++;
        }
        else //inactive mosquito can be reactivated
        {
            unsigned int iM = inactiveMosquitoes.back();
            mosquitoes->at(iM).kill();
            mosquitoes->at(iM).active = true;
            inactiveMosquitoes.pop_back();
            activeMosquitoes.push_back(iM);
            numMosquitoes++;
        }
    }
}

//Just selects remove_mosquito / add_mosquito as appropriate.
void MosquitoManager::modify_population(int numToChange)
{
    if (numToChange > 0)
        add_mosquito(numToChange);
    else if (numToChange < 0)
        remove_mosquito(std::abs(numToChange));
    else
        return;
}

unsigned int MosquitoManager::random_active_mos() const
{
    return activeMosquitoes[utilities::random(0, activeMosquitoes.size())];
}
