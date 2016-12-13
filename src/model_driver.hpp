#pragma once
#include "global_parameters.hpp"
#include "host.hpp"
#include "output.hpp"
#include "mosquito.hpp"


class MosquitoManager
{
private:
    unsigned int numMosquitoes=0;
    std::vector<Mosquito>* mosquitoes;
    std::vector<unsigned int> inactiveMosquitoes;
    std::vector<unsigned int> activeMosquitoes;

public:
    void initialise(std::vector<Mosquito>* mosquitoesArray);
    unsigned int get_count() const { return numMosquitoes; }
    void remove_mosquito(unsigned int numToRemove = 1);
    void add_mosquito(unsigned int numToAdd = 1);
    unsigned int random_active_mos() const;

};


class ModelDriver
{
private:
    PTABLE pDeathHosts;
    PTABLE pDeathMosquitoes;
    PTABLE cdfHosts;
    PTABLE cdfMosquitoes;
    BITE_FREQUENCY_TABLE cumulativeBiteFrequencyDistribution;

    std::vector<Host> hosts;
    std::vector<Mosquito> mosquitoes;
    MosquitoManager mManager;
    float mosChangeRemainder = 0.0;

    Output output;
    int burnInPeriod;

    void age_hosts();
    void age_mosquitoes();
    void update_host_infections();
    void update_mosquito_infections();
    void feed_mosquitoes();
    void update_parameters(const unsigned int time);
public:
    ModelDriver() : output(Output(this)) {  }
    void initialise_model();
    void run_model();
    const MosquitoManager* get_mos_manager() const {  return &mManager; }
};
