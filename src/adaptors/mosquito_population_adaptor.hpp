#pragma once
#include "adaptor.hpp"

class MosquitoPopulationAdaptor : public Adaptor
{
private:
    unsigned int tStart, tStop;
    std::string adaptorName;
    int targetPopulation;
    MosquitoManager* mManager;
    float fractionalChange = 0.0f; //Keeps track of fractions of a mosquito between updates in order to prevent rounding errors.
public:
    MosquitoPopulationAdaptor(unsigned int tStart, unsigned int tStop, unsigned int targetPopulation, MosquitoManager* mManager);
    void update(unsigned int time);
    unsigned int get_start_t() const { return tStart; }
    unsigned int get_stop_t() const { return tStop; }
    std::string get_adaptor_name() const { return adaptorName; }
    ~MosquitoPopulationAdaptor() {  }
};
