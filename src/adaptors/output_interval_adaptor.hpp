#pragma once
#include "adaptor.hpp"

//Unlike other adaptors this is a step function that changes the output_interval between (tStart, tStop] then returns it to the base value (whatever the model was initialised with).
class OutputIntervalAdaptor : public Adaptor
{
private:
    unsigned int tStart, tStop;
    std::string adaptorName;
    unsigned int targetValue;
    unsigned int baseValue; //Value to return to is always the base value used at the start of model initialisation.
public:
    OutputIntervalAdaptor(unsigned int tStart, unsigned int tStop, unsigned int targetValue);
    void update(unsigned int time);
    unsigned int get_start_t() const { return tStart; }
    unsigned int get_stop_t() const { return tStop; }
    unsigned int get_target_output_interval() const { return targetValue; }
    std::string get_adaptor_name() const { return adaptorName; }
    ~OutputIntervalAdaptor() {  }
};
