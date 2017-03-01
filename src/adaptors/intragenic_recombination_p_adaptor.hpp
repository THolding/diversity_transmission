#pragma once
#include "adaptor.hpp"

class IntragenicRecombinationPAdaptor : public Adaptor
{
private:
    unsigned int tStart, tStop;
    std::string adaptorName;
    float targetValue;
public:
    IntragenicRecombinationPAdaptor(unsigned int tStart, unsigned int tStop, float targetValue);
    void update(unsigned int time);
    unsigned int get_start_t() const { return tStart; }
    unsigned int get_stop_t() const { return tStop; }
    std::string get_adaptor_name() const { return adaptorName; }
    ~IntragenicRecombinationPAdaptor() {  }
};

