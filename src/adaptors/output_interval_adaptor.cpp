#include "output_interval_adaptor.hpp"
#include "../param_manager.hpp"

#include <iostream>

OutputIntervalAdaptor::OutputIntervalAdaptor(unsigned int tStart, unsigned int tStop, unsigned int targetValue)
    : tStart(tStart), tStop(tStop), adaptorName("OutputIntervalAdaptor"), targetValue(targetValue)
{
    baseValue = ParamManager::output_interval;
    ParamManager::recalculate_output_array_size_needed();
}

void OutputIntervalAdaptor::update(unsigned int time)
{
    if (time == tStart)
    {
        ParamManager::output_interval = targetValue;
    }

    if (time == tStop)
    {
        ParamManager::output_interval = baseValue;
    }
}
