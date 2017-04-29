
#include "intragenic_recombination_p_adaptor.hpp"
#include "../param_manager.hpp"

#include <iostream>

IntragenicRecombinationPAdaptor::IntragenicRecombinationPAdaptor(unsigned int tStart, unsigned int tStop, float targetValue)
    : tStart(tStart), tStop(tStop), adaptorName("IntragenicRecombinationPAdaptor"), targetValue(targetValue)
{
    ParamManager::dyn_intragenic_recombination_p = true;
    if (targetValue < 0.0f || targetValue > 1.0f)
        throw std::runtime_error("IntragenicRecombinationPAdaptor cannot have intragenic recombination probability of less than 0 or greater than 1.");
}

void IntragenicRecombinationPAdaptor::update(unsigned int time)
{
    //std::cout << "***bite_rate (current): " << ParamManager::instance().get_float("bite_rate") << "\n";
    if (time >= tStart && time < tStop)
    {
        float changeLeft = targetValue - ParamManager::intragenic_recombination_p;
        unsigned int timeLeft = tStop - time;
        float changeToDo = changeLeft / (float) timeLeft;

        /*std::cout << "\nCurrent value: " << ParamManager::instance().get_float("intragenic_recombination_p") << "\n";
        std::cout << "Target value: " << targetValue << "\n";
        std::cout << "change left: " << changeLeft << "\n";
        std::cout << "change to do: " << changeToDo << "\n";*/

        //std::cout << "Previous val: " << ParamManager::instance().get_float("intragenic_recombination_p");
        ParamManager::intragenic_recombination_p += changeToDo;
        //std::cout << "\tNew val: " << ParamManager::instance().get_float("intragenic_recombination_p") << "\n";
        ParamManager::recalculate_recombination_distributions();
    }
}
