#include "bite_rate_adaptor.hpp"
#include "../param_manager.hpp"

#include <iostream>

BiteRateAdaptor::BiteRateAdaptor(unsigned int tStart, unsigned int tStop, float targetValue)
    : tStart(tStart), tStop(tStop), adaptorName("BiteRateAdaptor"), targetValue(targetValue)
{
    ParamManager::instance().set_bool("dyn_bite_rate", true);
    if (targetValue < 0.0f)
        throw std::runtime_error("BiteRateAdaptor cannot have target bite rate of less than 1.");
}

void BiteRateAdaptor::update(unsigned int time)
{
    //std::cout << "***bite_rate (current): " << ParamManager::instance().get_float("bite_rate") << "\n";
    if (time >= tStart && time < tStop)
    {
        float changeLeft = targetValue - ParamManager::instance().get_float("bite_rate");
        unsigned int timeLeft = tStop - time;
        float changeToDo = changeLeft / (float) timeLeft;

        ParamManager::instance().set_float("bite_rate", ParamManager::instance().get_float("bite_rate")+changeToDo);
        ParamManager::instance().recalculate_cumulative_bite_frequency_distribution();
    }
}
