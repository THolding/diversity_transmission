#pragma once
#include "adaptors/adaptor.hpp"
#include "global_typedefs.hpp"
#include <array>
#include <list>
#include <string>
#include <unordered_map>

class Adaptor;

//Based on http://stackoverflow.com/questions/1008019/c-singleton-design-pattern
class ParamManager
{
private:
    std::string runName = "default";
    std::string filePath = "";
    std::unordered_map<std::string, float> paramsFloat;
    std::unordered_map<std::string, unsigned int> paramsInt;
    std::unordered_map<std::string, bool> paramsBool;
    std::array<float, 2> recombination_cumu_p;

    BITE_FREQUENCY_TABLE cumulativeBiteFrequencyDistribution;

    std::list<Adaptor*> adaptors;

    ParamManager();

public:
    static ParamManager& instance()
    {
        static ParamManager paramManager; //Implicit if flag ensures initialisation occurs only once!
        return paramManager;
    }

    void set_run_name(std::string _runName) { runName = _runName; }
    void set_file_path(std::string _filePath) { filePath = _filePath; }
    std::string run_name() { return runName; }
    std::string file_path() { return filePath; }
    float get_float(const std::string paramName) const;
    unsigned int get_int(const std::string paramName) const;
    bool get_bool(const std::string paramName) const;
    void set_param(const std::string name, const std::string value);
    void set_float(const std::string paramName, const float param);
    void set_int(const std::string paramName, const int param);
    void set_bool(const std::string paramName, const bool param);

    void add_adaptor(Adaptor* adaptor);
    bool is_compatable_adaptor(const Adaptor& adaptor) const;

    void update_adaptors(const unsigned int t);

    bool recalculate_derived_parameters();
    bool recalculate_recombination_distributions();
    void recalculate_cumulative_bite_frequency_distribution();
    const BITE_FREQUENCY_TABLE& get_cumulative_bite_frequency_distribution() const { return cumulativeBiteFrequencyDistribution;};
    void recalculate_output_array_size_needed();

    ParamManager(ParamManager const&) = delete; //disable copy construction
    void operator=(ParamManager const&) = delete; //disable copy assignment

    ~ParamManager();
};
