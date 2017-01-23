#pragma once
#include <string>

class MosquitoManager;

//Virtual base class for all adaptors.
class Adaptor
{
protected:
    //std::string paramName;
    //unsigned int startT;
    //unsigned int stopT;
public:
    virtual void update(unsigned int t) = 0;
    virtual unsigned int get_start_t() const = 0;
    virtual unsigned int get_stop_t() const = 0;
    virtual std::string get_adaptor_name() const = 0;
    virtual ~Adaptor() {  };
};
