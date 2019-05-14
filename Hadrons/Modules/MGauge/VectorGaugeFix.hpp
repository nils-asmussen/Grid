#ifndef Hadrons_MGauge_VectorGaugeFix_hpp_
#define Hadrons_MGauge_VectorGaugeFix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         VectorGaugeFix                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class VectorGaugeFixPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorGaugeFixPar,
                                    std::string, gaugeXform,
                                    std::string, field);
};

template <typename FImpl>
class TVectorGaugeFix: public Module<VectorGaugeFixPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
    // constructor
    TVectorGaugeFix(const std::string name);
    // destructor
    virtual ~TVectorGaugeFix(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(VectorGaugeFix, TVectorGaugeFix<FIMPL>, MGauge);

/******************************************************************************
 *                 TVectorGaugeFix implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TVectorGaugeFix<FImpl>::TVectorGaugeFix(const std::string name)
: Module<VectorGaugeFixPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TVectorGaugeFix<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field, par().gaugeXform};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TVectorGaugeFix<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TVectorGaugeFix<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TVectorGaugeFix<FImpl>::execute(void)
{
    auto in = envGet(PropagatorField, par().field);
    auto out = envGet(PropagatorField, getName());
    auto g = envGet(GaugeMat, par().gaugeXform);
    out=g*in;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_VectorGaugeFix_hpp_
