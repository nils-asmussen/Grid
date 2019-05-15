#ifndef Hadrons_MUtilities_JacobiSmear_hpp_
#define Hadrons_MUtilities_JacobiSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         JacobiSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class JacobiSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(JacobiSmearPar,
                                    std::string, U,
                                    double, width,
                                    int, iterations,
                                    int, orthog,
                                    std::string, source);
};

template <typename FImpl>
class TJacobiSmear: public Module<JacobiSmearPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TJacobiSmear(const std::string name);
    // destructor
    virtual ~TJacobiSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(JacobiSmear, TJacobiSmear<FIMPL>, MUtilities);

/******************************************************************************
 *                 TJacobiSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TJacobiSmear<FImpl>::TJacobiSmear(const std::string name)
: Module<JacobiSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().U};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TJacobiSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TJacobiSmear<FImpl>::execute(void)
{
    auto &out = envGet(PropagatorField, getName());
    auto &src = envGet(PropagatorField, par().source);
    auto &U = envGet(std::vector<LatticeColourMatrix>, par().U);
    CovariantSmearing<FImpl> covsmear;
    out=src;
    startTimer("Jacobi iteration");
    covsmear.GaussianSmear(U, out, par().width, par().iterations, par().orthog);
    stopTimer("Jacobi iteration");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_JacobiSmear_hpp_
