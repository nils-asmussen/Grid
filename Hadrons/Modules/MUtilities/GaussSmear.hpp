#ifndef Hadrons_MUtilities_GaussSmear_hpp_
#define Hadrons_MUtilities_GaussSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         GaussSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class GaussSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaussSmearPar,
                                    std::string, distribution,
                                    std::string, source);
};

template <typename FImpl>
class TGaussSmear: public Module<GaussSmearPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGaussSmear(const std::string name);
    // destructor
    virtual ~TGaussSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaussSmear, TGaussSmear<FIMPL>, MUtilities);

/******************************************************************************
 *                 TGaussSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaussSmear<FImpl>::TGaussSmear(const std::string name)
: Module<GaussSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaussSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().distribution};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaussSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaussSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmp(FFT, "fft", 1, env().getGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaussSmear<FImpl>::execute(void)
{
    auto &dist = envGet(std::vector<LatticeComplex>, par().distribution);
    auto &out = envGet(PropagatorField, getName());
    auto &src = envGet(PropagatorField, par().source);
    envGetTmp(FFT, fft);
    startTimer("Fourier transform");
    std::vector<int> mask(env().getNd(), 1);
    mask.back()=0; //transform only the spatial dimensions
    fft.FFT_dim_mask(out, src, mask, FFT::forward);
    stopTimer("Fourier transform");
    startTimer("applying distribution");
    out=dist.at(0)*out;
    stopTimer("applying distribution");
    startTimer("Fourier transform");
    fft.FFT_dim_mask(out, out, mask, FFT::backward);
    stopTimer("Fourier transform");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_GaussSmear_hpp_
