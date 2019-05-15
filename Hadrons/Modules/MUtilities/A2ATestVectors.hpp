#ifndef Hadrons_MUtilities_A2ATestVectors_hpp_
#define Hadrons_MUtilities_A2ATestVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ATestVectors                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ATestVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ATestVectorsPar,
                                    unsigned int, numVec,
                                    std::string, contents,
                                    unsigned int, dim);
};

template <typename FImpl>
class TA2ATestVectors: public Module<A2ATestVectorsPar>
{
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ATestVectors(const std::string name);
    // destructor
    virtual ~TA2ATestVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<Complex> contents_;
};

MODULE_REGISTER_TMP(A2ATestVectors, TA2ATestVectors<FIMPL>, MUtilities);

/******************************************************************************
 *                 TA2ATestVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ATestVectors<FImpl>::TA2ATestVectors(const std::string name)
: Module<A2ATestVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ATestVectors<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ATestVectors<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()+"_v", getName()+"_w"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ATestVectors<FImpl>::setup(void)
{
    envCreate(std::vector<FermionField>, getName()+"_v", 1, par().numVec, envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName()+"_w", 1, par().numVec, envGetGrid(FermionField));
    contents_=strToVec<Complex>(par().contents);

    std::vector<int> lattSize=env().getGrid()->FullDimensions();
    assert(lattSize.size()>par().dim);
    assert(lattSize[par().dim]>=contents_.size());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ATestVectors<FImpl>::execute(void)
{
    auto &v=envGet(std::vector<FermionField>, getName()+"_v");
    auto &w=envGet(std::vector<FermionField>, getName()+"_w");
    int dim=par().dim;
    std::vector<int> dims {0,1,2,3};
    dims.erase(dims.begin()+dim);

    Complex val (1.0, 0);
    ColourVector cv=zero;
    SpinColourVector scv=zero;
    pokeColour(cv, val, 0);
    pokeSpin(scv, cv, 0);

    auto set_vector=[&](std::vector<FermionField> &vector) {
       std::vector<int> lattSize=env().getGrid()->FullDimensions();
       for(auto &i : vector) {
          i=zero;
          std::vector<int> pos(4,0);
          for(pos[dims[0]]=0; pos[dims[0]]<lattSize.at(dims[0]); pos[dims[0]]++)
          for(pos[dims[1]]=0; pos[dims[1]]<lattSize.at(dims[1]); pos[dims[1]]++)
          for(pos[dims[2]]=0; pos[dims[2]]<lattSize.at(dims[2]); pos[dims[2]]++) {
             for(auto &c : contents_) {
                pokeSite(c*scv, i, pos);
                pos.at(dim)++;
             }
          }
       }
    };
    set_vector(v);
    set_vector(w);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ATestVectors_hpp_
