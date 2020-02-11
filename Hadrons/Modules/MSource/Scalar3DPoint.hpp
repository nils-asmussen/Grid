#ifndef Hadrons_MSource_Scalar3DPoint_hpp_
#define Hadrons_MSource_Scalar3DPoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Scalar3DPoint                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class Scalar3DPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Scalar3DPointPar,
                                    std::string, position,
                                    Integer,     tA,
                                    Integer,     tB);
};

template <typename FImpl>
class TScalar3DPoint: public Module<Scalar3DPointPar>
{
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TScalar3DPoint(const std::string name);
    // destructor
    virtual ~TScalar3DPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<int> position_;
};

MODULE_REGISTER_TMP(Scalar3DPoint, TScalar3DPoint<ScalarImplCR>, MSource);

/******************************************************************************
 *                 TScalar3DPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TScalar3DPoint<FImpl>::TScalar3DPoint(const std::string name)
: Module<Scalar3DPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TScalar3DPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TScalar3DPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScalar3DPoint<FImpl>::setup(void)
{
    auto parse_vector = [](const std::string &vec, int dim,
            const std::string &desc)
    {
        std::vector<int> res = strToVec<int>(vec);
        if(res.size() != dim) {
            HADRONS_ERROR(Size, desc + " has "
                    + std::to_string(res.size()) + " instead of "
                    + std::to_string(dim) + " components");
        }
        return res;
    };
    position_ = parse_vector(par().position, env().getNd()-1, "position");

    envCreateLat(PropagatorField, getName());
    //envTmpLat(ComplexField, "component");
    envTmp(LatticeInteger, "latcoord", 1, envGetGrid(ComplexField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScalar3DPoint<FImpl>::execute(void)
{
    auto &rho = envGet(PropagatorField, getName());
    envGetTmp(LatticeInteger, latcoord);

    const int dim=env().getNd()-1;

    rho=1.0;

    //select position
    for(int mu=0; mu<dim; mu++) {
        assert(position_[mu]>=0 && position_[mu]<env().getDim(mu));

        LatticeCoordinate(latcoord, mu);
        rho=where((latcoord==Integer(position_[mu])), rho, 0.*rho);
    }

    //select time slices
    LatticeCoordinate(latcoord, dim);
    rho=where((latcoord>=par().tA && latcoord<=par().tB), rho, 0.*rho);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Scalar3DPoint_hpp_
