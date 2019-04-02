/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: A2ASmearedMesonField.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Nils Asmussen <n.asmussen@soton.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2ASmearedMesonField_hpp_
#define Hadrons_MContraction_A2ASmearedMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ASmearedMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ASmearedMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    bool, allowModifyingLeftRight,
                                    std::string, left,
                                    std::string, right,
                                    std::string, gaugeMatrix,
                                    std::string, smearingLeft,
                                    std::string, smearingRight,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, smearMom,
                                    std::vector<std::string>, mom);
};

class A2ASmearedMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ASmearedMesonFieldMetadata,
                                    std::vector<int>, smearMom,
                                    std::vector<int>, mom,
                                    Gamma::Algebra, gamma);
};

template <typename T, typename FImpl>
class SmearedMesonFieldKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    SmearedMesonFieldKernel(const std::vector<Gamma::Algebra> &gamma,
                     const std::vector<LatticeComplex> &mom,
                     GridBase *grid)
    : gamma_(gamma), mom_(mom), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~SmearedMesonFieldKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left, 
                            const FermionField *right,
                            const unsigned int orthogDim, double &t)
    {
        A2Autils<FImpl>::MesonField(m, left, right, gamma_, mom_, orthogDim, &t);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(2*8.0+6.0+8.0*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
               +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }
private:
    const std::vector<Gamma::Algebra> &gamma_;
    const std::vector<LatticeComplex> &mom_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TA2ASmearedMesonField : public Module<A2ASmearedMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
    typedef A2AMatrixBlockComputation<Complex, 
                                      FermionField, 
                                      A2ASmearedMesonFieldMetadata, 
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef SmearedMesonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2ASmearedMesonField(const std::string name);
    // destructor
    virtual ~TA2ASmearedMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<Gamma::Algebra>       gamma_;
    std::vector<std::vector<int>>     smearMom_;
    std::vector<std::vector<int>>     mom_;
    template<typename Lattice, typename Iterable>
    void multidim_Cshift_inplace(Lattice &lat, const Iterable &shifts);
    void smearing_weight(std::vector<ComplexField> &out,
            const std::vector<ComplexField> &smearingLeft,
            const std::vector<ComplexField> &smearingRight,
            const std::vector<std::vector<int>> &smearMom,
            const std::vector<int> &mom);
};

MODULE_REGISTER(A2ASmearedMesonField, ARG(TA2ASmearedMesonField<FIMPL>), MContraction);

/******************************************************************************
*                  TA2ASmearedMesonField implementation                             
*                  *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ASmearedMesonField<FImpl>::TA2ASmearedMesonField(const std::string name)
: Module<A2ASmearedMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ASmearedMesonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right, par().gaugeMatrix,
        par().smearingLeft, par().smearingRight};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ASmearedMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ASmearedMesonField<FImpl>::setup(void)
{
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }

    auto parse_momenta = [](const std::vector<std::string> &momenta, int dim,
          const std::string &desc)
    {
        std::vector<std::vector<int>> res;
        for (auto &pstr: momenta)
        {
            auto p = strToVec<int>(pstr);
            if (p.size() != dim)
            {
                HADRONS_ERROR(Size, desc + " has " + std::to_string(p.size())
                                    + " components instead of " 
                                    + std::to_string(dim));
            }
            res.push_back(p);
        }
        return res;
    };
    smearMom_ = parse_momenta(par().smearMom, env().getNd()-1,
            "Smearing momentum");
    mom_ = parse_momenta(par().mom, env().getNd()-1,
            "Momentum");

    auto &smearingLeft=envGet(std::vector<LatticeComplex>,
            par().smearingLeft);
    auto &smearingRight=envGet(std::vector<LatticeComplex>,
            par().smearingRight);
    if(smearingLeft.size() != smearingRight.size()) {
        HADRONS_ERROR(Size, "The parameter '" + par().smearingLeft + "' has "
                            + std::to_string(smearingLeft.size())
                            + " elements while parameter '"
                            + par().smearingRight + "' has "
                            + std::to_string(smearingRight.size())
                            + " elements. They must be of equal length.");
    }
    const auto nsmearing=smearingLeft.size();
    const auto smear_size=nsmearing*smearMom_.size();
    envTmp(std::vector<ComplexField>, "smear_weight", 1, smear_size,
            envGetGrid(ComplexField));
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
            env().getNd() - 1, smear_size, gamma_.size(), par().block, 
            par().cacheBlock, this);
    envTmp(FFT, "fft", 1, env().getGrid());

    if(!par().allowModifyingLeftRight)
    {
        auto &left_orig=envGet(std::vector<LatticeFermion>, par().left);
        auto &right_orig=envGet(std::vector<LatticeFermion>, par().right);
        const auto size_l=left_orig.size();
        const auto size_r=right_orig.size();
        envTmp(std::vector<FermionField>, "left", 1, size_l,
            envGetGrid(FermionField));
        envTmp(std::vector<FermionField>, "right", 1, size_r,
            envGetGrid(FermionField));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ASmearedMesonField<FImpl>::execute(void)
{
    auto getA2AField = [&](const std::string &field, const std::string &tmpName)
        -> std::vector<FermionField>&
    {
        if(par().allowModifyingLeftRight)
        {
            return envGet(std::vector<FermionField>, field);
        }
        else
        {
            auto &orig=envGet(std::vector<FermionField>, field);
            //envGetTmp(std::vector<FermionField>, tmpname);
            auto &field=envGet(std::vector<FermionField>,
                    getName() + "_tmp_" + tmpName);
            assert(field.size()==orig.size());
            for(int i=0, size=field.size(); i<size; i++)
            {
                field[i]=orig[i];
            }
            return field;
        }
    };
    auto &left = getA2AField(par().left, "left");
    auto &right = getA2AField(par().right, "right");
    auto &g = envGet(GaugeMat, par().gaugeMatrix);
    auto &smearingLeft
        = envGet(std::vector<LatticeComplex>, par().smearingLeft);
    auto &smearingRight
        = envGet(std::vector<LatticeComplex>, par().smearingRight);
    envGetTmp(FFT, fft);
    envGetTmp(std::vector<ComplexField>, smear_weight);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int ngamma     = gamma_.size();
    int nMomSmear  = smearMom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing smeared all-to-all meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right
        << "'" << std::endl;
    LOG(Message) << "Gauge fixing matrix: '" << par().gaugeMatrix
        << "'" << std::endl;
    LOG(Message) << "Smearing momenta:" << std::endl;
    for (auto &p: smearMom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Relative momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize "
                 << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    //---Fixing gauge---
    LOG(Message) << "Fixing gauge in A2A vectors" << std::endl;
    startTimer("fixing gauge in A2A vectors");
    for(auto &i: left)
    {
        i=g*i;
    }
    for(auto &i: right)
    {
        i=g*i;
    }
    stopTimer("fixing gauge in A2A vectors");

    //---Fourier transform A2A vectors---
    LOG(Message) << "Fourier transforming A2A vectors" << std::endl;
    {
        startTimer("Fourier transform A2A vectors");
        std::vector<int> mask(env().getNd(), 1);
        mask.back()=0; //transform only the spatial dimensions
        for(auto &i: left)
        {
            fft.FFT_dim_mask(i, i, mask, FFT::forward);
        }
        for(auto &i: right)
        {
            fft.FFT_dim_mask(i, i, mask, FFT::forward);
        }
        stopTimer("Fourier transform A2A vectors");
    }

    for(const auto &m: mom_)
    {

        auto ionameFn = [this,&m,nMomSmear](const unsigned int ms,
                const unsigned int g)
        {
            std::stringstream ss;

            ss << gamma_[g] << "_smear";
            for(auto pmu: smearMom_[ms%nMomSmear])
            {
                ss << '_' << pmu;
            }
            ss << "_mom";
            for(auto pmu: m)
            {
                ss << '_' << pmu;
            }
            ss << "_smearing_" << ms/nMomSmear;

            return ss.str();
        };

        auto filenameFn = [this, &ionameFn](const unsigned int ms, const unsigned int g)
        {
            return par().output + "." + std::to_string(vm().getTrajectory()) 
                   + "/" + ionameFn(ms, g) + ".h5";
        };

        auto metadataFn = [this,&m,nMomSmear](const unsigned int ms,
                const unsigned int g)
        {
            A2ASmearedMesonFieldMetadata md;

            for (auto pmu: smearMom_[ms%nMomSmear])
            {
                md.smearMom.push_back(pmu);
            }
            for (auto pmu: m)
            {
                md.mom.push_back(pmu);
            }
            md.gamma = gamma_[g];
            
            return md;
        };

        startTimer("Cshift momentum A2A vector");
        for(auto &field: right)
        {
            multidim_Cshift_inplace(field, m);
        }
        stopTimer("Cshift momentum A2A vector");

        startTimer("compute smearing weight");
        smearing_weight(smear_weight, smearingLeft, smearingRight,
                smearMom_, m);
        stopTimer("compute smearing weight");
        Kernel kernel(gamma_, smear_weight, envGetGrid(FermionField));

        envGetTmp(Computation, computation);
        computation.execute(left, right, kernel,
                ionameFn, filenameFn, metadataFn);

        //restore A2A vectors
        {
            startTimer("Cshift momentum A2A vector");
            auto mNeg(m);
            for(auto &i: mNeg)
            {
                i*=-1;
            }
            for(auto &field: right)
            {
                multidim_Cshift_inplace(field, mNeg);
            }
            stopTimer("Cshift momentum A2A vector");
        }
    }
}

template<typename FImpl>
template<typename Lattice, typename Iterable>
void TA2ASmearedMesonField<FImpl>::multidim_Cshift_inplace(Lattice &lat,
        const Iterable &shifts)
{
    int dim=0;
    for(auto shift: shifts)
    {
        if(shift!=0)
        {
            lat=Cshift(lat, dim, shift);
        }
        dim++;
    }
}

//compute the smearing weight
//out[n*Nm+m](q)=smearingLeft[n](q)*smearingRight[n](q+p-smearMom[m])
//where
//Nm=smearMom.size()
//q: coordinate of the (momentum space) field
//p=mom
template<typename FImpl>
void TA2ASmearedMesonField<FImpl>::smearing_weight(
        std::vector<ComplexField> &out,
        const std::vector<ComplexField> &smearingLeft,
        const std::vector<ComplexField> &smearingRight,
        const std::vector<std::vector<int>> &smearMom,
        const std::vector<int> &mom)
{
    assert(smearingLeft.size() == smearingRight.size());
    assert(out.size() == smearingLeft.size()*smearMom.size());
    const int numSpatialDims=env().getNd()-1;

    const int spatialVol=[this]{
        int ret=1;
        std::vector<int> lattSizeSpatial=env().getGrid()->FullDimensions();
        lattSizeSpatial.pop_back();
        for(auto i : lattSizeSpatial)
        {
            ret*=i;
        }
        return ret;
    }();
    const Complex invSpVol(1./spatialVol);

    std::vector<int> totalShift(numSpatialDims);
    auto out_it  =out.begin();
    auto left_it =smearingLeft.begin();
    auto right_it=smearingRight.begin();
    auto left_end=smearingLeft.end();
    for( ; left_it!=left_end ; ++left_it, ++right_it)
    {
        for(const auto &msmear: smearMom)
        {
            assert(msmear.size() == numSpatialDims);
            assert(mom.size() == numSpatialDims);
            for(int i=0; i<numSpatialDims; i++)
            {
                totalShift[i]=mom[i]-msmear[i];
            }
            *out_it=*right_it;
            multidim_Cshift_inplace(*out_it, totalShift);
            *out_it=*out_it * *left_it * invSpVol;
            out_it++;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ASmearedMesonField_hpp_
