/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: GaussDistribution.hpp

Copyright (C) 2015-2019

Author: Nils Asmussen <n.asmussen@soton.ac.uk>

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
#ifndef Hadrons_MUtilities_GaussDistribution_hpp_
#define Hadrons_MUtilities_GaussDistribution_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         GaussDistribution                                  *
 * if fourier=false:                                                          *
 *   result[i,x] = 1/(sqrt(2*pi)*widths[i])^dim*exp(-|x|^2/(2*width[i]^2))    *
 * if fourier=true:                                                           *
 *   result[i,x] = exp(-|x|^2/(2*(1/width[i])^2))                             *
 *
 * where:
 *   x=(x[0],x[1],...,x[dim-1])
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class GaussDistributionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaussDistributionPar,
                                    std::string, widths,
                                    bool,        fourier,
                                    int,         dim);
};

template <typename FImpl>
class TGaussDistribution: public Module<GaussDistributionPar>
{
public:
    // constructor
    TGaussDistribution(const std::string name);
    // destructor
    virtual ~TGaussDistribution(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<Real> widths_;
    void GaussianHelper(LatticeComplex &f, Real fact, int dim);
};

MODULE_REGISTER_TMP(GaussDistribution, TGaussDistribution<FIMPL>, MUtilities);

/******************************************************************************
 *                 TGaussDistribution implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaussDistribution<FImpl>::TGaussDistribution(const std::string name)
: Module<GaussDistributionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaussDistribution<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaussDistribution<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaussDistribution<FImpl>::setup(void)
{
    widths_=strToVec<Real>(par().widths);
    envCreate(std::vector<LatticeComplex>, getName(), 1, widths_.size(),
          envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "component");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaussDistribution<FImpl>::execute(void)
{
    auto &rho = envGet(std::vector<LatticeComplex>, getName());
    assert(widths_.size() == rho.size());
    const bool fourier=par().fourier;
    const int size=widths_.size();
    const int dim=par().dim;
    for(int i=0; i<size; i++) {
        const Real sig=widths_[i];
        auto &r=rho[i];
        const Real var=sig*sig;
        const Real fact=-0.5*( fourier ? var : 1/var );
        GaussianHelper(r, fact, dim);
        if(!fourier) {
            r*=static_cast<Complex>(std::pow(sqrt(2*M_PI)*sig,dim));
        }
    }
}

//exp(fact*|x|^2)
template <typename FImpl>
void TGaussDistribution<FImpl>::GaussianHelper(LatticeComplex &f,
      const Real fact, const int dim) {
   envGetTmp(LatticeComplex, component);
   f=zero;
   for(int mu=0; mu<dim; mu++) {
      LatticeCoordinate(component, mu);
      f+=component*component*fact;
   }
   f=exp(f);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_GaussDistribution_hpp_
