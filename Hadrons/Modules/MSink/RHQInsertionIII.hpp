/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSink/Smear.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MSink_RHQInsertionIII_hpp_
#define Hadrons_MSink_RHQInsertionIII_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 RHQInsertionIII                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class RHQInsertionIIIPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIIIPar,
                                    std::string, q,
                                    std::string, gauge);
};

template <typename FImpl, typename GImpl>
class TRHQInsertionIII: public Module<RHQInsertionIIIPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,)
    SINK_TYPE_ALIASES();
public:
    // constructor
    TRHQInsertionIII(const std::string name);
    // destructor
    virtual ~TRHQInsertionIII(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionIII, ARG(TRHQInsertionIII<FIMPL, GIMPL>), MSink);

/******************************************************************************
 *                          TSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionIII<FImpl, GImpl>::TRHQInsertionIII(const std::string name)
: Module<RHQInsertionIIIPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIII<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIII<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIII<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());//, 1, env().getDim(Tp));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIII<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Applying Improvement term III to'" << par().q
                 << std::endl;

    auto &field    = envGet(PropagatorField, par().q);
    const auto &gaugefield  = envGet(GaugeField, par().gauge);
    const auto gauge_x = peekLorentz(gaugefield, 0);
    const auto gauge_y = peekLorentz(gaugefield, 1);
    const auto gauge_z = peekLorentz(gaugefield, 2);

    Gamma gx(Gamma::Algebra::GammaX);
    Gamma gy(Gamma::Algebra::GammaY);
    Gamma gz(Gamma::Algebra::GammaZ);
    PropagatorField  insertion = 
        gx * (GImpl::CovShiftForward( gauge_x,0,field) - GImpl::CovShiftBackward( gauge_x,0,field))
     + gy * (GImpl::CovShiftForward( gauge_y,1,field) - GImpl::CovShiftBackward( gauge_y,1,field))
     + gz * (GImpl::CovShiftForward( gauge_z,2,field) - GImpl::CovShiftBackward( gauge_z,2,field));
     
    auto &out  = envGet(PropagatorField, getName());
    out = insertion;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_RHQInsertionIII_hpp_
