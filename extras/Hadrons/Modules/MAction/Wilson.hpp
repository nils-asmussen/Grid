/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/TWilson.hpp

Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#ifndef Hadrons_Wilson_hpp_
#define Hadrons_Wilson_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            TWilson quark action                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class WilsonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonPar,
                                    std::string, gauge,
                                    double     , mass);
};

template <typename FImpl>
class TWilson: public Module<WilsonPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWilson(const std::string name);
    // destructor
    virtual ~TWilson(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

/******************************************************************************
 *                     TWilson template implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWilson<FImpl>::TWilson(const std::string name)
: Module<WilsonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWilson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWilson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilson<FImpl>::setup(void)
{
    unsigned int size;
    
    size = 3*env().template lattice4dSize<typename FImpl::DoubledGaugeField>();
    env().registerObject(getName(), size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilson<FImpl>::execute()
{
    LOG(Message) << "Setting up TWilson fermion matrix with m= " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
    auto &U      = *env().template getObject<LatticeGaugeField>(par().gauge);
    auto &grid   = *env().getGrid();
    auto &gridRb = *env().getRbGrid();
    FMat *fMatPt = new WilsonFermion<FImpl>(U, grid, gridRb, par().mass);
    env().setObject(getName(), fMatPt);
}

typedef TWilson<FIMPL> Wilson;

END_MODULE_NAMESPACE

MODULE_REGISTER_NS(Wilson, MAction);

END_HADRONS_NAMESPACE

#endif // Hadrons_Wilson_hpp_
