/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::combustionModels::baseFSD

Description
    Partially stirred reactor turbulent combustion model.

    This model calculates a finite rate, based on both turbulence and chemistry
    time scales.  Depending on mesh resolution, the Cmix parameter can be used
    to scale the turbulence mixing time scale.

SourceFiles
    baseFSD.C

\*---------------------------------------------------------------------------*/

#ifndef baseFSD_H
#define baseFSD_H

#include "../laminar/laminar.H"
#include "uniformDimensionedFields.H"
#include "CanteraMixture.H"
// #include "fvOptions.H"
#include "LESfilter.H"
//#include "basicThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class baseFSD Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class baseFSD
:
    public laminar<ReactionThermo>
{
    // Protected Data
 protected:

    Switch buffer_;

    //- incompressible pressure reference
    scalar incompPref_;

    //- ignition switch
    Switch ignition_;

    //- combustion switch
    Switch combustion_;

    //- enthalpy solve witch
    Switch solveEnthalpy_;
    
    //- T-lookup switch
    Switch flameletT_;

    word tablePath_;

    //- psi=\rho/p
    volScalarField& psi_;

    //- fluid density
    volScalarField& rho_;
    volScalarField& rho_inRhoThermo_;

    //- fluid pressure
    const volScalarField& p_;

    //- fluid temperature
    volScalarField& T_;

    PtrList<volScalarField>& Y_;

    //- fluid velocity
    const volVectorField& U_;    

    //- air mass/mole(28.96kg/mol)
    volScalarField Wt_;
    
    //- heat capacity(m^2*s^{-2}*K^{-1}) 
    volScalarField Cp_;

    //- mixture fraction Z
    volScalarField Z_;

    volScalarField Y_fuel_;

    //- mixture fraction variance Zvar
    volScalarField Zvar_;    

    //- progress variable c variance
    volScalarField cvar_;

    //- co-variance Z, c
    volScalarField Zcvar_;
    
    //- enthalpy 
    volScalarField He_;

    //- formation enthalpy(m^2*s^{-2})
    volScalarField Hf_;

    volScalarField mu_tab_;

    //- progress variable c
    volScalarField c_;

    //- flame surface density
    volScalarField fsd_;

    //- flame surface density
    volScalarField SdA_;

    //- strain factor characterizing the effects of flame straining on the flame speed
    volScalarField I_s_FSD_;

    //- laminar flame speed
    volScalarField SL0_FSD_;

    //- termal thickness of the unstrained laminar flame
    volScalarField thermoThickness_L0_FSD_;

    //- Kalovitz factor
    volScalarField Ka_;

    //- temp field to save Ka_, s^{-1}
    volScalarField Ka_temp_;

    volScalarField Ka_temp1_;
    volScalarField Ka_temp2_;
    volScalarField Ka_temp3_;

    //- chemistry source reaction rate(kg*m^{-3}*s^{-1})
    volScalarField omega_c_;

    volScalarField omega_c_chiZ_;

    // volScalarField omega_Y_fuel_;

    //- mixture fraction dissipation rate (s^{-1})
    volScalarField chi_Z_; 

    //- progress variable dissipation rate (s^{-1})
    volScalarField chi_c_;

    volScalarField chi_Zc_; 

    //- filtered SDR of mixture fraction
    volScalarField chi_Zfltd_; 

    //- betac for LES
    volScalarField betac_;

    volScalarField He_s_;
    volScalarField c_s_;
    volScalarField Zvar_s_;

    //- unburnt density
    volScalarField rho_u_;

    //- filter size
    volScalarField delta_;

    //- efficiency function
    volScalarField gamma_k_;

    //- the 1st source term in FSD equation, representing the flame curvature effects
    volScalarField source_fsd_FCE1_;

    //- the 2nd source term in FSD equation, representing the flame curvature effects
    volScalarField source_fsd_FCE2_;

    //- source term in FSD equation, representing the flame propagation
    volScalarField source_fsd_flaProp_;

    //- source term in FSD equation, representing the strain-rate effects
    volScalarField source_fsd_SRE1_;
    volScalarField source_fsd_SRE2_;

    //- flux term for FSD equation
    surfaceScalarField phiFSD_;

    //- unit tensor
    volTensorField IField_;

    //- unit tensor with m^{-1}
    volScalarField IField_m_1_;

    //- unit tensor with m
    volScalarField IField_m1_;

    //- unit tensor with m/s
    volScalarField IField_m1_s_1_;

    //- modeled surface-averaged normal vector
    volVectorField n_FSD_;

    volScalarField cOmega_c_;
    volScalarField ZOmega_c_;

    //- internal fields
    scalarField& WtCells_; 
    scalarField& CpCells_;
    scalarField& ZCells_;
    scalarField& ZvarCells_;
    scalarField& cvarCells_;
    scalarField& ZcvarCells_;
    scalarField& HCells_;
    scalarField& HfCells_;
    scalarField& cCells_;
    scalarField& fsdCells_;
    scalarField& KaCells_;
    scalarField& omega_cCells_;
    scalarField& cOmega_cCells_;
    scalarField& ZOmega_cCells_;
    scalarField& omega_c_chiZCells_;
    // scalarField& omega_Y_fuelCells_;
    scalarField& chi_ZCells_;
    scalarField& chi_cCells_;
    scalarField& chi_ZcCells_;
    scalarField& chi_ZfltdCells_;
    scalarField& betacCells_; 
    scalarField& rho_uCells_;
    scalarField& SL0_FSDCells_;
    scalarField& thermoThickness_L0_FSDCells_;
    scalarField& I_s_FSDCells_;
    vectorField& n_FSDCells_;
    scalarField& gamma_kCells_;
    scalarField& source_fsd_FCE2Cells_;
    scalarField& SdACells_;
    scalarField& Y_fuelCells_;

    //- filter size
    scalarField deltaCells_; 

    //- max and min Z
    scalar ZMax_;
    scalar ZMin_;

    //- max and min c
    scalar cMax_;
    scalar cMin_;

    //- max and min Z variance
    scalar ZvarMax_;
    scalar ZvarMin_;

    //- max and min c variance
    scalar cvarMax_;
    scalar cvarMin_;

    //- max and min Zc variance 
    scalar ZcvarMax_;
    scalar ZcvarMin_;   

    //- max and min T
    scalar TMax_;
    scalar TMin_; 

    //- error limitation for FSD
    scalar small_FSD_;

    // pressure derivate
    const volScalarField& dpdt_;

    //- fluid flux
    const surfaceScalarField& phi_;

    scalarField& TCells_;    

    //- ignition begin time
    scalar ignBeginTime_;

    //- ignition duration time
    scalar ignDurationTime_;    
 
    //- reactFlowTime_=currentTime - ignBeginTime_
    scalar reactFlowTime_;

    //- igntion origin point
    scalar x0_;
    scalar y0_;
    scalar z0_;

    //- ignition radius
    scalar R0_;

    //- turbulent schmidt number
    scalar Sct_;

    //- laminar schmidt number
    scalar Sc_;

    //- FSD Schmidt number
    scalar Sc_fsd_;

    //- FSD model constants
    scalar beta_FSD_;

    scalar phi_fac_FSD_;

    scalar bufferTime_;
 
    //- relaxation switch
    Switch relaxation_;

    //- if SdA consider the flame curcature term
    Switch curv_SdA_;  

    //- Dpdt switch 
    Switch DpDt_;    

    //- index for omega speices
    labelList omega_Yi_index_;

    wordList omega_YiNames_base_;

    // //- reaction rate for fuel species
    PtrList<volScalarField> omega_Yis_;

    //- LES
    //- if in LES mode
    Switch isLES_;

    //SGS
    volScalarField magUPrime_;
    scalarField& magUPrimeCells_; 

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

 public:
    //- Runtime type information
    TypeName("baseFSD");

    // Constructors

        //- Construct from components
        baseFSD
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~baseFSD();


    // Member Functions

        //- transport equation 
        virtual void transport();

        //- initialize flame kernel 
        virtual void initialiseFalmeKernel();    

        //- return mag SGS
        virtual void magUPrime();    

        //- reaturn buffer_
        virtual bool buffer();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "baseFSD.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
