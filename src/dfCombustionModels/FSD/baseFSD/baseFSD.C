/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "baseFSD.H"
#include "basicSprayCloud.H" 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFSD<ReactionThermo>::baseFSD
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    buffer_(this->coeffs().lookupOrDefault("buffer", false)),
    incompPref_(this->coeffs().lookupOrDefault("incompPref", -10.0)),
    ignition_(this->coeffs().lookupOrDefault("ignition", false)),
    combustion_(this->coeffs().lookupOrDefault("combustion", false)),
    solveEnthalpy_(this->coeffs().lookupOrDefault("solveEnthalpy", false)),
    flameletT_(this->coeffs().lookupOrDefault("flameletT", false)),
    tablePath_(this->coeffs().lookup("tablePath")),
    psi_(const_cast<volScalarField&>(dynamic_cast<rhoThermo&>(thermo).psi())),
    rho_(const_cast<volScalarField&>(this->mesh().objectRegistry::lookupObject<volScalarField>("rho"))),
    rho_inRhoThermo_(dynamic_cast<rhoThermo&>(thermo).rho()),
    p_(this->thermo().p()),
    T_(this->thermo().T()),
    Y_(this->chemistryPtr_->Y()),
    U_(this->mesh().objectRegistry::lookupObject<volVectorField>("U")),
    Wt_ 
    (
        IOobject
        (
            "Wt",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Wt",dimensionSet(1,0,0,0,-1,0,0),28.96)  
    ),
    Cp_ 
    (
        IOobject
        (
            "Cp",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Cp",dimensionSet(0,2,-2,-1,0,0,0),1010.1)
    ),
    Z_
    (
        IOobject
        (
            "Z",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    Y_fuel_
    (
        IOobject
        (
            "Y_fuel",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    Zvar_
    (
        IOobject
        (
            "Zvar",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    cvar_
    (
        IOobject
        (
            "cvar",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),  
    Zcvar_
    (
        IOobject
        (
            "Zcvar",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    He_
    (
        IOobject
        (
            "Ha",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ha",dimensionSet(0,2,-2,0,0,0,0),0.0)
    ),
    Hf_
    (
        IOobject
        (
            "Hf",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Hf",dimensionSet(0,2,-2,0,0,0,0),1907.0)
    ),
    mu_tab_
    (
        IOobject
        (
            "mu_tab",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("mu_tab",dimensionSet(1,-1,-1,0,0,0,0),0.0)
    ),
    c_
    (
        IOobject
        (
            "c_PV",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("c_PV",dimensionSet(0,0,0,0,0,0,0),0.0) 
    ), 
    fsd_
    (
        IOobject
        (
            "fsd",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("fsd",dimensionSet(0,-1,0,0,0,0,0),0.0)
    ), 
    SdA_
    (
        IOobject
        (
            "SdA",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("SdA",dimensionSet(0,1,-1,0,0,0,0),0.0)
    ), 
    I_s_FSD_
    (
        IOobject
        (
            "I_s_FSD",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("I_s_FSD",dimless,0.0)
    ), 
    SL0_FSD_
    (
        IOobject
        (
            "SL0_FSD",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("SL0_FSD",dimensionSet(0,1,-1,0,0,0,0),0.0)
    ), 
    thermoThickness_L0_FSD_
    (
        IOobject
        (
            "thermoThickness_L0_FSD",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("thermoThickness_L0_FSD",dimensionSet(0,1,0,0,0,0,0),0.0)
    ), 
    Ka_
    (
        IOobject
        (
            "Ka",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ka",dimensionSet(0,0,0,0,0,0,0),0.0)
    ), 
    Ka_temp_
    (
        IOobject
        (
            "Ka_temp",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ka_temp",dimensionSet(0,0,-1,0,0,0,0),0.0)
    ), 
    Ka_temp1_
    (
        IOobject
        (
            "Ka_temp1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ka_temp1",dimensionSet(0,0,-1,0,0,0,0),0.0)
    ), 
    Ka_temp2_
    (
        IOobject
        (
            "Ka_temp2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ka_temp2",dimensionSet(0,0,-1,0,0,0,0),0.0)
    ), 
    Ka_temp3_
    (
        IOobject
        (
            "Ka_temp3",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ka_temp3",dimensionSet(0,0,-1,0,0,0,0),0.0)
    ), 
    omega_c_ 
    (
        IOobject
        (
            "omega_c",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("omega_c",dimensionSet(1,-3,-1,0,0,0,0),0.0)
    ),
    omega_c_chiZ_
    (
        IOobject
        (
            "omega_c_chiZ",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("omega_c_chiZ",dimensionSet(1,-3,-1,0,0,0,0),0.0)
    ),
    chi_Z_
    (
        IOobject
        (
            "chi_Z",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("chi_Z",dimensionSet(0,0,-1,0,0,0,0),0.0) 
    ),  
    chi_c_
    (
        IOobject
        (
            "chi_c",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("chi_c",dimensionSet(0,0,-1,0,0,0,0),0.0)
    ),   
    chi_Zc_
    (
        IOobject
        (
            "chi_Zc",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("chi_Zc",dimensionSet(0,0,-1,0,0,0,0),0.0) 
    ), 
    chi_Zfltd_
    (
        IOobject
        (
            "chi_Zfltd",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("chi_Zfltd",dimensionSet(0,0,-1,0,0,0,0),0)
    ),
    betac_  
    (
        IOobject
        (
            "betac",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("betac",dimensionSet(0,0,0,0,0,0,0),4.0) 
    ), 
    He_s_
    (
        IOobject
        (
            "Ha_s",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Ha_s",dimensionSet(0,2,-2,0,0,0,0),0.0)
    ),
    c_s_
    (
        IOobject
        (
            "c_s",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("c_s",dimensionSet(0,0,0,0,0,0,0),0.0) 
    ), 
    Zvar_s_
    (
        IOobject
        (
            "Zvar_s",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("Zvar_s",dimensionSet(0,0,0,0,0,0,0),0.0) 
    ),
    rho_u_
    (
        IOobject
        (
            "rho_u",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("rho_u",dimensionSet(1,-3,0,0,0,0,0),0.0) 
    ),
    delta_
    (
        IOobject
        (
            "delta",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("delta",dimensionSet(0,1,0,0,0,0,0),0.0) 
    ),
    gamma_k_
    (
        IOobject
        (
            "gamma_k",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("gamma_k",dimless,0.0) 
    ),
    source_fsd_FCE1_
    (
        IOobject
        (
            "source_fsd_FCE1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("source_fsd_FCE1",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
    ),
    source_fsd_FCE2_
    (
        IOobject
        (
            "source_fsd_FCE2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("source_fsd_FCE2",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
    ),
    source_fsd_flaProp_
    (
        IOobject
        (
            "source_fsd_flaProp",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("source_fsd_flaProp",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
    ),
    source_fsd_SRE1_
    (
        IOobject
        (
            "source_fsd_SRE1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("source_fsd_SRE1",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
    ),
    source_fsd_SRE2_
    (
        IOobject
        (
            "source_fsd_SRE2",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("source_fsd_SRE2",dimensionSet(0,-1,-1,0,0,0,0),0.0) 
    ),
    phiFSD_
    (
        IOobject
        (
            "phiFSD",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & this->mesh().Sf()
    ),
    IField_
    (
        IOobject
        (
            "IField",         // 名称
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedTensor("IField", dimless, tensor::I)
    ),
    IField_m_1_
    (
        IOobject
        (
            "IField_m_1",         // 名称
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("IField_m_1", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0)
    ),
    IField_m1_
    (
        IOobject
        (
            "IField_m1",         // 名称
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("IField_m1", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1.0)
    ),
    IField_m1_s_1_
    (
        IOobject
        (
            "IField_m1_s_1",         // 名称
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("IField_m1_s_1", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1.0)
    ),
    n_FSD_
    (
        IOobject
        (
            "n_FSD",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedVector(dimensionSet(0, 0, 0, 0, 0, 0, 0), vector::zero)
    ),
    cOmega_c_(omega_c_),
    ZOmega_c_(omega_c_),
    WtCells_ (Wt_.primitiveFieldRef()),
    CpCells_ (Cp_.primitiveFieldRef()),
    ZCells_(Z_.primitiveFieldRef()),
    ZvarCells_(Zvar_.primitiveFieldRef()), 
    cvarCells_(cvar_.primitiveFieldRef()), 
    ZcvarCells_(Zcvar_.primitiveFieldRef()), 
    HCells_(He_.primitiveFieldRef()),
    HfCells_(Hf_.primitiveFieldRef()),
    cCells_(c_.primitiveFieldRef()),
    fsdCells_(fsd_.primitiveFieldRef()),
    KaCells_(Ka_.primitiveFieldRef()),
    omega_cCells_(omega_c_.primitiveFieldRef()), 
    cOmega_cCells_(cOmega_c_.primitiveFieldRef()),
    ZOmega_cCells_(ZOmega_c_.primitiveFieldRef()),
    omega_c_chiZCells_(omega_c_chiZ_.primitiveFieldRef()), 
    // omega_Y_fuelCells_(omega_Y_fuel_.primitiveFieldRef()), 
    chi_ZCells_(chi_Z_.primitiveFieldRef()),
    chi_cCells_(chi_c_.primitiveFieldRef()),
    chi_ZcCells_(chi_Zc_.primitiveFieldRef()),
    chi_ZfltdCells_(chi_Zfltd_.primitiveFieldRef()),  
    betacCells_(betac_.primitiveFieldRef()),   
    rho_uCells_(rho_u_.primitiveFieldRef()),
    SL0_FSDCells_(SL0_FSD_.primitiveFieldRef()),
    thermoThickness_L0_FSDCells_(thermoThickness_L0_FSD_.primitiveFieldRef()),
    I_s_FSDCells_(I_s_FSD_.primitiveFieldRef()),
    n_FSDCells_(n_FSD_.primitiveFieldRef()),
    gamma_kCells_(gamma_k_.primitiveFieldRef()),
    source_fsd_FCE2Cells_(source_fsd_FCE2_.primitiveFieldRef()),
    SdACells_(SdA_.primitiveFieldRef()),
    Y_fuelCells_(Y_fuel_.primitiveFieldRef()),
    ZMax_(1.0),
    ZMin_(0.0),
    cMax_(1.0),
    cMin_(0.0),
    ZvarMax_(0.25),
    ZvarMin_(0.0),
    cvarMax_(0.25),
    cvarMin_(0.0),
    ZcvarMax_(0.25),
    ZcvarMin_(-0.25),
    TMax_(this->coeffs().lookupOrDefault("TMax", 5000.0)),
    TMin_(this->coeffs().lookupOrDefault("TMin", 200.0)),
    small_FSD_(this->coeffs().lookupOrDefault("small_FSD_", 1.0e-4)),
    dpdt_(this->mesh().objectRegistry::lookupObject<volScalarField>("dpdt")),        
    phi_(this->mesh().objectRegistry::lookupObject<surfaceScalarField>("phi")),
    TCells_(T_.primitiveFieldRef()),
    ignBeginTime_(this->coeffs().lookupOrDefault("ignBeginTime", 0.0)),  
    ignDurationTime_(this->coeffs().lookupOrDefault("ignDurationTime", 0.0)),
    reactFlowTime_(0.0),
    x0_(this->coeffs().lookupOrDefault("x0", 0.0)),   
    y0_(this->coeffs().lookupOrDefault("y0", 0.0)), 
    z0_(this->coeffs().lookupOrDefault("z0", 0.05)),  
    R0_(this->coeffs().lookupOrDefault("R0", 0.04)),
    Sct_(this->coeffs().lookupOrDefault("Sct", 0.7)),
    Sc_(this->coeffs().lookupOrDefault("Sc", 1.0)),   
    Sc_fsd_(this->coeffs().lookupOrDefault("Sc_fsd", 0.7)), 
    beta_FSD_(this->coeffs().lookupOrDefault("beta_FSD", 3.0)),
    phi_fac_FSD_(this->coeffs().lookupOrDefault("phi_fac_FSD", 1.0)),
    bufferTime_(this->coeffs().lookupOrDefault("bufferTime", 0.0)),
    relaxation_(this->coeffs().lookupOrDefault("relaxation", false)),
    curv_SdA_(this->coeffs().lookupOrDefault("curv_SdA", false)),
    DpDt_(this->coeffs().lookupOrDefault("DpDt", false)),
    omega_Yis_(this->coeffs().lookupOrDefault("nOmega_Yis", 1)),  
    magUPrime_(U_.component(vector::X)),
    magUPrimeCells_(magUPrime_.primitiveFieldRef())
{
    if(this->incompPref_ < 0.0)  //local pressure used to calculate EOS
    {
      Info<< "Equation of State: local pressure used" << nl << endl;
    }
    else //constant pressure used to calculate EOS
    {
      Info<< "Equation of State: constant pressure used: "
        << this->incompPref_ << " Pa" << nl << endl;
    }

    omega_YiNames_base_ = wordList();

    if (this->coeffs().found("omega_YiNames"))
    {
        wordList additionalNames = this->coeffs().lookup("omega_YiNames");
        forAll(additionalNames, speciesI)
        {
            omega_YiNames_base_.append(additionalNames[speciesI]);
        }

        Info << "omega_YiNames_base_ is : " << omega_YiNames_base_ <<endl;

        forAll(omega_YiNames_base_, speciesI)
        {
            word specieName2Update = this->omega_YiNames_base_[speciesI];
            const label& specieLabel2Update = this->chemistryPtr_->species()[specieName2Update];
            omega_Yi_index_.append(specieLabel2Update);

            omega_Yis_.set
            (
                speciesI,
                new volScalarField
                (
                    IOobject
                    (
                        "omega_" + this->omega_YiNames_base_[speciesI],
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar("omega_" + this->omega_YiNames_base_[speciesI],dimensionSet(1,-3,-1,0,0,0,0),0.0)
                )
            );
        }

        Info << "end initialize omega_Yi"<<endl;
    }


    //- LES
    this->isLES_ = this->mesh().objectRegistry::foundObject<compressible::LESModel>(turbulenceModel::propertiesName);

    if (this->isLES_)
    {
        const Foam::compressible::LESModel& lesModel = this->mesh().objectRegistry::lookupObject<compressible::LESModel>(turbulenceModel::propertiesName);
            
        const scalarField& LESdeltaCells=lesModel.delta().internalField();

        this->delta_ = lesModel.delta();

        this->deltaCells_ = LESdeltaCells;

    }
    else
    {
        this->deltaCells_ = Foam::pow(this->mesh().V(), 1.0/3.0);
        this->delta_.primitiveFieldRef() = Foam::pow(this->mesh().V(), 1.0/3.0);
    }

    volScalarField::Boundary& deltaBF = this->delta_.boundaryFieldRef();
    forAll(this->delta_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pdelta = deltaBF[patchi];

        forAll(pdelta, facei)
        {
            label cellID = this->mesh().boundary()[patchi].faceCells()[facei];
            pdelta[facei] = this->deltaCells_[cellID];
        }
    }

    // add fields
    fields_.add(Z_);
    fields_.add(Zvar_);
    fields_.add(c_); 
    fields_.add(cvar_);    
    fields_.add(fsd_);  
    fields_.add(He_);  
    fields_.add(Y_fuel_);  

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFSD<ReactionThermo>::~baseFSD()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::baseFSD<ReactionThermo>::transport()
{
    buffer();

    tmp<volScalarField> tmut(this->turbulence().mut());
    const volScalarField& mut = tmut();

    tmp<volScalarField> tnut(this->turbulence().nut());
    const volScalarField& nut = tnut();

    tmp<volScalarField> tnu(this->turbulence().nu());
    const volScalarField& nu = tnu();

    tmp<volScalarField> tmu(this->turbulence().mu());
    const volScalarField& mu = tmu();

    //scalarUW used for cEqn, cvarEqn, ZEqn, ZvarEqn,ZcvarEqn to ensure the convergence when buffer_ is true
    tmp<fv::convectionScheme<scalar> > scalarUWConvection   
    (
        fv::convectionScheme<scalar>::New
        (
            this->mesh(),
            fields_,
            phi_,
            this->mesh().divScheme("div(phi,scalarUW)") 
        )
    );

    // if (this->mesh().objectRegistry::foundObject<sprayCloud>(Foam::sprayCloud::typeName))
    // {
 
    // }
    // else
    {

        // Solve the mixture fraction transport equation
        if(buffer_) 
        {
            Info<< "UW schemes used for scalars" << endl;

            fvScalarMatrix ZEqn
            (
                fvm::ddt(rho_,Z_)
                +scalarUWConvection->fvmDiv(phi_, Z_)  
                -fvm::laplacian( mut/Sct_ + mu/Sc_ , Z_)     
            );
            if(relaxation_)
            {
                ZEqn.relax();
            }
            ZEqn.solve();
        }

        else
        {
            Info<< "TVD schemes used for scalars" << endl;

            fvScalarMatrix ZEqn
            (
                fvm::ddt(rho_,Z_)
                +fvm::div(phi_, Z_)
                -fvm::laplacian(mut/Sct_ + mu/Sc_ , Z_) 
            );

            if(relaxation_)
            {
                ZEqn.relax();
            }
            ZEqn.solve();
        }
        Z_.min(ZMax_);  
        Z_.max(ZMin_);



        // fvScalarMatrix Y_fuelEqn
        // (
        //     fvm::ddt(rho_,Y_fuel_)
        //     +fvm::div(phi_, Y_fuel_)
        //     -fvm::laplacian(mut/Sct_ + mu/Sc_ , Y_fuel_) 
        //     - omega_Y_fuel_
        // );

        // if(relaxation_)
        // {
        //     Y_fuelEqn.relax();
        // }
        // Y_fuelEqn.solve();

        // Y_fuel_.min(ZMax_);  
        // Y_fuel_.max(ZMin_);



        if (combustion_)
        { 

            // Solve the progress variable transport equation
            fvScalarMatrix cEqn
            (
                fvm::ddt(rho_, c_)
                +(
                    buffer_
                    ? scalarUWConvection->fvmDiv(phi_, c_)
                    : fvm::div(phi_, c_)
                )
                -fvm::laplacian( mut/Sct_ + mu/Sc_, c_)
                // - fvm::laplacian( mut/Sct_, c_)
                // - rho_ * SdA_ * fsd_ - omega_c_chiZ_
                -omega_c_
            );  
            
            if(relaxation_)
            {
                cEqn.relax();
            }
            cEqn.solve();
            c_.min(cMax_);
            c_.max(cMin_); 

            c_.correctBoundaryConditions();
        } // end if (combustion_)
    }



    // solve the total enthalpy transport equation
    if(solveEnthalpy_)
    {
        // if (this->mesh().objectRegistry::foundObject<sprayCloud>(Foam::sprayCloud::typeName))
        // {

        // }
        // else
        {
            fvScalarMatrix HEqn
            (
                fvm::ddt(rho_,He_)
                +(
                    buffer_
                    ? scalarUWConvection->fvmDiv(phi_, He_)
                    : fvm::div(phi_, He_)
                )
                +(
                    DpDt_
                    ? - dpdt_ - ( U_ & fvc::grad(p_) ) - fvm::laplacian( mut/Sct_ + mu/Sc_, He_)
                    : - fvm::laplacian( mut/Sct_ + mu/Sc_, He_) 
                )
            ); 

            if(relaxation_)
            {
                HEqn.relax();
            }
            HEqn.solve();
        }
 
    }    

    // Solve the mixture fraction variance transport equation

    fvScalarMatrix ZvarEqn
    (
        fvm::ddt(rho_,Zvar_)
        +(
            buffer_
            ?  scalarUWConvection->fvmDiv(phi_, Zvar_)
            :  fvm::div(phi_, Zvar_)
        )
        -fvm::laplacian( mut/Sct_+mu/Sc_, Zvar_)
        -(2.0*mut/Sct_*(fvc::grad(Z_) & fvc::grad(Z_)))
        +(2.0*rho_*chi_Z_)  
    );

    if(relaxation_)
    {
        ZvarEqn.relax();
    }

    ZvarEqn.solve();

    Zvar_.min(ZvarMax_);
    Zvar_.max(ZvarMin_); 


    // Solve the progress variable variance transport equation
    fvScalarMatrix cvarEqn
    (
        fvm::ddt(rho_,cvar_)
        +(
            buffer_
            ? scalarUWConvection->fvmDiv(phi_, cvar_)
            : fvm::div(phi_, cvar_)
        )
        -fvm::laplacian( mut/Sct_ + mu/Sc_, cvar_)
        -(2.0*mut/Sct_*(fvc::grad(c_) & fvc::grad(c_)))
        +2.0*(rho_*chi_c_)
        -2.0*(cOmega_c_-omega_c_*c_)  
    ); 

    if(relaxation_)
    {
        cvarEqn.relax();
    }
    cvarEqn.solve();
    cvar_.min(cvarMax_);
    cvar_.max(cvarMin_);    


    // Solve the covariance transport equation
    fvScalarMatrix ZcvarEqn
    (
        fvm::ddt(rho_,Zcvar_)
        +(
            buffer_
            ?  scalarUWConvection->fvmDiv(phi_, Zcvar_)
            :  fvm::div(phi_, Zcvar_)
        )
        -fvm::laplacian( mut/Sct_+mu/Sc_, Zcvar_)
        -(2.0*mut/Sct_*(fvc::grad(Z_) & fvc::grad(c_)))
        +(2.0*rho_*chi_Zc_)  
        -1.0*(ZOmega_c_-omega_c_*Z_)  
    );

    if(relaxation_)
    {
        ZcvarEqn.relax();
    }

    ZcvarEqn.solve();

    Zcvar_.min(ZcvarMax_);
    Zcvar_.max(ZcvarMin_);


    // Solve the FSD transport equation
    this->phiFSD_ = linearInterpolate(U_) & this->mesh().Sf();


    if(combustion_ && reactFlowTime_ > 0.0)  
    {
        fvScalarMatrix fsdEqn
        (
            fvm::ddt(fsd_)
            + fvm::div(phiFSD_, fsd_)
            - fvm::laplacian(nut / Sc_fsd_, fsd_)     // dispersion
            - source_fsd_SRE1_ - source_fsd_SRE2_                      // strain-rate
            - source_fsd_flaProp_              // flame propagation
            - source_fsd_FCE1_
            + source_fsd_FCE2_                      // flame curvature
        );

        if(relaxation_)
        {
            fsdEqn.relax();
        }

        fsdEqn.solve();

        fsd_.max(0.0);  
    }

    // forAll(this->rho_, celli)  
    // {
    //     if ( Y_fuelCells_[celli] < 0.0025 )
    //     {
    //         fsdCells_[celli] = 0.0;
    //     }
    // }

    // volScalarField::Boundary& fsdBF = this->fsd_.boundaryFieldRef();
    // volScalarField::Boundary& Y_fuelBF = this->Y_fuel_.boundaryFieldRef();
    // forAll(this->rho_.boundaryField(), patchi)
    // {
    //     fvPatchScalarField& pfsd = fsdBF[patchi];
    //     fvPatchScalarField& pY_fuel = Y_fuelBF[patchi];

    //     forAll(pfsd, facei)
    //     {
    //         if (pY_fuel[facei] < 0.0025)
    //         {
    //             pfsd[facei] = 0.0;
    //         }
    //     }
    // }

}


template<class ReactionThermo>
void Foam::combustionModels::baseFSD<ReactionThermo>::initialiseFalmeKernel()
{

    reactFlowTime_= this->mesh().time().value()-ignBeginTime_; 

    Info<< "Time = " << this->mesh().time().timeName()
            << "   reactFlowTime = " << reactFlowTime_ << nl << endl;

    if(ignition_ && reactFlowTime_ > 0.0 && reactFlowTime_ < ignDurationTime_)
    {
        const vectorField& centre = this->mesh().cellCentres();    
        const scalarField x = centre.component(vector::X);  
        const scalarField y = centre.component(vector::Y); 
        const scalarField z = centre.component(vector::Z); 

        forAll(cCells_,celli) 
        {
            scalar R = Foam::sqrt(magSqr(x[celli]-x0_) + magSqr(y[celli]-y0_)
                                + magSqr(z[celli]-z0_));
            if(R <= R0_) {cCells_[celli] = 1.0;}
        }

        Info<< "Flame initialisation done"<< endl;
    }    
}

template<class ReactionThermo>
void Foam::combustionModels::baseFSD<ReactionThermo>::magUPrime()
{
    autoPtr<LESfilter> filterPtr_(LESfilter::New(this->mesh(), this->turbulence().coeffDict()));
    LESfilter& filter_(filterPtr_());

    magUPrime_=mag(U_-filter_(U_));

   // return magUPrime_;
}

template<class ReactionThermo>
bool Foam::combustionModels::baseFSD<ReactionThermo>::buffer()
{
    if(this->mesh().time().value() - this->mesh().time().startTime().value() < bufferTime_) buffer_ = true;  
    else buffer_ = false;

    return buffer_;
}