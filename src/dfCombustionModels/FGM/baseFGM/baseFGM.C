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

#include "baseFGM.H"
#include "basicSprayCloud.H" 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::baseFGM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    fvOptions(fv::options::New(this->mesh())),
    buffer_(this->coeffs().lookupOrDefault("buffer", false)),
    scaledPV_(this->coeffs().lookupOrDefault("scaledPV", false)),
    incompPref_(this->coeffs().lookupOrDefault("incompPref", -10.0)),
    ignition_(this->coeffs().lookupOrDefault("ignition", false)),
    combustion_(this->coeffs().lookupOrDefault("combustion", false)),
    solveEnthalpy_(this->coeffs().lookupOrDefault("solveEnthalpy", false)),
    flameletT_(this->coeffs().lookupOrDefault("flameletT", false)),
    solveYiEqn_(this->coeffs().lookupOrDefault("solveYiEqn", false)),
    droplet_stat_(this->coeffs().lookupOrDefault("droplet_stat", false)),
    droplet_stat_startTime_(this->coeffs().lookupOrDefault("droplet_stat_startTime", 0.0)),
    tablePath_(this->coeffs().lookup("tablePath")),
    psi_(const_cast<volScalarField&>(dynamic_cast<rhoThermo&>(thermo).psi())),
    alpha_(const_cast<volScalarField&>(thermo.alpha())),
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
    Cp_e_ 
    (
        IOobject
        (
            "Cp_e",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Cp_e",dimensionSet(0,2,-2,-1,0,0,0),1010.1)
    ),
    Cp_mass_ 
    (
        IOobject
        (
            "Cp_mass",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("Cp_mass",dimensionSet(0,2,-2,-1,0,0,0),1010.1)
    ),
    kappa_ 
    (
        IOobject
        (
            "kappa",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("kappa",dimensionSet(1,1,-3,-1,0,0,0),0.223)
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
    c_
    (
        IOobject
        (
            "c",
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
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
    cvar_s_
    (
        IOobject
        (
            "cvar_s",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("cvar_s",dimensionSet(0,0,0,0,0,0,0),0.0) 
    ),    
    Zcvar_s_
    (
        IOobject
        (
            "Zcvar_s",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),  
        dimensionedScalar("Zcvar_s",dimensionSet(0,0,0,0,0,0,0),0.0) 
    ),
    dropletVF_
    (
        IOobject
        (
            "dropletVF",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("dropletVF", dimensionSet(0, 4, -1, 0, 0, 0, 0), 0.0)
    ),
    dropletU_
    (
        IOobject
        (
            "dropletU",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedVector(dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
    ),
    dropletCounter_
    (
        IOobject
        (
            "dropletCounter",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("dropletCounter", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
    ),
    cOmega_c_(omega_c_),
    ZOmega_c_(omega_c_), 
    WtCells_ (Wt_.primitiveFieldRef()),
    Cp_eCells_ (Cp_e_.primitiveFieldRef()),
    Cp_massCells_ (Cp_mass_.primitiveFieldRef()),
    kappaCells_ (kappa_.primitiveFieldRef()),
    ZCells_(Z_.primitiveFieldRef()),
    ZvarCells_(Zvar_.primitiveFieldRef()), 
    HCells_(He_.primitiveFieldRef()),
    HfCells_(Hf_.primitiveFieldRef()),
    cCells_(c_.primitiveFieldRef()),
    cvarCells_(cvar_.primitiveFieldRef()),
    ZcvarCells_(Zcvar_.primitiveFieldRef()), 
    omega_cCells_(omega_c_.primitiveFieldRef()), 
    cOmega_cCells_(cOmega_c_.primitiveFieldRef()),
    ZOmega_cCells_(ZOmega_c_.primitiveFieldRef()),
    chi_cCells_(chi_c_.primitiveFieldRef()),
    chi_ZCells_(chi_Z_.primitiveFieldRef()),
    chi_ZcCells_(chi_Zc_.primitiveFieldRef()),     
    chi_ZfltdCells_(chi_Zfltd_.primitiveFieldRef()),  
    betacCells_(betac_.primitiveFieldRef()),     
    ZMax_(1.0),
    ZMin_(0.0),
    ZvarMax_(0.25),
    ZvarMin_(0.0),
    cMax_(1.0),
    cMin_(0.0),
    Ycmaxall_(1.0),
    cvarMax_(0.25),
    cvarMin_(0.0),
    ZcvarMax_(0.25),
    ZcvarMin_(-0.25),
    TMax_(this->coeffs().lookupOrDefault("TMax", 5000.0)),
    TMin_(this->coeffs().lookupOrDefault("TMin", 200.0)),
    rho_(const_cast<volScalarField&>(this->mesh().objectRegistry::lookupObject<volScalarField>("rho"))),
    rho_inRhoThermo_(dynamic_cast<rhoThermo&>(thermo).rho()),
    p_(this->thermo().p()),
    T_(this->thermo().T()),
    Y_(this->chemistryPtr_->Y()),
    U_(this->mesh().objectRegistry::lookupObject<volVectorField>("U")),
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
    omega_Yis_(this->coeffs().lookupOrDefault("nOmega_Yis", 1)),
    bufferTime_(this->coeffs().lookupOrDefault("bufferTime", 0.0)),
    relaxation_(this->coeffs().lookupOrDefault("relaxation", false)),
    DpDt_(this->coeffs().lookupOrDefault("DpDt", false)),
    magUPrime_(U_.component(vector::X)),
    magUPrimeCells_(magUPrime_.primitiveFieldRef())
{
    if(incompPref_ < 0.0)  //local pressure used to calculate EOS
    {
      Info<< "Equation of State: local pressure used" << nl << endl;
    }
    else //constant pressure used to calculate EOS
    {
      Info<< "Equation of State: constant pressure used: "
        << incompPref_ << " Pa" << nl << endl;
    }


    // this->omega_YiNames_base_ = wordList({"C2H5OH", "O2"});
    // this->omega_YiNames_base_ = wordList({"C2H5OH"});
    // this->omega_YiNames_base_ = wordList({"CH4"});
    // Info << "omega_YiNames_base_ is : " << omega_YiNames_base_ <<endl;
    
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

        this->deltaCells_ = LESdeltaCells;

    }


//     //initialize species fields
//     forAll(Y_, speciesI)
//     {

//      Y_.set
//      (
//       speciesI,
//         new volScalarField
//         (
//             IOobject
//             (
//                 "Y_"  + speciesNames_[speciesI],
//                 this->mesh().time().timeName(),
//                 this->mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             this->mesh(),
//             dimensionedScalar(dimless, Zero)
//         )
//      );

//    }

   // add fields
     fields_.add(Z_);
     fields_.add(Zvar_);
     fields_.add(c_);     
     fields_.add(cvar_);  
     fields_.add(Zcvar_);  
     fields_.add(He_);  

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::~baseFGM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::transport()
{
    buffer();

    tmp<volScalarField> tmut(this->turbulence().mut());
    const volScalarField& mut = tmut();

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

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            this->mesh(),
            fields_,
            phi_,
            this->mesh().divScheme("div(phi,Yi_h)")
        )
    );


    if (this->mesh().objectRegistry::foundObject<sprayCloud>(Foam::sprayCloud::typeName))
    {
        Info << "Solve Spray Z" << endl;
 
        //- create spray object 
        basicSprayCloud& spray = dynamic_cast<basicSprayCloud&>
        (this->mesh().objectRegistry::lookupObjectRef<sprayCloud>(Foam::sprayCloud::typeName));        

        tmp<fvScalarMatrix> tZSource(new fvScalarMatrix(this->Z_, dimMass/dimTime));
        fvScalarMatrix& SZ = tZSource.ref();
        forAll(spray.rhoTrans(), i)
        {
            SZ.source() += spray.rhoTrans()[i]/this->mesh().time().deltaT();
            spray.rhoTrans()[i].writeOpt() = IOobject::NO_WRITE;
        }

        // Solve the mixture fraction transport equation
        if(buffer_) 
        {
            Info<< "UW schemes used for scalars" << endl;

            fvScalarMatrix ZEqn
            (
                fvm::ddt(rho_,Z_)
                +scalarUWConvection->fvmDiv(phi_, Z_)  
                -fvm::laplacian( mut/Sct_ + mu/Sc_ , Z_)   
                ==
                - SZ  
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
                ==
                - SZ                  
            );

            if(relaxation_)
            {
                ZEqn.relax();
            }
            ZEqn.solve();
        }
        Z_.min(ZMax_);  
        Z_.max(ZMin_);

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


        if (combustion_ && reactFlowTime_ > 0.0)
        {

            tmp<fvScalarMatrix> tcSource(new fvScalarMatrix(this->c_, dimMass/dimTime));
            fvScalarMatrix& S_c = tcSource.ref();
            forAll(spray.rhoTrans(), i)
            {
                S_c.source() += spray.rhoTrans()[i]/this->mesh().time().deltaT();
                spray.rhoTrans()[i].writeOpt() = IOobject::NO_WRITE;
            }

            // At initial time, cMax is set as Ycmaxall when unscaled PV employed
            if(this->mesh().time().timeIndex() == 1 && !scaledPV_) cMax_ = Ycmaxall_;  

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
                -omega_c_
                // + fvOptions(rho_, c_)
                ==
                c_ * S_c
            );  
            // fvOptions.correct(c_);
            
            if(relaxation_)
            {
                cEqn.relax();
            }
            cEqn.solve();
            c_.min(cMax_);
            c_.max(cMin_); 


            // Solve the progress variable variance transport equation
            tmp<fvScalarMatrix> tcvarSource(new fvScalarMatrix(this->cvar_, dimMass/dimTime));
            fvScalarMatrix& S_cvar = tcvarSource.ref();
            forAll(spray.rhoTrans(), i)
            {
                S_cvar.source() += spray.rhoTrans()[i]/this->mesh().time().deltaT();
                spray.rhoTrans()[i].writeOpt() = IOobject::NO_WRITE;
            }
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
                ==
                2.0 * cvar_ * S_cvar
            ); 

            if(relaxation_)
            {
                cvarEqn.relax();
            }
            cvarEqn.solve();
            cvar_.min(cvarMax_);
            cvar_.max(cvarMin_);      


            // Solve the covariance transport equation
            tmp<fvScalarMatrix> tZcvarSource(new fvScalarMatrix(this->Zcvar_, dimMass/dimTime));
            fvScalarMatrix& S_Zcvar = tZcvarSource.ref();
            forAll(spray.rhoTrans(), i)
            {
                S_Zcvar.source() += spray.rhoTrans()[i]/this->mesh().time().deltaT();
                spray.rhoTrans()[i].writeOpt() = IOobject::NO_WRITE;
            }
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
                ==
                Zcvar_ * S_Zcvar
            );

            if(relaxation_)
            {
                ZcvarEqn.relax();
            }

            ZcvarEqn.solve();

            Zcvar_.min(ZcvarMax_);
            Zcvar_.max(ZcvarMin_);
        } // end if (combustion_)
        

        // solve the total enthalpy transport equation
        if(solveEnthalpy_)
        {
            Info << "Solve Spray He" << endl;
            //- create spray object 
            basicSprayCloud& spray = dynamic_cast<basicSprayCloud&>
            (this->mesh().objectRegistry::lookupObjectRef<sprayCloud>(Foam::sprayCloud::typeName));     
                      
            tmp<fvScalarMatrix> thSource(new fvScalarMatrix(this->He_, dimEnergy/dimTime));
            fvScalarMatrix& hSource = thSource.ref();

            forAll(this->chemistryPtr_->species(), i)
            {
                hSource.source() += spray.rhoTrans(i)*this->chemistryPtr_->hci(i)/this->mesh().time().deltaT();
            }

            hSource.source() += spray.hsTrans()/this->mesh().time().deltaT();  // equivalent to parcels.Sh(he)


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
                ==
                - hSource
            ); 

            if(relaxation_)
            {
                HEqn.relax();
            }
            HEqn.solve();
        }

        if ( (!this->omega_YiNames_base_.empty()) 
            and this->solveYiEqn_ )
        {

            forAll(omega_Yi_index_, i)
            {
                volScalarField& Yi = Y_[omega_Yi_index_[i]];

                Info << "Solve Spray Y_" << omega_YiNames_base_[i] << endl;

                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho_, Yi)
                    + (
                        buffer_
                        ? scalarUWConvection->fvmDiv(phi_, Yi)
                        : mvConvection->fvmDiv(phi_, Yi)
                    )
                    -fvm::laplacian( mut/Sct_ + mu/Sc_, Yi)
                    - omega_Yis_[i]
                    + fvOptions(rho_, Yi)
                    ==
                    spray.SYi(omega_Yi_index_[i], Yi)  
                );  
                fvOptions.correct(Yi);
                
                if(relaxation_)
                {
                    YiEqn.relax();
                }
                YiEqn.solve("Yi");

                Yi.min(1.0);
                Yi.max(0.0); 
            }
                    
        }

        if (droplet_stat_ && this->mesh().time().value() >= droplet_stat_startTime_)
        {
            Info << "droplet statistic"<< endl;
            forAllConstIter(basicSprayCloud, spray, iter)
            {
                const auto& droplet = iter();
                label cellId = spray.mesh().findCell(droplet.position());

                dropletCounter_[cellId] += 1.0;
                scalar beta = 1.0 / dropletCounter_[cellId];

                dropletVF_[cellId] = (1.0 - beta)*dropletVF_[cellId] 
                                    + beta* ( M_PI / 6.0 * pow3(droplet.d()) * mag(droplet.U()) );
                dropletU_[cellId] = (1.0 - beta)*dropletU_[cellId] + beta* droplet.U();

                // std::cout << "droplet statistic: dropletU = " << droplet.U()[0] << ", d = "<< droplet.d() << std::endl;
            }

            volScalarField::Boundary& dropletVFBF = this->dropletVF_.boundaryFieldRef();
            volVectorField::Boundary& dropletUBF = this->dropletU_.boundaryFieldRef();
            forAll(this->dropletVF_.boundaryFieldRef(), patchi)
            {
                fvPatchScalarField& pdropletVF = dropletVFBF[patchi];
                fvPatchVectorField& pdropletU = dropletUBF[patchi];

                forAll(pdropletVF, facei)
                {
                    label cellID = this->mesh().boundary()[patchi].faceCells()[facei];
                    pdropletVF[facei] = this->dropletVF_[cellID];
                    pdropletU[facei] = this->dropletU_[cellID];
                }
            }
        }
    }
    else  // not spray
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

        if(combustion_ && reactFlowTime_ > 0.0)
        {

            // At initial time, cMax is set as Ycmaxall when unscaled PV employed
            if(this->mesh().time().timeIndex() == 1 && !scaledPV_) cMax_ = Ycmaxall_;  

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
                -omega_c_
                // + fvOptions(rho_, c_)
            );  
            // fvOptions.correct(c_);
            
            if(relaxation_)
            {
                cEqn.relax();
            }
            cEqn.solve();
            c_.min(cMax_);
            c_.max(cMin_); 


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
        } // end if (combustion_)

        if(solveEnthalpy_)
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


        if ( (!this->omega_YiNames_base_.empty()) 
            and this->solveYiEqn_ )
        {

            forAll(omega_Yi_index_, i)
            {
                volScalarField& Yi = Y_[omega_Yi_index_[i]];

                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho_, Yi)
                    + (
                        buffer_
                        ? scalarUWConvection->fvmDiv(phi_, Yi)
                        : mvConvection->fvmDiv(phi_, Yi)
                    )
                    -fvm::laplacian( mut/Sct_ + mu/Sc_, Yi)
                    - omega_Yis_[i]
                    + fvOptions(rho_, Yi)
                );  
                fvOptions.correct(Yi);
                
                if(relaxation_)
                {
                    YiEqn.relax();
                }
                YiEqn.solve("Yi");

                Yi.min(1.0);
                Yi.max(0.0); 
            }
                    
        }
    }


    if (this->isLES_)
    {
        chi_Zfltd_ = this->turbulence().mu()/rho_*(fvc::grad(Z_) & fvc::grad(Z_)); 
    }

}

template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::initialiseFalmeKernel()
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
void Foam::combustionModels::baseFGM<ReactionThermo>::magUPrime()
{
    autoPtr<LESfilter> filterPtr_(LESfilter::New(this->mesh(), this->turbulence().coeffDict()));
    LESfilter& filter_(filterPtr_());

    magUPrime_=mag(U_-filter_(U_));

   // return magUPrime_;
}


template<class ReactionThermo>
bool Foam::combustionModels::baseFGM<ReactionThermo>::buffer()
{
    if(this->mesh().time().value() - this->mesh().time().startTime().value() < bufferTime_) buffer_ = true;  
    else buffer_ = false;

    return buffer_;
}


// ************************************************************************* //
