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

#include "FSD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FSD<ReactionThermo>::FSD
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    baseFSD<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    tableSolver_FSD(baseFSD<ReactionThermo>::tablePath_)
{
    //- retrieval data from table
    retrieval();

    Info<< "At FSD, min/max(T) = " << min(this->T_).value() << ", " << max(this->T_).value() << endl;    

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::FSD<ReactionThermo>::~FSD()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::FSD<ReactionThermo>::correct()
{
    //- initialize flame kernel
    baseFSD<ReactionThermo>::initialiseFalmeKernel();

    //- solve transport equation
    baseFSD<ReactionThermo>::transport();

    //update enthalpy using lookup data
    if(!(this->solveEnthalpy_))
    {
        this->He_ = this->Z_*(this->H_fuel-this->H_ox) + this->H_ox;
    }

    if (this->isLES_)
    {
        baseFSD<ReactionThermo>::magUPrime();
    }
  
    //- retrieval data from table
    retrieval();
}

template<class ReactionThermo>
void Foam::combustionModels::FSD<ReactionThermo>::retrieval()
{
    tmp<volScalarField> tk(this->turbulence().k());
    volScalarField& k = const_cast<volScalarField&>(tk());
    scalarField& kCells = k.primitiveFieldRef();

    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    volScalarField& epsilon = const_cast<volScalarField&>(tepsilon());
    const scalarField& epsilonCells =epsilon.primitiveFieldRef();

    tmp<volScalarField> tmu = this->turbulence().mu();  
    volScalarField& mu = const_cast<volScalarField&>(tmu());
    scalarField& muCells = mu.primitiveFieldRef();     

    tmp<volScalarField> tmut = this->turbulence().mut();  
    volScalarField& mut = const_cast<volScalarField&>(tmut());
    scalarField& mutCells = mut.primitiveFieldRef();   

    // for FSD source term calculation
    volVectorField gradC_ = fvc::grad(this->c_);
    vectorField& gradCCells = gradC_.primitiveFieldRef(); 

    volScalarField gradC_mod_ = Foam::mag(gradC_);  
    scalarField& gradC_modCells = gradC_mod_.primitiveFieldRef(); 

    volScalarField u_fluct_ = Foam::sqrt(2.0*k/3.0); 
    scalarField& u_fluctCells = u_fluct_.primitiveFieldRef();

    this->chi_Zc_ = mu/(this->Sc_ * this->rho_) * ( fvc::grad(this->Z_) & fvc::grad(this->c_) );

    scalar Zl, Zr, ZKl, ZKr;

    scalar delta_max = Foam::max(this->deltaCells_);

    forAll(this->rho_, celli)  
    {
        if(this->isLES_)
        {
            this->chi_ZCells_[celli] = this->sdrLRXmodel(2.0,mutCells[celli]
                /this->rho_[celli],this->deltaCells_[celli],this->ZvarCells_[celli]); 
        }
        else
        {
            this->chi_ZCells_[celli] = 1.0*epsilonCells[celli]/kCells[celli] *this->ZvarCells_[celli]; 
        }


        // look up phiTab
        Zl = this->z_Tb5[0];
        Zr = this->z_Tb5[this->NZL-1];

        if(this->ZCells_[celli] >= Zl && this->ZCells_[celli] <= Zr
            && this->combustion_
            // && this->cCells_[celli] > this->small
            ) 
        {
            this->SL0_FSDCells_[celli] = this->lookup1d(this->NZL, this->z_Tb5, this->ZCells_[celli], this->sl_Tb5);
            this->thermoThickness_L0_FSDCells_[celli] = this->lookup1d(this->NZL, this->z_Tb5, this->ZCells_[celli], this->th_Tb5);
            this->rho_uCells_[celli] = this->lookup1d(this->NZL, this->z_Tb5, this->ZCells_[celli], this->rho_u_Tb5);
        }
        else
        {
            this->SL0_FSDCells_[celli] = 0.0;
            this->thermoThickness_L0_FSDCells_[celli] = 0.0;
            this->rho_uCells_[celli] = this->rho_[celli];
        }

        if (this->SL0_FSDCells_[celli] > this->small && this->thermoThickness_L0_FSDCells_[celli] > this->small
            && this->fsdCells_[celli] > 1.0 && this->cCells_[celli] > this->small
            )
        {
            this->gamma_kCells_[celli] = 0.75 * Foam::exp( -1.2 / Foam::pow(u_fluctCells[celli] / this->SL0_FSDCells_[celli], 0.3) )
                                * Foam::pow(5.0 * (this->deltaCells_[celli]) / this->thermoThickness_L0_FSDCells_[celli], 2.0/3.0);
        }
        else
        {
            this->gamma_kCells_[celli] = 0.0;
        }

        double gz{this->cal_gvar(this->ZCells_[celli],this->ZvarCells_[celli])}, Ycmax{-1.0};   
        double gc = this->cal_gvar(this->cCells_[celli],this->cvarCells_[celli],Ycmax);  

        double cNorm = this->cCells_[celli];

        if(this->ZCells_[celli] >= Zl && this->ZCells_[celli] <= Zr
            && this->combustion_ && this->cCells_[celli] > this->small)  
        {
            double kc_s = this->lookup1d(this->NZL,this->z_Tb5,this->ZCells_[celli],this->kctau_Tb5);      
            double tau = this->lookup1d(this->NZL,this->z_Tb5,this->ZCells_[celli],this->tau_Tb5);    
            double sl = this->lookup1d(this->NZL,this->z_Tb5,this->ZCells_[celli],this->sl_Tb5);     
            double dl = this->lookup1d(this->NZL,this->z_Tb5,this->ZCells_[celli],this->th_Tb5);   

            if (this->isLES_)
            {
                this->chi_cCells_[celli] = this->sdrFLRmodel(this->cvarCells_[celli],this->magUPrimeCells_[celli],
                            this->deltaCells_[celli],sl,dl,tau,kc_s,this->betacCells_[celli]);

                this->omega_cCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                                    this->NC, this->c_Tb3, this->cCells_[celli],
                                                    this->NGZ, this->gz_Tb3, gz,
                                                    this->NGC, this->gc_Tb3, gc,
                                                    this->tableValues_[0])
                                            + (this->chi_ZCells_[celli] + this->chi_ZfltdCells_[celli])*this->cCells_[celli]
                                                *this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2);

                this->omega_c_chiZCells_[celli] = (this->chi_ZCells_[celli] + this->chi_ZfltdCells_[celli])*this->cCells_[celli]
                                                *this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2)
                                                + 2.0*this->chi_ZcCells_[celli]*this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d1Yeq_Tb2);
            }
            else
            {
                this->chi_cCells_[celli] =
                    this->RANSsdrFLRmodel(this->cvarCells_[celli],epsilonCells[celli],
                        kCells[celli],muCells[celli]/this->rho_[celli],
                        sl,dl,tau,kc_s,this->rho_[celli]);  

                this->omega_cCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                                            this->NC, this->c_Tb3, this->cCells_[celli],
                                                            this->NGZ, this->gz_Tb3, gz,
                                                            this->NGC, this->gc_Tb3, gc,
                                                            this->tableValues_[0])
                                            + this->chi_ZCells_[celli]*this->cCells_[celli]
                                                *this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2);

                this->omega_c_chiZCells_[celli] = this->chi_ZCells_[celli]*this->cCells_[celli]
                                                *this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2)
                                                + 2.0*this->chi_ZcCells_[celli]*this->lookup2d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                                this->NGZ,this->gz_Tb3,gz,this->d1Yeq_Tb2);
            }

            this->cOmega_cCells_[celli] = this->lookup4d(this->NZ,this->z_Tb3,this->ZCells_[celli],
                                                this->NC,this->c_Tb3,cNorm,
                                                this->NGZ,this->gz_Tb3,gz,
                                                this->NGC,this->gc_Tb3,gc,
                                                this->tableValues_[1]);    

            this->omega_Y_fuelCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                    this->NC, this->c_Tb3, this->cCells_[celli],
                                    this->NGZ, this->gz_Tb3, gz,
                                    this->NGC, this->gc_Tb3, gc,
                                    this->tableValues_[7]);
        }           
        else
        {
            if(this->isLES_)
            {
                this->chi_cCells_[celli] = this->sdrLRXmodel(2.0,mutCells[celli]
                           /this->rho_[celli],this->deltaCells_[celli],this->cvarCells_[celli]);     
            }
            else
            {
                this->chi_cCells_[celli] = 1.0*epsilonCells[celli]/kCells[celli]*this->cvarCells_[celli]; 
            }

            this->omega_cCells_[celli] = 0.0;
            this->cOmega_cCells_[celli] = 0.0;
            this->omega_c_chiZCells_[celli] = 0.0;
            this->omega_Y_fuelCells_[celli] = 0.0;
        }
        this->omega_cCells_[celli] = this->omega_cCells_[celli]*this->rho_[celli];  
        this->omega_c_chiZCells_[celli] = this->omega_c_chiZCells_[celli] *this->rho_[celli];  
        this->cOmega_cCells_[celli] = this->cOmega_cCells_[celli] *this->rho_[celli];  
        this->omega_Y_fuelCells_[celli] = this->omega_Y_fuelCells_[celli]*this->rho_[celli];  


        this->Zvar_s_[celli] = gz;

        this->CpCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                        this->NC, this->c_Tb3, this->cCells_[celli],
                                        this->NGZ, this->gz_Tb3, gz,
                                        this->NGC, this->gc_Tb3, gc,
                                        this->tableValues_[2]);

        this->HfCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                        this->NC, this->c_Tb3, this->cCells_[celli],
                                        this->NGZ, this->gz_Tb3, gz,
                                        this->NGC, this->gc_Tb3, gc,
                                        this->tableValues_[3]);

        this->WtCells_[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                        this->NC, this->c_Tb3, this->cCells_[celli],
                                        this->NGZ, this->gz_Tb3, gz,
                                        this->NGC, this->gc_Tb3, gc,
                                        this->tableValues_[4]);

        muCells[celli] = this->lookup4d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                    this->NC, this->c_Tb3, this->cCells_[celli],
                                    this->NGZ, this->gz_Tb3, gz,
                                    this->NGC, this->gc_Tb3, gc,
                                    this->tableValues_[5])*this->rho_[celli];          

        if (baseFSD<ReactionThermo>::flameletT_)
        {
            this->TCells_[celli] = this->lookup3d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                            this->NC, this->c_Tb3, this->cCells_[celli],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->tableValues_[6]);
        }
        else
        {
            this->TCells_[celli] = (this->He_[celli]-this->Hf_[celli])/this->Cp_[celli] + this->T0;
        }

        this->TCells_[celli] = max(this->TCells_[celli], this->TMin_);
        this->TCells_[celli] = min(this->TCells_[celli], this->TMax_);


        // // -------------------- Yis begin ------------------------------
        for (int yi=0; yi<this->NY; yi++)
        {
            word specieName2Update = this->speciesNames_table_[yi];
            const label& specieLabel2Update = this->chemistryPtr_->species()[specieName2Update];
            this->Y_[specieLabel2Update].primitiveFieldRef()[celli] = this->lookup3d(this->NZ, this->z_Tb3, this->ZCells_[celli],
                                            this->NC, this->c_Tb3, this->cCells_[celli],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->tableValues_[this->NS+yi]);  
        }

        // -------------------- Yis end ------------------------------

        //- surface-averaged normal vector
        if (
            // gradC_modCells[celli] > 0.05 / delta_max
            gradC_modCells[celli] <= this->fsdCells_[celli]
            && this->omega_cCells_[celli] > 1.0
            // && this->cCells_[celli] > this->small 
            && this->fsdCells_[celli] > 0.05 / delta_max
            && this->ZCells_[celli] >= 0.99*Zl && this->ZCells_[celli] <= 1.01*Zr
            && this->combustion_ && this->cCells_[celli] > this->small)
        {
            this->n_FSDCells_[celli] = - gradCCells[celli] / (this->fsdCells_[celli]);
        }
        else if (gradC_modCells[celli] > 0.05 / delta_max
            && this->omega_cCells_[celli] > 1.0
            // && this->fsdCells_[celli] > (0.05 - this->small) / delta_max
            && this->ZCells_[celli] >= 0.99*Zl && this->ZCells_[celli] <= 1.01*Zr
            && this->combustion_ && this->cCells_[celli] > this->small)
        {
            this->n_FSDCells_[celli] = - gradCCells[celli] / gradC_modCells[celli];
        }
        else
        {
            this->n_FSD_.primitiveFieldRef()[celli] = vector::zero;
        }

    }

    //----------------------------- update boundary ---------------------------------//

    volScalarField::Boundary& rhoBF = this->rho_.boundaryFieldRef();
    volScalarField::Boundary& HeBF = this->He_.boundaryFieldRef();
    volScalarField::Boundary& ZBF = this->Z_.boundaryFieldRef();
    volScalarField::Boundary& cBF = this->c_.boundaryFieldRef();
    volScalarField::Boundary& ZvarBF = this->Zvar_.boundaryFieldRef();
    volScalarField::Boundary& cvarBF = this->cvar_.boundaryFieldRef();
    volScalarField::Boundary& KaBF = this->Ka_.boundaryFieldRef();
    volScalarField::Boundary& chi_ZBF = this->chi_Z_.boundaryFieldRef();
    volScalarField::Boundary& chi_cBF = this->chi_c_.boundaryFieldRef();
    volScalarField::Boundary& chi_ZcBF = this->chi_Zc_.boundaryFieldRef();
    volScalarField::Boundary& chi_ZfltdBF = this->chi_Zfltd_.boundaryFieldRef();
    volScalarField::Boundary& muBF = mu.boundaryFieldRef();
    volScalarField::Boundary& mutBF = mut.boundaryFieldRef();
    volScalarField::Boundary& kBF = k.boundaryFieldRef();
    volScalarField::Boundary& epsilonBF = epsilon.boundaryFieldRef();
    volScalarField::Boundary& I_s_FSDBF = this->I_s_FSD_.boundaryFieldRef();
    volScalarField::Boundary& rho_uBF = this->rho_u_.boundaryFieldRef();
    volScalarField::Boundary& omega_cBF = this->omega_c_.boundaryFieldRef();
    volScalarField::Boundary& cOmega_cBF = this->cOmega_c_.boundaryFieldRef();
    volScalarField::Boundary& omega_c_chiZBF = this->omega_c_chiZ_.boundaryFieldRef();
    volScalarField::Boundary& omega_Y_fuelBF = this->omega_Y_fuel_.boundaryFieldRef();
    volScalarField::Boundary& CpBF = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& HfBF = this->Hf_.boundaryFieldRef();
    volScalarField::Boundary& WtBF = this->Wt_.boundaryFieldRef();
    volScalarField::Boundary& SL0_FSDBF = this->SL0_FSD_.boundaryFieldRef();
    volScalarField::Boundary& thermoThickness_L0_FSDBF = this->thermoThickness_L0_FSD_.boundaryFieldRef();
    volScalarField::Boundary& TBF = this->T_.boundaryFieldRef();
    volScalarField::Boundary& fsdBF = this->fsd_.boundaryFieldRef();
    volScalarField::Boundary& u_fluctBF = u_fluct_.boundaryFieldRef();
    volVectorField::Boundary& n_FSDBF = this->n_FSD_.boundaryFieldRef();
    volVectorField::Boundary& grad_CBF = gradC_.boundaryFieldRef();
    volScalarField::Boundary& gamma_kBF = this->gamma_k_.boundaryFieldRef();
    volScalarField::Boundary& gradC_modBF = gradC_mod_.boundaryFieldRef();
    volScalarField::Boundary& magUPrimeBF = this->magUPrime_.boundaryFieldRef();


    forAll(this->rho_.boundaryField(), patchi)
    {
        fvPatchScalarField& prho = rhoBF[patchi];
        fvPatchScalarField& pHe = HeBF[patchi];
        fvPatchScalarField& pZ = ZBF[patchi];
        fvPatchScalarField& pc = cBF[patchi];
        fvPatchScalarField& pZvar = ZvarBF[patchi];
        fvPatchScalarField& pcvar = cvarBF[patchi];
        fvPatchScalarField& pchi_Z = chi_ZBF[patchi];
        fvPatchScalarField& pchi_c = chi_cBF[patchi];
        fvPatchScalarField& pchi_Zc = chi_ZcBF[patchi];
        fvPatchScalarField& pchi_Zfltd = chi_ZfltdBF[patchi];
        fvPatchScalarField& pmu = muBF[patchi];
        fvPatchScalarField& pmut = mutBF[patchi];
        fvPatchScalarField& pk = kBF[patchi];
        fvPatchScalarField& pepsilon = epsilonBF[patchi];
        // fvPatchScalarField& pI_s_FSD = I_s_FSDBF[patchi];
        fvPatchScalarField& prho_u = rho_uBF[patchi];
        fvPatchScalarField& pomega_c = omega_cBF[patchi];
        fvPatchScalarField& pcOmega_c = cOmega_cBF[patchi];
        fvPatchScalarField& pomega_c_chiZ = omega_c_chiZBF[patchi];
        fvPatchScalarField& pomega_Y_fuel = omega_Y_fuelBF[patchi];
        fvPatchScalarField& pCp = CpBF[patchi];
        fvPatchScalarField& pHf = HfBF[patchi];
        fvPatchScalarField& pWt = WtBF[patchi];
        fvPatchScalarField& pSL0_FSD = SL0_FSDBF[patchi];
        fvPatchScalarField& pthermoThickness_L0_FSD = thermoThickness_L0_FSDBF[patchi];
        fvPatchScalarField& pT = TBF[patchi];
        fvPatchScalarField& pfsd = fsdBF[patchi];
        fvPatchScalarField& pu_fluct = u_fluctBF[patchi];
        fvPatchVectorField& pn_FSD = n_FSDBF[patchi];
        fvPatchVectorField& pgrad_C = grad_CBF[patchi];
        fvPatchScalarField& pgamma_k = gamma_kBF[patchi];
        fvPatchScalarField& pgradC_mod = gradC_modBF[patchi];
        fvPatchScalarField& pmagUPrime = magUPrimeBF[patchi];

        forAll(prho, facei)
        {
            label cellID = this->mesh().boundary()[patchi].faceCells()[facei]; 

            if(this->isLES_)
            {
                pchi_Z[facei] = sdrLRXmodel(2.0,pmut[facei]
                    /prho[facei],this->deltaCells_[cellID],pZvar[facei]);
            }
            else
            {
                pchi_Z[facei] = 1.0*pepsilon[facei]/pk[facei] *pZvar[facei];
            }

            Zl = this->z_Tb5[0];
            Zr = this->z_Tb5[this->NZL-1];

            if(pZ[facei] >= Zl && pZ[facei] <= Zr
                && this->combustion_ 
                // && pc[facei] > this->small
                )
            {   
                pSL0_FSD[facei] = this->lookup1d(this->NZL, this->z_Tb5, pZ[facei], this->sl_Tb5);
                pthermoThickness_L0_FSD[facei] = this->lookup1d(this->NZL, this->z_Tb5, pZ[facei], this->th_Tb5);
                prho_u[facei] = this->lookup1d(this->NZL, this->z_Tb5, pZ[facei], this->rho_u_Tb5);
            }
            else
            {
                pSL0_FSD[facei] = 0.0;
                pthermoThickness_L0_FSD[facei] = 0.0;
                prho_u[facei] = prho[facei];
                // if (pZ[facei] < Zl)
                // {
                //     prho_u[facei] = this->rho_u_Tb5[0];
                // }
                // else if (pZ[facei] > Zr)
                // {
                //     prho_u[facei] = this->rho_u_Tb5[this->NZL-1];
                // }
                // else
                // {
                //     prho_u[facei] = prho[facei];
                // }
            }

            if (pSL0_FSD[facei] > this->small && pthermoThickness_L0_FSD[facei] > this->small
                && pfsd[facei] > 1.0 && pc[facei] > this->small
                )
            {
                pgamma_k[facei] = 0.75 * Foam::exp( -1.2 / Foam::pow(pu_fluct[facei] / pSL0_FSD[facei], 0.3) )
                                * Foam::pow(5.0 * (this->deltaCells_[cellID]) / pthermoThickness_L0_FSD[facei], 2.0/3.0);
            }
            else
            {
                pgamma_k[facei] = 0.0;
            }

            double cNorm = pc[facei];

            double gz{this->cal_gvar(pZ[facei],pZvar[facei])}, Ycmax{-1.0}; 
            double gc = this->cal_gvar(pc[facei],pcvar[facei],Ycmax);   

            if(pZ[facei] >= Zl && pZ[facei] <= Zr
                && this->combustion_ && pc[facei] > this->small) 
            {
                double kc_s = this->lookup1d(this->NZL,this->z_Tb5,pZ[facei],this->kctau_Tb5);     
                double tau = this->lookup1d(this->NZL,this->z_Tb5,pZ[facei],this->tau_Tb5);    
                double sl = this->lookup1d(this->NZL,this->z_Tb5,pZ[facei],this->sl_Tb5);      
                double dl = this->lookup1d(this->NZL,this->z_Tb5,pZ[facei],this->th_Tb5);      

                if(this->isLES_)
                {
                    pchi_c[facei] =  sdrFLRmodel(pcvar[facei],pmagUPrime[facei],
                                      this->deltaCells_[cellID],sl,dl,tau,kc_s,this->betacCells_[cellID]);   
                }
                else
                {
                    pchi_c[facei] =
                        this->RANSsdrFLRmodel(pcvar[facei],pepsilon[facei],
                            pk[facei],pmu[facei]/prho[facei],
                            sl,dl,tau,kc_s,prho[facei]);
                }


                if (this->isLES_)
                {
                    pomega_c[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[0])
                                    + (pchi_Z[facei]+ pchi_Zfltd[facei])*pc[facei]
                                    *this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                    this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2);

                    pomega_c_chiZ[facei] = (pchi_Z[facei]+ pchi_Zfltd[facei])*pc[facei]
                                    *this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                    this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2)
                                    + 2.0*pchi_Zc[facei]*this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                    this->NGZ,this->gz_Tb3,gz,this->d1Yeq_Tb2);
                }
                else
                {
                    pomega_c[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[0])
                                    + pchi_Z[facei]*pc[facei]
                                *this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2);

                    pomega_c_chiZ[facei] = pchi_Z[facei]*pc[facei]
                                *this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                this->NGZ,this->gz_Tb3,gz,this->d2Yeq_Tb2)
                                + 2.0*pchi_Zc[facei]*this->lookup2d(this->NZ,this->z_Tb3,pZ[facei],
                                                    this->NGZ,this->gz_Tb3,gz,this->d1Yeq_Tb2);
                }

                pcOmega_c[facei] = this->lookup4d(this->NZ,this->z_Tb3,pZ[facei],
                                    this->NC,this->c_Tb3,cNorm,
                                    this->NGZ,this->gz_Tb3,gz,
                                    this->NGC,this->gc_Tb3,gc,
                                    this->tableValues_[1]);

                pomega_Y_fuel[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                this->NC, this->c_Tb3, pc[facei],
                                this->NGZ, this->gz_Tb3, gz,
                                this->NGC,this->gc_Tb3,gc,
                                this->tableValues_[7]);
            }
            else
            {
                if(this->isLES_)
                {
                    pchi_c[facei] = sdrLRXmodel(2.0,pmut[facei]
                        /this->rho_[facei],this->deltaCells_[cellID],pcvar[facei]);                    
                }
                else
                {
                    pchi_c[facei] = 1.0*pepsilon[facei]/(pk[facei]+SMALL)*pcvar[facei];    
                }

                pomega_c[facei] = 0.0;
                pcOmega_c[facei] = 0.0;
                pomega_Y_fuel[facei] = 0.0;
                pomega_c_chiZ[facei] = 0.0;
            }
            pomega_c[facei] = pomega_c[facei]*prho[facei];  
            pcOmega_c[facei] = pcOmega_c[facei] * prho[facei];  
            pomega_c_chiZ[facei] = pomega_c_chiZ[facei]*prho[facei];  
            pomega_Y_fuel[facei] = pomega_Y_fuel[facei]*prho[facei];  
            

            pCp[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[2]);

            pHf[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[3]);

            pWt[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[4]);

            pmu[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[5])*prho[facei];

            if (baseFSD<ReactionThermo>::flameletT_)
            {
                pT[facei] = this->lookup4d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->NGC,this->gc_Tb3,gc,
                                            this->tableValues_[6]);
            }
            else
            {
                pT[facei] = (pHe[facei]-pHf[facei])/pCp[facei] + this->T0;  
            }

            pT[facei] = max(pT[facei], this->TMin_);
            pT[facei] = min(pT[facei], this->TMax_);


            // // -------------------- Yis begin ------------------------------

            for (int yi=0; yi<this->NY; yi++)
            {
                word specieName2Update = this->speciesNames_table_[yi];
                const label& specieLabel2Update = this->chemistryPtr_->species()[specieName2Update];
                this->Y_[specieLabel2Update].boundaryFieldRef()[patchi][facei]  = this->lookup3d(this->NZ, this->z_Tb3, pZ[facei],
                                            this->NC, this->c_Tb3, pc[facei],
                                            this->NGZ, this->gz_Tb3, gz,
                                            this->tableValues_[this->NS+yi]);  
            }

            // -------------------- Yis end ------------------------------

            if (
                // pgradC_mod[facei] > 0.05 / delta_max
                pgradC_mod[facei] <= pfsd[facei] 
                // && pomega_c[facei] > 1.0
                // && pc[facei] > this->small
                && pfsd[facei] > 0.05 / delta_max
                && pZ[facei] >= 0.99*Zl && pZ[facei] <= 1.01*Zr
                && this->combustion_ && pc[facei] > this->small)
            {
                pn_FSD[facei] = -pgrad_C[facei] / pfsd[facei];
            }
            else if (pgradC_mod[facei] > 0.05 / delta_max
                // && pomega_c[facei] > 1.0
                // && pfsd[facei] > (0.05 - this->small) / delta_max
                && pZ[facei] >= 0.99*Zl && pZ[facei] <= 1.01*Zr
                && this->combustion_ && pc[facei] > this->small) 
            {
                pn_FSD[facei] = -pgrad_C[facei] / pgradC_mod[facei];
            }
            else
            {
                pn_FSD[facei] = vector::zero;
            }
        }
    }

    volScalarField n_FSD_div_ = fvc::div(this->n_FSD_);

    tmp<volTensorField> t_SRE_fac_ = - this->n_FSD_ * this->n_FSD_ + 2.0 * this->IField_ / 3.0 
                + this->IField_ * ( this->n_FSD_ & this->n_FSD_)  / 3.0;

    //- source term in FSD equation, representing the strain-rate effects 
    this->source_fsd_SRE1_ = this->fsd_ * ( t_SRE_fac_ && fvc::grad(this->U_) ) ;
    this->source_fsd_SRE2_ = this->fsd_ * this->phi_fac_FSD_ * this->gamma_k_ * Foam::sqrt(k)
                            / (5.0*(this->delta_ + this->small * this->IField_m1_)) ;

    this->source_fsd_SRE1_.max(0.0);

    // orientation factor
    volScalarField orientation_fac_ = 1.0 - ( this->n_FSD_ & this->n_FSD_ );

    orientation_fac_.max(0.0);

    //-  tensor modelling the strain rate
    volTensorField NTensorField_ = this->n_FSD_ * this->n_FSD_ 
                    + (1.0 - (this->n_FSD_ & this->n_FSD_)) / 3.0 * this->IField_;

    this->Ka_temp_ = fvc::div(this->U_) - ( NTensorField_ && fvc::grad(this->U_) )
                        + this->gamma_k_ * Foam::sqrt(k) / (50.0*(this->delta_ + this->small * this->IField_m1_)) ;

    this->Ka_temp_.max(0.0);

    this->Ka_temp1_ = fvc::div(this->U_);
    this->Ka_temp2_ =  - ( NTensorField_ && fvc::grad(this->U_) );
    this->Ka_temp3_ = this->gamma_k_ * Foam::sqrt(k) / (50.0*(this->delta_ + this->small * this->IField_m1_)) ;

    scalar maxSL0_ = max(this->SL0_FSDCells_);

    // this->Ka_temp_ = this->Ka_temp_ + this->gamma_k_ * Foam::sqrt(k) / (5.0*(this->delta_)) ;

    // this->source_fsd_SRE_ = this->fsd_ * this->Ka_temp_;

    // volScalarField Sd_temp = this->omega_c_ + ( this->n_FSD_ & fvc::grad(mu/this->Sc_ * this->n_FSD_ & gradC_) );

    volScalarField n_FSD_mod = Foam::mag(this->n_FSD_);

    forAll(this->rho_, celli)  
    {
        if (n_FSD_mod.primitiveFieldRef()[celli] < this->small)
        {
            this->source_fsd_SRE1_.primitiveFieldRef()[celli] = 0.0;
        }

        Zl = this->z_Tb5[0];
        Zr = this->z_Tb5[this->NZL-1];

        // look up ZK-K-I0 table
        ZKl = this->ZK_Tb3[0];
        ZKr = this->ZK_Tb3[this->NZK-1];

        if ( this->cCells_[celli] < 0.99  && this->cCells_[celli] > this->small)  // finish calculate source_fsd_FCE2_
        {
            this->source_fsd_FCE2_[celli] = orientation_fac_[celli] * this->fsdCells_[celli] * this->fsdCells_[celli] 
                        / (1.0 - this->cCells_[celli]) 
                        * this->beta_FSD_ * this->SL0_FSDCells_[celli];
        }
        else if (this->cCells_[celli] > this->small)
        {
            this->source_fsd_FCE2_[celli] = orientation_fac_.primitiveFieldRef()[celli] * this->fsdCells_[celli] * this->fsdCells_[celli]
                        / (0.01)
                        * this->beta_FSD_ * this->SL0_FSDCells_[celli];
        }
        else
        {
            this->source_fsd_FCE2_[celli] = 0.0;
        }

        if ( (this->SL0_FSDCells_[celli] > this->small) 
            // and (maxSL0_ > this->small)
            and this->ZCells_[celli] >= ZKl && this->ZCells_[celli] <= ZKr  )
        {
            this->KaCells_[celli] = this->Ka_temp_.primitiveFieldRef()[celli] 
                // * this->thermoThickness_L0_FSDCells_[celli] /maxSL0_;
                * this->thermoThickness_L0_FSDCells_[celli] / this->SL0_FSDCells_[celli];
        }
        else
        {
            this->KaCells_[celli] = 0.0;
        }
        

        if (this->ZCells_[celli] >= ZKl && this->ZCells_[celli] <= ZKr
            && this->combustion_ && this->cCells_[celli] > this->small)  
        {
            this->I_s_FSDCells_[celli] = 1.0;
            // this->I_s_FSDCells_[celli] = this->lookup2d(this->NZK, this->ZK_Tb3, this->ZCells_[celli],
            //                                     this->NK, this->K_Tb3, this->KaCells_[celli],
            //                                     this->I0_Tab);

            // if (this->I_s_FSDCells_[celli] < 1.0)  
            // {
            //     this->I_s_FSDCells_[celli] = 1.0;
            // }                           
        }
        else
        {
            this->I_s_FSDCells_[celli] = 0.0;  
        }

        this->SdACells_[celli] = this->rho_uCells_[celli] * this->SL0_FSDCells_[celli]
                * this->I_s_FSDCells_[celli] / this->rho_[celli]
                - ( this->curv_SdA_ ?
                    ( muCells[celli] / (this->rho_[celli] * this->Sc_) * n_FSD_div_[celli] )
                    : 0.0 );

        // if (gradC_modCells[celli] > 0.05 / delta_max
        //     // && this->omega_cCells_[celli] > 1.0
        //     )
        // {
        //     this->SdACells_[celli] = Sd_temp[celli] / (this->rho_[celli] * gradC_modCells[celli]) 
        //             - muCells[celli] / (this->rho_[celli] * this->Sc_) * n_FSD_div_[celli];
        // }
        // // else if (n_FSD_mod.primitiveFieldRef()[celli] > this->small)
        // // {
        // //     this->SdACells_[celli] = - muCells[celli] / (this->rho_[celli] * this->Sc_) * n_FSD_div_[celli];
        // // }
        // else
        // {
        //     this->SdACells_[celli] = 0.0;
        // }

    }   

    volScalarField::Boundary& source_fsd_FCE2BF = this->source_fsd_FCE2_.boundaryFieldRef();
    volScalarField::Boundary& source_fsd_SRE1BF = this->source_fsd_SRE1_.boundaryFieldRef();
    volScalarField::Boundary& Ka_tempBF = this->Ka_temp_.boundaryFieldRef();
    volScalarField::Boundary& n_FSD_divBF = n_FSD_div_.boundaryFieldRef();
    volScalarField::Boundary& SdABF = this->SdA_.boundaryFieldRef();
    volScalarField::Boundary& orientation_facBF = orientation_fac_.boundaryFieldRef();

    // volScalarField::Boundary& Sd_tempBF = Sd_temp.boundaryFieldRef();
    volScalarField::Boundary& n_FSD_modBF = n_FSD_mod.boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {   
        fvPatchScalarField& prho = rhoBF[patchi];
        fvPatchScalarField& pZ = ZBF[patchi];
        fvPatchScalarField& pc = cBF[patchi];
        fvPatchScalarField& pKa = KaBF[patchi];
        fvPatchScalarField& pKa_temp = Ka_tempBF[patchi];
        fvPatchScalarField& psource_fsd_FCE2 = source_fsd_FCE2BF[patchi];
        fvPatchScalarField& pI_s_FSD = I_s_FSDBF[patchi];
        fvPatchScalarField& pSL0_FSD = SL0_FSDBF[patchi];
        fvPatchScalarField& pthermoThickness_L0_FSD = thermoThickness_L0_FSDBF[patchi];
        fvPatchScalarField& pn_FSD_div = n_FSD_divBF[patchi];
        fvPatchScalarField& pSdA = SdABF[patchi];
        fvPatchScalarField& prho_u = rho_uBF[patchi];
        fvPatchScalarField& pmu = muBF[patchi];
        fvPatchScalarField& pfsd = fsdBF[patchi];
        fvPatchScalarField& psource_fsd_SRE1 = source_fsd_SRE1BF[patchi];
        fvPatchScalarField& porientation_fac = orientation_facBF[patchi];

        // fvPatchScalarField& pSd_temp = Sd_tempBF[patchi];
        fvPatchVectorField& pgrad_C = grad_CBF[patchi];
        fvPatchScalarField& pgradC_mod = gradC_modBF[patchi];

        fvPatchScalarField& pn_FSD_mod = n_FSD_modBF[patchi];
        fvPatchScalarField& pomega_c = omega_cBF[patchi];

        forAll(pZ, facei)
        {
            if (pn_FSD_mod[facei] < this->small)
            {
                psource_fsd_SRE1[facei] = 0.0;
            }

            Zl = this->z_Tb5[0];
            Zr = this->z_Tb5[this->NZL-1];

            ZKl = this->ZK_Tb3[0];
            ZKr = this->ZK_Tb3[this->NZK-1];
            
            if (pc[facei] < 0.99 && pc[facei] > this->small)
            {
                psource_fsd_FCE2[facei] = porientation_fac[facei] * pfsd[facei] * pfsd[facei] / (1.0 - pc[facei])
                        * this->beta_FSD_ * pSL0_FSD[facei];
            }
            else if (pc[facei] > this->small)
            {
                psource_fsd_FCE2[facei] = porientation_fac[facei] * pfsd[facei] * pfsd[facei] / (0.01)
                        * this->beta_FSD_ * pSL0_FSD[facei];
            }
            else
            {
                psource_fsd_FCE2[facei] = 0.0;
            }

            if ( ( pSL0_FSD[facei] > this->small ) 
                // and (maxSL0_ > this->small) 
                and pZ[facei] >= ZKl && pZ[facei] <= ZKr)
            {
                pKa[facei] = pKa_temp[facei] * pthermoThickness_L0_FSD[facei] 
                            // / maxSL0_;
                            / pSL0_FSD[facei];
            }
            else
            {
                pKa[facei] = 0.0;
            }

            if(pZ[facei] >= ZKl && pZ[facei] <= Zr
                && this->combustion_ && pc[facei] > this->small) 
            {
                pI_s_FSD[facei] = 1.0;
                // pI_s_FSD[facei] = this->lookup2d(this->NZK, this->ZK_Tb3, pZ[facei],
                //                                 this->NK, this->K_Tb3, pKa[facei],
                //                                 this->I0_Tab);

                // if (pI_s_FSD[facei] < 1.0)
                // {
                //     pI_s_FSD[facei] = 1.0;
                // }                                
            }
            else
            {
                pI_s_FSD[facei] = 0.0; 
            }

            pSdA[facei] = prho_u[facei] * pSL0_FSD[facei]
                            * pI_s_FSD[facei] / prho[facei]
                            - ( this->curv_SdA_ ?
                                ( pmu[facei] / (prho[facei] * this->Sc_) * pn_FSD_div[facei] )
                                : 0.0 );



            // label cellID = this->mesh().boundary()[patchi].faceCells()[facei]; 

            // if (
            //     pgradC_mod[facei] > 0.05 / delta_max
            //     // && pomega_c[facei] > 1.0
            //     )
            // {
            //     pSdA[facei] = pSd_temp[facei] / (prho[facei]*pgradC_mod[facei])
            //                 - pmu[facei] / (prho[facei] * this->Sc_) * pn_FSD_div[facei];
            // }
            // // else if (pn_FSD_mod[facei] > this->small)
            // // {
            // //     pSdA[facei] = - pmu[facei] / (prho[facei] * this->Sc_) * pn_FSD_div[facei];
            // // }
            // else
            // {
            //     pSdA[facei] = 0.0;
            // }
        }
    }

    // this->SdA_.max(0.0);


    this->source_fsd_flaProp_ = fvc::laplacian(this->SdA_, this->c_);

    // this->source_fsd_flaProp_.max(0.0);  // ？？

    forAll(this->rho_, celli)  
    {
        if ( ( (this->cCells_[celli] > 0.5 && this->fsdCells_[celli] < 1.0e-4) 
            or this->cCells_[celli] < 1.0e-4))
        {
            this->source_fsd_flaProp_.primitiveFieldRef()[celli] = 0.0;
        }
    }
    
    volScalarField::Boundary& source_fsd_flaPropBF = this->source_fsd_flaProp_.boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {
        fvPatchScalarField& psource_fsd_flaProp = source_fsd_flaPropBF[patchi];
        fvPatchScalarField& pc = cBF[patchi];
        fvPatchScalarField& pfsd = fsdBF[patchi];

        forAll(pc, facei)
        {
            if ( (pc[facei] > 0.5 && pfsd[facei] < 1.0e-4) or pc[facei] < 1.0e-4)
            {
                psource_fsd_flaProp[facei] = 0.0;
            }
        }
    }


    this->source_fsd_FCE1_ = this->SdA_ * this->fsd_ * fvc::div(this->n_FSD_);
    


    // update rho
    double R_uniGas_ = 8.314e3;
    double p_operateDim_ = this->coeffs().lookupOrDefault("p_operateDim", this->incompPref_);

    forAll(this->rho_.boundaryFieldRef(), patchi)   
    {
        fvPatchScalarField& ppsi_ = this->psi_.boundaryFieldRef()[patchi];
        fvPatchScalarField& prho_ = this->rho_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pWt = this->Wt_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryFieldRef()[patchi];
        fvPatchScalarField pp_ = this->p_.boundaryField()[patchi];

        forAll(prho_, facei)   
        {
            ppsi_[facei] = pWt[facei] / (R_uniGas_*pT[facei]);

            if(this->incompPref_ > 0.0) 
            {  
                prho_[facei] = p_operateDim_*ppsi_[facei]; 
            }
            else 
            {
                prho_[facei] = pp_[facei]*ppsi_[facei];
            }
        }
    }

    dimensionedScalar R_uniGas("R_uniGas",dimensionSet(1,2,-2,-1,-1,0,0),8.314e3);
    this->psi_ = this->Wt_/(R_uniGas*this->T_);
  
    if(this->mesh().time().timeIndex() > 0)   
    {
        dimensionedScalar p_operateDim("p_operateDim", dimensionSet(1,-1,-2,0,0,0,0),this->incompPref_);  

        if(this->incompPref_ > 0.0) 
        {  
            this->rho_ = p_operateDim*this->psi_; 
        }
        else 
        {
            this->rho_ = this->p_*this->psi_;
        }
    }
    this->rho_inRhoThermo_ = this->rho_;
}