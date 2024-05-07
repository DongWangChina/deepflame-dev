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
    Foam::combustionModels::FSD

Description
    Partially stirred reactor turbulent combustion model.

    This model calculates a finite rate, based on both turbulence and chemistry
    time scales.  Depending on mesh resolution, the Cmix parameter can be used
    to scale the turbulence mixing time scale.

SourceFiles
    tableLookup_FSD.C

\*---------------------------------------------------------------------------*/

#ifndef FSD_H
#define FSD_H

#include "baseFSD.H"
#include "tableSolver_FSD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class FSD Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FSD
:
    public baseFSD<ReactionThermo>,
    public tableSolver_FSD
{

public:

    //- Runtime type information
    TypeName("FSD");


    // Constructors

        //- Construct from components
        FSD
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        FSD(const FSD&);


    //- Destructor
    virtual ~FSD();


    // Member Operators

        //- Correct combustion rate
        virtual void correct();

        //- retrieval data from table
        virtual void retrieval();



        //- Disallow default bitwise assignment
        void operator=(const FSD&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FSD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
