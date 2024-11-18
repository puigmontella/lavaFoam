/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::addTImmiscibleIncompressibleTwoPhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    addTImmiscibleIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef addTImmiscibleIncompressibleTwoPhaseMixture_H
#define addTImmiscibleIncompressibleTwoPhaseMixture_H

#include "addTIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class addTImmiscibleIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class addTImmiscibleIncompressibleTwoPhaseMixture
:
    public addTIncompressibleTwoPhaseMixture,
    public interfaceProperties
{

public:

    const dimensionedScalar emissivity_;
    const dimensionedScalar sigma_SB_;
    const dimensionedScalar heatTransferCoeff_;
    const dimensionedScalar fractionalAreaExposed_;
    const dimensionedScalar T_env_;


    //- Construct from components
    addTImmiscibleIncompressibleTwoPhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );


    //- Destructor
    virtual ~addTImmiscibleIncompressibleTwoPhaseMixture() = default;


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            addTIncompressibleTwoPhaseMixture::correct();
            interfaceProperties::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read();


        // Add Rad
        // Calculates the radiative source term
        tmp<fvScalarMatrix> calcSourceRadiation(
          const volVectorField &U,
          volScalarField &T,
          volScalarField &RadiativeCoeff
        );
        // - Calculates convective source term , returned as volume fraction rate
        tmp<fvScalarMatrix>calcSourceForcedConvection(
          const volVectorField &U,
          volScalarField &T
        );

        // Calculates the difference between the fourth powers
        // of temperature and environmental temperature
        void calcT4mTenv4(
          volScalarField &T,
          volScalarField &T4mTenv4
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
