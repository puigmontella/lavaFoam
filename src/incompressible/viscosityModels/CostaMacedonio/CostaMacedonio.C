/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "CostaMacedonio.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CostaMacedonio, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        CostaMacedonio,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::CostaMacedonio::calcNu() const
{
const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T");

// volScalarField exp_term = exp(T * b_ );
    // volScalarField viscosity = nu_ref_ * T * b_ ;

    // Debug prints
    // std::cout << "T: " << "\n";
//     forAll ( T, celli) //loop through cell centres
// {
//   std::cout << "T: " << T[celli] << "\n";
//   // std::cout << "exp(T * b + 1.0): " << exp_term[celli] << "\n";
//   std::cout << "viscosity before clamp: " << viscosity[celli] << "\n";
// }

    // std::cout << "b: " << b_ << "\n";


    // return max
    // (
    //     nuMin_,
    //     min
    //     (
    //         nuMax_,
    //         nu_ref_*T*b_
    //     )
    // );
    //


// return nu_ref_;
// return nu_ref_*exp(b_ * (T - T_ref_));
// return max(min(nu_ref_*exp(b_ * (T - T_ref_)),nuMax_),nuMin_);
// return max(nuMin_,min(nuMax_,(b*nu_ref_*(T-T_ref_))));
// return max(nuMin_,min(nuMax_,(k_*(T-Tbase_))*pow(max(dimensionedScalar("one", dimTime, 1.0)*strainRate(),dimensionedScalar("VSMALL", dimless, VSMALL)),n_.value() - scalar(1.0))));
return max(nuMin_,min(nuMax_,(nu_ref_-nu_ref_*b_*(T-T_ref_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::CostaMacedonio::CostaMacedonio
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CostaMacedonioCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),



    nu_ref_("nu_ref",  dimensionSet( 0, 2, -1, 0, 0, 0, 0), CostaMacedonioCoeffs_),
    b_("b",  dimensionSet( 0, 0, 0, -1, 0, 0, 0), CostaMacedonioCoeffs_),
    // kslope_("kslope",  dimensionSet( 0, 2, -1, -1, 0, 0, 0), CostaMacedonioCoeffs_),
    T_ref_("T_ref",  dimensionSet( 0, 0, 0, 1, 0, 0, 0), CostaMacedonioCoeffs_),
    nuMin_("nuMin",  dimensionSet( 0, 2,-1, 0, 0, 0, 0), CostaMacedonioCoeffs_),
    nuMax_("nuMax",  dimensionSet( 0, 2,-1, 0, 0, 0, 0), CostaMacedonioCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::CostaMacedonio::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CostaMacedonioCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    CostaMacedonioCoeffs_.readEntry("nu_ref", nu_ref_);
    CostaMacedonioCoeffs_.readEntry("b", b_);
    // CostaMacedonioCoeffs_.readEntry("kslope", kslope_);
    CostaMacedonioCoeffs_.readEntry("T_ref",T_ref_);
    CostaMacedonioCoeffs_.readEntry("nuMin",nuMin_);
    CostaMacedonioCoeffs_.readEntry("nuMax",nuMax_);
    // CostaMacedonioCoeffs_.readEntry("Tbase",Tbase_);
    return true;
}


// ************************************************************************* //
