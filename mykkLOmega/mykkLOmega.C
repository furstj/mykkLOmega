/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "mykkLOmega.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mykkLOmega, 0);
addToRunTimeSelectionTable(RASModel, mykkLOmega, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> mykkLOmega::fv(const volScalarField& Ret) const
{
    return(1.0 - exp(-sqrt(Ret)/Av_));
}


tmp<volScalarField> mykkLOmega::fINT() const
{
    return
    (
        min
        (
            kt_/(Cint_*(kl_ + kt_)),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );
}


tmp<volScalarField> mykkLOmega::fSS(const volScalarField& Omega) const
{
    return(exp(-sqr(Css_*nu()*Omega/kt_)));
}


tmp<volScalarField> mykkLOmega::Cmu(const volScalarField& S) const
{
    return(1.0/(A0_ + As_*(S/omega_)));
}


tmp<volScalarField> mykkLOmega::BetaTS(const volScalarField& Rew) const
{
    return(scalar(1) - exp(-sqr(max(Rew - CtsCrit_, scalar(0)))/Ats_));
}


tmp<volScalarField> mykkLOmega::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -CtauL_*ktL
          /
            (
                sqr
                (
		 max(lambdaEff*Omega,
		     dimensionedScalar
		     (
		      "ROOTVSMALL",
		      dimLength*inv(dimTime),
		      ROOTVSMALL
		      )
		     )
		 )
	     )
	 )
     );
}
  

tmp<volScalarField> mykkLOmega::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}


tmp<volScalarField> mykkLOmega::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -0.41
           * pow4
            (
                lambdaEff
		/ max(lambdaT, 
		      dimensionedScalar
		      (
		       "ROTVSMALL",
		       lambdaT.dimensions(),
		       ROOTVSMALL
		       )
		      )
	     )
	 )
     );
}
  

tmp<volScalarField> mykkLOmega::gammaBP(const volScalarField& Omega) const
{
    return
    (
        max
        (
            kt_/nu()
          / max(
                Omega,
		dimensionedScalar("ROOTVSMALL", Omega.dimensions(), ROOTVSMALL)
		)
	    - CbpCrit_,
            scalar(0)
	 )
     );
}
  

tmp<volScalarField> mykkLOmega::gammaNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
	    / max(
		  fNatCrit, dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mykkLOmega::mykkLOmega
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            coeffDict_,
            2.12
        )
    ),
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            coeffDict_,
            6.75
        )
    ),
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
            coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            coeffDict_,
            0.035
        )
    ),
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
            coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.92
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            coeffDict_,
            2.495
        )
    ),
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
            coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
            coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            coeffDict_,
            1.17
        )
    ),
    kt_
    (
        IOobject
        (
            "kt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("kt", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    kl_
    (
        IOobject
        (
            "kl",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("kl", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    y_(mesh_)
{
  kMin_ = dimensionedScalar("kMin", sqr(dimVelocity), VSMALL);
  omegaMin_ = dimensionedScalar("omegaMin", inv(dimTime), VSMALL);

    bound(kt_, kMin_);
    bound(kl_, kMin_);
    bound(omega_, omegaMin_);

    nut_ = kt_/(omega_ + omegaMin_);
    nut_.correctBoundaryConditions();

    Info << ROOTVSMALL << endl;
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> mykkLOmega::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*(kt_) - nut_*twoSymm(fvc::grad(U_)),
            kt_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> mykkLOmega::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> mykkLOmega::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


bool mykkLOmega::read()
{
    if (RASModel::read())
    {
        A0_.readIfPresent(coeffDict());
        As_.readIfPresent(coeffDict());
        Av_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Anat_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Ats_.readIfPresent(coeffDict());
        CbpCrit_.readIfPresent(coeffDict());
        Cnc_.readIfPresent(coeffDict());
        CnatCrit_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CtsCrit_.readIfPresent(coeffDict());
        CrNat_.readIfPresent(coeffDict());
        C11_.readIfPresent(coeffDict());
        C12_.readIfPresent(coeffDict());
        CR_.readIfPresent(coeffDict());
        CalphaTheta_.readIfPresent(coeffDict());
        Css_.readIfPresent(coeffDict());
        CtauL_.readIfPresent(coeffDict());
        Cw1_.readIfPresent(coeffDict());
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        CwR_.readIfPresent(coeffDict());
        Clambda_.readIfPresent(coeffDict());
        CmuStd_.readIfPresent(coeffDict());
        Prtheta_.readIfPresent(coeffDict());
        Sigmak_.readIfPresent(coeffDict());
        Sigmaw_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void mykkLOmega::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        y_.boundaryField() = max(y_.boundaryField(), VSMALL);
    }


    const volScalarField kT(kt_ + kl_);

    const volScalarField lambdaT( sqrt(kt_)/max(omega_,omegaMin_) );

    const volScalarField lambdaEff(min(Clambda_*y_, lambdaT));

    const volScalarField fw
    (
        pow(
            lambdaEff/max(lambdaT,dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2./3.)
    );

    const volTensorField gradU(fvc::grad(U_));

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2.0*magSqr(symm(gradU)));

    const volScalarField ktS(fSS(Omega)*fw*kt_);

    const volScalarField nuts
    (
     fv(sqr(fw)*kt_/nu()/omega_)
     *fINT()
       *Cmu(sqrt(S2))*sqrt(ktS)*lambdaEff
    );
    const volScalarField Pkt(nuts*S2);

    const volScalarField ktL(kt_ - ktS);
    const volScalarField ReOmega(sqr(y_)*Omega/nu());
    const volScalarField nutl
    (
        min
        (
	 C11_*fTaul(lambdaEff, ktL, Omega)*Omega*sqr(lambdaEff)
          * sqrt(ktL)*lambdaEff/nu()
          + C12_*BetaTS(ReOmega)*ReOmega*sqr(y_)*Omega
        ,
            0.5*(kl_ + ktL)/sqrt(S2)
        )
    );

    const volScalarField Pkl(nutl*S2);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv(sqr(fw)*kt_/nu()/max(omega_, omegaMin_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);

    const volScalarField Rbp
    (
        CR_*(1.0 - exp(-gammaBP(Omega)()/Abp_))*omega_
	/ max(fw,fwMin)
    );

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y_/nu()));
    // Natural source term divided by kl_
    const volScalarField Rnat
    (
        CrNat_*(1.0 - exp(-gammaNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );


    // Turbulence specific dissipation rate equation
    omega_.boundaryField().updateCoeffs();
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian
        (
            DomegaEff(alphaTEff),
            omega_,
            "laplacian(alphaTEff,omega)"
        )
     ==
        Cw1_*Pkt*omega_/(kt_ + kMin_)
      + fvm::SuSp
        (
            (CwR_/(fw + fwMin) - 1.0)*kl_*(Rbp + Rnat)/(kt_ + kMin_)
          , omega_
        )
        - fvm::Sp(Cw2_*omega_*sqr(fw), omega_)
      + Cw3_*fOmega(lambdaEff, lambdaT)*alphaTEff*sqr(fw)*sqrt(kt_)/pow3(y_)
    );


    omegaEqn().relax();
    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    // Laminar kinetic energy equation
    const volScalarField Dl(nu()*magSqr(fvc::grad(sqrt(kl_))));
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(kl_)
      + fvm::div(phi_, kl_)
      - fvm::Sp(fvc::div(phi_), kl_)
      - fvm::laplacian(nu(), kl_, "laplacian(nu,kl)")
     ==
        Pkl
        - fvm::Sp(Rbp, kl_)
        - fvm::Sp(Rnat, kl_)
        - fvm::Sp(Dl/max(kl_,kMin_), kl_)
    );

    klEqn().relax();
    klEqn().boundaryManipulate(kl_.boundaryField());

    solve(klEqn);
    bound(kl_, kMin_);

    // Turbulent kinetic energy equation
    const volScalarField Dt(nu()*magSqr(fvc::grad(sqrt(kt_))));
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(kt_)
      + fvm::div(phi_, kt_)
      - fvm::Sp(fvc::div(phi_), kt_)
      - fvm::laplacian(DkEff(alphaTEff), kt_, "laplacian(alphaTEff,kt)")
     ==
        Pkt
        + (Rbp + Rnat)*kl_
        - fvm::Sp(Dt/max(kt_,kMin_), kt_)
        - fvm::Sp(omega_, kt_)
    );

    ktEqn().relax();
    ktEqn().boundaryManipulate(kt_.boundaryField());

    solve(ktEqn);
    bound(kt_, kMin_);



    // Re-calculate viscosity
    nut_ = nuts + nutl;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
