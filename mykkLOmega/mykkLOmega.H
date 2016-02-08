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

    Class
    Foam::incompressible::RASModels::mykkLOmega

    Description
    Low Reynolds-number k-kl-omega turbulence model for
    incompressible flows.

    Turbulence model described in:
    \verbatim
    D. Keith Walters, Davor Cokljat
    "A Three-Equation Eddy-Viscosity Model for Reynold-Averaged
    Navier-Stokes Simulations of Transitional Flow"
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
    mykkLOmegaCoeffs
    {
    A0             4.04
    As             2.12
    Av             6.75
    Abp            0.6
    Anat           200
    Ats            200
    CbpCrit        1.2
    Cnc            0.1
    CnatCrit       1250
    Cint           0.75
    CtsCrit        1000
    CrNat          0.02
    C11            3.4e-6
    C12            1.0e-10
    CR             0.12
    CalphaTheta    0.035
    Css            1.5
    CtauL          4360
    Cw1            0.44
    Cw2            0.92
    Cw3            0.3
    CwR            1.5
    Clambda        2.495
    CmuStd         0.09
    Prtheta        0.85
    Sigmak         1
    Sigmaw         1.17
    }
    \endverbatim

    SourceFiles
    mykkLOmega.C

    \*---------------------------------------------------------------------------*/

#ifndef mykkLOmega_H
#define mykkLOmega_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------* \
    Class mykkLOmega Declaration
    \*---------------------------------------------------------------------------*/

class mykkLOmega
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private memmber functions
    
    tmp<volScalarField> fv(const volScalarField& Ret) const;
    
    tmp<volScalarField> fINT() const;
    
    tmp<volScalarField> fSS(const volScalarField& Omega) const;
    
    tmp<volScalarField> Cmu(const volScalarField& S) const;
    
    tmp<volScalarField> BetaTS(const volScalarField& Rew, const volScalarField& L) const;

    tmp<volScalarField> L(const volScalarField& Rew) const;

    tmp<volScalarField> CTSCrit(const volScalarField& L) const;

    tmp<volScalarField> CnatCrit(const volScalarField& L) const;

    tmp<volScalarField> fTaul
    (
        const volScalarField& lambdaEff,
        const volScalarField& ktL,
        const volScalarField& Omega
    ) const;
    
    tmp<volScalarField> alphaT
    (
        const volScalarField& lambdaEff,
        const volScalarField& fv,
        const volScalarField& ktS
    ) const;
    
    tmp<volScalarField> fOmega
    (
        const volScalarField& lambdaEff,
        const volScalarField& lambdaT
    ) const;
    
    tmp<volScalarField> gammaBP(const volScalarField& Omega) const;
    
    tmp<volScalarField> gammaNAT
    (
        const volScalarField& ReOmega,
        const volScalarField& fNatCrit,
        const volScalarField& L
    ) const;
    
    void boundMinMax
    (
        volScalarField& vsf,
        const dimensionedScalar& vsf0,
        const dimensionedScalar& vsf1
    ) const;

    tmp<volScalarField> D(const volScalarField& k) const;

protected:
    
    // Protected data
    
    // Model coefficients
    
    dimensionedScalar A0_;
    dimensionedScalar As_;
    dimensionedScalar Av_;
    dimensionedScalar Abp_;
    dimensionedScalar Anat_;
    dimensionedScalar Ats_;
    dimensionedScalar CbpCrit_;
    dimensionedScalar Cnc_;
    dimensionedScalar CnatCrit_;
    dimensionedScalar CnatApg_;
    dimensionedScalar Cint_;
    dimensionedScalar CtsCrit_;
    dimensionedScalar CtsApg_;
    dimensionedScalar CrNat_;
    dimensionedScalar C11_;
    dimensionedScalar C12_;
    dimensionedScalar CR_;
    dimensionedScalar CalphaTheta_;
    dimensionedScalar Css_;
    dimensionedScalar CtauL_;
    dimensionedScalar Cw1_;
    dimensionedScalar Cw2_;
    dimensionedScalar Cw3_;
    dimensionedScalar CwR_;
    dimensionedScalar Clambda_;
    dimensionedScalar CmuStd_;
    dimensionedScalar Prtheta_;
    dimensionedScalar Sigmak_;
    dimensionedScalar Sigmaw_;
    
    dimensionedScalar kMax_;
    dimensionedScalar omegaMax_;

    // Fields
    
    volScalarField kt_;
    volScalarField kl_;
    volScalarField omega_;
    volScalarField epsilon_;
    const volScalarField& y_;
    
    Switch medinaPhiBP_;
    Switch furstPG_;
    dimensionedScalar CDT_;
    dimensionedScalar viscosityRatioLimit_;

    // Temporary fields for model analysis
    volScalarField fv_;
    volScalarField fINT_;
    volScalarField fSS_;
    volScalarField Cmu_;
    volScalarField BetaTS_;
    volScalarField BetaNAT_;
    volScalarField BetaBP_;
    volScalarField fTaul_;
    volScalarField fOmega_;
    volScalarField gammaBP_;
    volScalarField gammaNAT_;
    volScalarField ReOmega_;
    volScalarField RBP_;
    volScalarField RNAT_;
    volScalarField Pkt_;
    volScalarField Pkl_;
    volScalarField L_;

    virtual void correctNut();

public:
    
    //- Runtime type information
    TypeName("mykkLOmega");

    // Constructors

    //- Construct from components
    mykkLOmega
    (
        const geometricOneField& alpha,
        const geometricOneField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );


    //- Destructor
    virtual ~mykkLOmega()
        {}


    // Member Functions

    //- Return the turbulence viscosity
    virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DkEff", alphaT/Sigmak_ + nu())
                );
        }

    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("DomegaEff", alphaT/Sigmaw_ + nu())
                );
        }

    //- Return the laminar kinetic energy
    virtual tmp<volScalarField> kl() const
        {
            return kl_;
        }

    //- Return the turbulence kinetic energy
    virtual tmp<volScalarField> kt() const
        {
            return kt_;
        }

    //- Return the turbulence specific dissipation rate
    virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }
    
    virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject
                    (
                        "k",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                        kt_ + kl_,
                        omega_.boundaryField().types()
                    )
                );
        }
    
    //- Return the turbulence kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

    //- Read RASProperties dictionary
    virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
