/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::prototype_retract_motion_gen

Description
    Foam::prototype_retract_motion_gen

SourceFiles
    prototype_retract_motion_gen.C
	
calculate rev_pitch angle,  do not use specified rev_pitch
meant for 6dof

\*---------------------------------------------------------------------------*/

#ifndef prototype_retract_motion_gen_H
#define prototype_retract_motion_gen_H

#include "fixedValuePointPatchField.H"
//#include "rotatingMotion_2dof_pitch_y_translate_newmark.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class prototype_retract_motion_gen Declaration
\*---------------------------------------------------------------------------*/

class prototype_retract_motion_gen
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        vector axis_;
        scalar angle0_;
		scalar x1;
		scalar z1;
		scalar x1_outer;
		scalar z1_outer;
		scalar x2;
		scalar z2;
		scalar x2_outer;
		scalar z2_outer;
        scalar amplitude_;
		scalar rev_pitch_;
		scalar max_retract_angle_;
        scalar omega_;

        pointField p0_;
		
		//mutable rotatingMotion_2dof_pitch_y_translate_newmark angle_old;
		mutable bool firstRun_;
		mutable scalar saved_time_;


public:

    //- Runtime type information
    TypeName("prototype_retract_motion_gen");


    // Constructors

        //- Construct from patch and internal field
        prototype_retract_motion_gen
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prototype_retract_motion_gen
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        prototype_retract_motion_gen
        (
            const prototype_retract_motion_gen&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new prototype_retract_motion_gen
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        prototype_retract_motion_gen
        (
            const prototype_retract_motion_gen&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new prototype_retract_motion_gen
                (
                    *this,
                    iF
                )
            );
        }
		/*
		IOdictionary orientationdict
		(
			IOobject
			(
				"orientationdict",
				runTime.constant(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);
		*/

    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
