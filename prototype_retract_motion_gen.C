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

\*---------------------------------------------------------------------------*/

/*

skip calculation if angle > max_retract_angle
calculate rev_pitch angle,  do not use specified rev_pitch
meant for 6dof

*/

#include "prototype_retract_motion_gen.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

inline bool exists_check(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

using namespace std;

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

prototype_retract_motion_gen::
prototype_retract_motion_gen
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    axis_(Zero),
    angle0_(0.0),
    x1(0.0),
	z1(0.0),
	x1_outer(0.0),
	z1_outer(0.0),
	x2(0.0),
	z2(0.0),
	x2_outer(0.0),
	z2_outer(0.0),
    amplitude_(0.0),
	rev_pitch_(0.0),
	max_retract_angle_(0.0),
    omega_(0.0),
    p0_(p.localPoints()),
	firstRun_(true),
	saved_time_(-1)
{}


prototype_retract_motion_gen::
prototype_retract_motion_gen
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    axis_(dict.lookup("axis")),
    angle0_(dict.get<scalar>("angle0")),
    x1(dict.get<scalar>("x1")),
    z1(dict.get<scalar>("z1")),
    x1_outer(dict.get<scalar>("x1_outer")),
    z1_outer(dict.get<scalar>("z1_outer")),
    x2(dict.get<scalar>("x2")),
    z2(dict.get<scalar>("z2")),
    x2_outer(dict.get<scalar>("x2_outer")),
    z2_outer(dict.get<scalar>("z2_outer")),
	amplitude_(dict.get<scalar>("amplitude")),
	max_retract_angle_(dict.get<scalar>("maximum_retract_angle")),
    omega_(dict.get<scalar>("omega")),
	firstRun_(dict.getOrDefault<bool>("firstRun", true))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());

		/* char cc;

		Info<< " p.size = " << p.size() << endl;

		cin>>cc; */

    }
    else
    {
        p0_ = p.localPoints();

		/* char cc;

		Info<< " p0_.size = " << p0_.size() << endl;

		cin>>cc; */
    }
}


prototype_retract_motion_gen::
prototype_retract_motion_gen
(
    const prototype_retract_motion_gen& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    angle0_(ptf.angle0_),
    x1(ptf.x1),
    z1(ptf.z1),
    x1_outer(ptf.x1_outer),
    z1_outer(ptf.z1_outer),
    x2(ptf.x2),
    z2(ptf.z2),
    x2_outer(ptf.x2_outer),
    z2_outer(ptf.z2_outer),
	amplitude_(ptf.amplitude_),
	rev_pitch_(ptf.rev_pitch_),
	max_retract_angle_(ptf.max_retract_angle_),
    omega_(ptf.omega_),
    p0_(ptf.p0_, mapper),
	firstRun_(ptf.firstRun_),
	saved_time_(ptf.saved_time_)
{}


prototype_retract_motion_gen::
prototype_retract_motion_gen
(
    const prototype_retract_motion_gen& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    angle0_(ptf.angle0_),
    x1(ptf.x1),
    z1(ptf.z1),
    x1_outer(ptf.x1_outer),
    z1_outer(ptf.z1_outer),
    x2(ptf.x2),
    z2(ptf.z2),
    x2_outer(ptf.x2_outer),
    z2_outer(ptf.z2_outer),
	amplitude_(ptf.amplitude_),
	rev_pitch_(ptf.rev_pitch_),
	max_retract_angle_(ptf.max_retract_angle_),
    omega_(ptf.omega_),
    p0_(ptf.p0_),
	firstRun_(ptf.firstRun_),
	saved_time_(ptf.saved_time_)
{}
/*
prototype_retract_motion_gen::
prototype_retract_motion_gen
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
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
		
		//tensor orientation = orientationdict.lookup("orientation");
{}

*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void prototype_retract_motion_gen::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void prototype_retract_motion_gen::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const prototype_retract_motion_gen& aODptf =
        refCast<const prototype_retract_motion_gen>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void prototype_retract_motion_gen::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
	
	if (firstRun_ == true) {
		
		saved_time_ = -1;	//just any negative value
		firstRun_ = false;
	
		if (t.timeIndex() == 1) {
			
			Info << "Starting from new" << endl;
			
		}
		
		else {
			
			Info << "Starting from resumed" << endl;
			
		}
		
	}
	
		
	//constant rate instead of sine
	//omega not impt, key is amplitude and runtime duration
	
	scalar angle = angle0_ + amplitude_*t.value();
	
	
	
	
	if (angle < max_retract_angle_ && t.value() - saved_time_ > 1e-8)
    
	
	
	{
		
		
		/*separated into 3 regions. inner most region near root fixed, not moving
		mid region do rotation varying from 29deg near trailing edge to 45deg near leading edge (for 45deg rotation). slight change in length. center of rotation perpendicular to diagonal_chord x1x2
		outer region do full body translation according to leading edge x1 outer point + slight small angle rotation due to change in angle between new and old x1x2_outer
		*/
		
		Info << "prototype_retract_motion running, current retract angle = " << angle << endl;
		
		
		
		//need to pitch up rotation into consideration
		//not constant, input during setup
		//scalar x1 = 0.;
		//scalar z1 = 0.214;

		//scalar x2 = 0.2;
		//scalar z2 = 0.28;
		
		//scalar diagonal_chord = 0.217;
		
		//diagonal_chord = sqrt((z2 - z1)*(z2 - z1) + (x2 -x1)*(x2 -x1));
		
		//need to pitch up rotation into consideration
		//scalar x1_outer = -0.045;
		//scalar z1_outer = 0.357;

		//scalar x2_outer = 0.128;
		//scalar z2_outer = 0.52;

		
			
		scalar angle_12_outer = atan((z2_outer - z1_outer)/(x2_outer - x1_outer));
			
		scalar diagonal_outer_chord = sqrt((z2_outer - z1_outer)*(z2_outer - z1_outer) + (x2_outer -x1_outer)*(x2_outer -x1_outer));
		
		//vector pitch_axisHat = {0., 0., 1.};
		//vector pitch_origin({0.04825, 0, 0});
				
		//scalar pitch_amplitude_ = 0.0873;//0.175;//0.7854;
		//scalar pitch_omega_ = 3.14;
		//scalar pitch_angle(pitch_amplitude_*sin(pitch_omega_*t.value()));
		
		//scalar ratio = 0.;
		scalar which_side = 0.;
		scalar which_side_outer = 0.;
		//scalar pd = 0.;
		scalar xp = 0.;
		scalar yp = 0.;
		scalar zp = 0.;
		
		
		
			
		vectorField sd(p0_);
		
		//ratio of pt wrt to left/right/top/bottom
		scalar ratio_pt_wrt_l;
		scalar ratio_pt_wrt_r;
		scalar ratio_pt_wrt_t;
		scalar ratio_pt_wrt_b;
		
		
		
		
		//t=0
		
		
		//t=1
		//scalar x1_outer_S0R45 = -0.142;
		//scalar z1_outer_S0R45 = 0.275;

		//scalar x2_outer_S0R45 = 0.034;
		//scalar z2_outer_S0R45 = 0.443;	
		
		
		
		scalar x1_1_outer;
		scalar z1_1_outer;
		
		scalar x2_2_outer;
		scalar z2_2_outer;
		
		scalar x12;
		scalar z12;
		
		
		scalar x1_outer_cur = (x1_outer - x1)*cos(angle) - (z1_outer - z1)*sin(angle) + x1;
		scalar z1_outer_cur = (x1_outer - x1)*sin(angle) + (z1_outer - z1)*cos(angle) + z1;
		
		scalar x2_outer_cur = diagonal_outer_chord*cos(angle_12_outer) + x1_outer_cur;
		scalar z2_outer_cur = diagonal_outer_chord*sin(angle_12_outer) + z1_outer_cur;
		
		/*
		cout << "angle_12_outer,diagonal_outer_chord,x1_outer_cur,x1_outer_cur,x2_outer_cur,z2_outer_cur\n";
		
		cout << angle_12_outer <<  '\n';
		cout << diagonal_outer_chord <<  '\n';
		
		
		cout << x1_outer_cur << '\n' << z1_outer_cur << '\n';
		
		cout << x2_outer_cur <<  '\n';
		cout << z2_outer_cur <<  '\n';
		
		getchar();
		*/
		
		//pitch rotation center
		//scalar pitch_rotation_center_x = 0.00896;
		//scalar pitch_rotation_center_y = -0.04056;
		//scalar pitch_rotation_center_z = 0;
		
		scalar dd = 0;
		
		 
		
		scalar x12_outer;
		scalar z12_outer;
		
		scalar x_init;
		scalar y_init;
		scalar z_init;
		
		scalar xp_new;
		scalar yp_new;
		scalar zp_new;
		
		scalar xp_new2;
		scalar yp_new2;
		scalar zp_new2;
		
		
		/*	
		cout << "iter0,x0,y0,z0\n";
					
		cout << iter0 << '\n';
		cout << x0 << '\n';
		cout << y0 << '\n';
		cout << z0 << '\n';*/
		
			
		
		
		forAll(p0_,iter)
		{
			
			//do mid region 1st	
			
			//ratio = 0.;
			
			//get each pt coordinate
			
			x_init = p0_[iter].component(vector::X);
			y_init = p0_[iter].component(vector::Y);
			z_init = p0_[iter].component(vector::Z);
			
			//rotate /pitch back 1st - no need to
			
			/*xp = (x_init - pitch_origin[0])*cos(orient_pitch_old_rev) - (y_init - pitch_origin[1])*sin(orient_pitch_old_rev) + pitch_origin[0];
			yp = (x_init - pitch_origin[0])*sin(orient_pitch_old_rev) + (y_init - pitch_origin[1])*cos(orient_pitch_old_rev) + pitch_origin[1];
			zp = z_init;*/
			xp = x_init;
			yp = y_init;
			zp = z_init;
			
			
			//find out which side pt belongs to - moveable or non-moveable region 
			// > 0 = red + green circle region
			
			which_side = (zp - z1)*(x2 - x1) - (xp - x1)*(z2 - z1);
			
			which_side_outer = (zp - z1_outer)*(x2_outer - x1_outer) - (xp - x1_outer)*(z2_outer - z1_outer);
			
			// in green region
			
			if (which_side > 0 && which_side_outer < 0) {
				
				ratio_pt_wrt_r = -(z2_outer - (xp*z1_outer*z1_outer - x2*z1_outer*z1_outer - x1*z2_outer*z2_outer + xp*z2_outer*z2_outer - z1_outer*sqrt(x1*x1*z2_outer*z2_outer - 2*x1*x1*z2_outer*zp + x1*x1*zp*zp - 2*x1*x2*z1_outer*z2_outer + 2*x1*x2*z1_outer*zp + 2*x1*x2*z2_outer*zp - 2*x1*x2*zp*zp - 2*x1*x1_outer*z2*z2_outer + 2*x1*x1_outer*z2*zp + 2*x1*x1_outer*z2_outer*zp - 2*x1*x1_outer*zp*zp - 2*x1*x2_outer*z1*z2_outer + 2*x1*x2_outer*z1*zp + 4*x1*x2_outer*z2*z1_outer - 4*x1*x2_outer*z2*zp - 4*x1*x2_outer*z1_outer*zp + 2*x1*x2_outer*z2_outer*zp + 2*x1*x2_outer*zp*zp + 2*x1*xp*z1*z2_outer - 2*x1*xp*z1*zp - 4*x1*xp*z2*z1_outer + 2*x1*xp*z2*z2_outer + 2*x1*xp*z2*zp + 2*x1*xp*z1_outer*z2_outer + 2*x1*xp*z1_outer*zp - 2*x1*xp*z2_outer*z2_outer - 2*x1*xp*z2_outer*zp + x2*x2*z1_outer*z1_outer - 2*x2*x2*z1_outer*zp + x2*x2*zp*zp + 4*x2*x1_outer*z1*z2_outer - 4*x2*x1_outer*z1*zp - 2*x2*x1_outer*z2*z1_outer + 2*x2*x1_outer*z2*zp + 2*x2*x1_outer*z1_outer*zp - 4*x2*x1_outer*z2_outer*zp + 2*x2*x1_outer*zp*zp - 2*x2*x2_outer*z1*z1_outer + 2*x2*x2_outer*z1*zp + 2*x2*x2_outer*z1_outer*zp - 2*x2*x2_outer*zp*zp + 2*x2*xp*z1*z1_outer - 4*x2*xp*z1*z2_outer + 2*x2*xp*z1*zp + 2*x2*xp*z2*z1_outer - 2*x2*xp*z2*zp - 2*x2*xp*z1_outer*z1_outer + 2*x2*xp*z1_outer*z2_outer - 2*x2*xp*z1_outer*zp + 2*x2*xp*z2_outer*zp + x1_outer*x1_outer*z2*z2 - 2*x1_outer*x1_outer*z2*zp + x1_outer*x1_outer*zp*zp - 2*x1_outer*x2_outer*z1*z2 + 2*x1_outer*x2_outer*z1*zp + 2*x1_outer*x2_outer*z2*zp - 2*x1_outer*x2_outer*zp*zp + 2*x1_outer*xp*z1*z2 - 4*x1_outer*xp*z1*z2_outer + 2*x1_outer*xp*z1*zp - 2*x1_outer*xp*z2*z2 + 2*x1_outer*xp*z2*z1_outer + 2*x1_outer*xp*z2*z2_outer - 2*x1_outer*xp*z2*zp - 2*x1_outer*xp*z1_outer*zp + 2*x1_outer*xp*z2_outer*zp + x2_outer*x2_outer*z1*z1 - 2*x2_outer*x2_outer*z1*zp + x2_outer*x2_outer*zp*zp - 2*x2_outer*xp*z1*z1 + 2*x2_outer*xp*z1*z2 + 2*x2_outer*xp*z1*z1_outer + 2*x2_outer*xp*z1*z2_outer - 2*x2_outer*xp*z1*zp - 4*x2_outer*xp*z2*z1_outer + 2*x2_outer*xp*z2*zp + 2*x2_outer*xp*z1_outer*zp - 2*x2_outer*xp*z2_outer*zp + xp*xp*z1*z1 - 2*xp*xp*z1*z2 - 2*xp*xp*z1*z1_outer + 2*xp*xp*z1*z2_outer + xp*xp*z2*z2 + 2*xp*xp*z2*z1_outer - 2*xp*xp*z2*z2_outer + xp*xp*z1_outer*z1_outer - 2*xp*xp*z1_outer*z2_outer + xp*xp*z2_outer*z2_outer) + z2_outer*sqrt(x1*x1*z2_outer*z2_outer - 2*x1*x1*z2_outer*zp + x1*x1*zp*zp - 2*x1*x2*z1_outer*z2_outer + 2*x1*x2*z1_outer*zp + 2*x1*x2*z2_outer*zp - 2*x1*x2*zp*zp - 2*x1*x1_outer*z2*z2_outer + 2*x1*x1_outer*z2*zp + 2*x1*x1_outer*z2_outer*zp - 2*x1*x1_outer*zp*zp - 2*x1*x2_outer*z1*z2_outer + 2*x1*x2_outer*z1*zp + 4*x1*x2_outer*z2*z1_outer - 4*x1*x2_outer*z2*zp - 4*x1*x2_outer*z1_outer*zp + 2*x1*x2_outer*z2_outer*zp + 2*x1*x2_outer*zp*zp + 2*x1*xp*z1*z2_outer - 2*x1*xp*z1*zp - 4*x1*xp*z2*z1_outer + 2*x1*xp*z2*z2_outer + 2*x1*xp*z2*zp + 2*x1*xp*z1_outer*z2_outer + 2*x1*xp*z1_outer*zp - 2*x1*xp*z2_outer*z2_outer - 2*x1*xp*z2_outer*zp + x2*x2*z1_outer*z1_outer - 2*x2*x2*z1_outer*zp + x2*x2*zp*zp + 4*x2*x1_outer*z1*z2_outer - 4*x2*x1_outer*z1*zp - 2*x2*x1_outer*z2*z1_outer + 2*x2*x1_outer*z2*zp + 2*x2*x1_outer*z1_outer*zp - 4*x2*x1_outer*z2_outer*zp + 2*x2*x1_outer*zp*zp - 2*x2*x2_outer*z1*z1_outer + 2*x2*x2_outer*z1*zp + 2*x2*x2_outer*z1_outer*zp - 2*x2*x2_outer*zp*zp + 2*x2*xp*z1*z1_outer - 4*x2*xp*z1*z2_outer + 2*x2*xp*z1*zp + 2*x2*xp*z2*z1_outer - 2*x2*xp*z2*zp - 2*x2*xp*z1_outer*z1_outer + 2*x2*xp*z1_outer*z2_outer - 2*x2*xp*z1_outer*zp + 2*x2*xp*z2_outer*zp + x1_outer*x1_outer*z2*z2 - 2*x1_outer*x1_outer*z2*zp + x1_outer*x1_outer*zp*zp - 2*x1_outer*x2_outer*z1*z2 + 2*x1_outer*x2_outer*z1*zp + 2*x1_outer*x2_outer*z2*zp - 2*x1_outer*x2_outer*zp*zp + 2*x1_outer*xp*z1*z2 - 4*x1_outer*xp*z1*z2_outer + 2*x1_outer*xp*z1*zp - 2*x1_outer*xp*z2*z2 + 2*x1_outer*xp*z2*z1_outer + 2*x1_outer*xp*z2*z2_outer - 2*x1_outer*xp*z2*zp - 2*x1_outer*xp*z1_outer*zp + 2*x1_outer*xp*z2_outer*zp + x2_outer*x2_outer*z1*z1 - 2*x2_outer*x2_outer*z1*zp + x2_outer*x2_outer*zp*zp - 2*x2_outer*xp*z1*z1 + 2*x2_outer*xp*z1*z2 + 2*x2_outer*xp*z1*z1_outer + 2*x2_outer*xp*z1*z2_outer - 2*x2_outer*xp*z1*zp - 4*x2_outer*xp*z2*z1_outer + 2*x2_outer*xp*z2*zp + 2*x2_outer*xp*z1_outer*zp - 2*x2_outer*xp*z2_outer*zp + xp*xp*z1*z1 - 2*xp*xp*z1*z2 - 2*xp*xp*z1*z1_outer + 2*xp*xp*z1*z2_outer + xp*xp*z2*z2 + 2*xp*xp*z2*z1_outer - 2*xp*xp*z2*z2_outer + xp*xp*z1_outer*z1_outer - 2*xp*xp*z1_outer*z2_outer + xp*xp*z2_outer*z2_outer) + x1*z1_outer*z2_outer - 2*x1_outer*z1*z2_outer + x1_outer*z2*z1_outer + x2_outer*z1*z1_outer + x2*z1_outer*z2_outer + x1_outer*z2*z2_outer + x2_outer*z1*z2_outer - 2*x2_outer*z2*z1_outer + x1*z1_outer*zp - xp*z1*z1_outer - x1*z2_outer*zp - x2*z1_outer*zp + xp*z1*z2_outer + xp*z2*z1_outer + x2*z2_outer*zp - x1_outer*z1_outer*zp - xp*z2*z2_outer + x1_outer*z2_outer*zp + x2_outer*z1_outer*zp - 2*xp*z1_outer*z2_outer - x2_outer*z2_outer*zp)/(2*(x1*z1_outer - x1_outer*z1 - x1*z2_outer - x2*z1_outer + x1_outer*z2 + x2_outer*z1 + x2*z2_outer - x2_outer*z2)))/(z1_outer - z2_outer);
				
				ratio_pt_wrt_l = 1 - ratio_pt_wrt_r;
				
				ratio_pt_wrt_b = (z2 - (xp*z2*z2 - x1_outer*z2*z2 - x1*z2_outer*z2_outer + xp*z2_outer*z2_outer + z2*sqrt(x1*x1*z2_outer*z2_outer - 2*x1*x1*z2_outer*zp + x1*x1*zp*zp - 2*x1*x2*z1_outer*z2_outer + 2*x1*x2*z1_outer*zp + 2*x1*x2*z2_outer*zp - 2*x1*x2*zp*zp - 2*x1*x1_outer*z2*z2_outer + 2*x1*x1_outer*z2*zp + 2*x1*x1_outer*z2_outer*zp - 2*x1*x1_outer*zp*zp - 2*x1*x2_outer*z1*z2_outer + 2*x1*x2_outer*z1*zp + 4*x1*x2_outer*z2*z1_outer - 4*x1*x2_outer*z2*zp - 4*x1*x2_outer*z1_outer*zp + 2*x1*x2_outer*z2_outer*zp + 2*x1*x2_outer*zp*zp + 2*x1*xp*z1*z2_outer - 2*x1*xp*z1*zp - 4*x1*xp*z2*z1_outer + 2*x1*xp*z2*z2_outer + 2*x1*xp*z2*zp + 2*x1*xp*z1_outer*z2_outer + 2*x1*xp*z1_outer*zp - 2*x1*xp*z2_outer*z2_outer - 2*x1*xp*z2_outer*zp + x2*x2*z1_outer*z1_outer - 2*x2*x2*z1_outer*zp + x2*x2*zp*zp + 4*x2*x1_outer*z1*z2_outer - 4*x2*x1_outer*z1*zp - 2*x2*x1_outer*z2*z1_outer + 2*x2*x1_outer*z2*zp + 2*x2*x1_outer*z1_outer*zp - 4*x2*x1_outer*z2_outer*zp + 2*x2*x1_outer*zp*zp - 2*x2*x2_outer*z1*z1_outer + 2*x2*x2_outer*z1*zp + 2*x2*x2_outer*z1_outer*zp - 2*x2*x2_outer*zp*zp + 2*x2*xp*z1*z1_outer - 4*x2*xp*z1*z2_outer + 2*x2*xp*z1*zp + 2*x2*xp*z2*z1_outer - 2*x2*xp*z2*zp - 2*x2*xp*z1_outer*z1_outer + 2*x2*xp*z1_outer*z2_outer - 2*x2*xp*z1_outer*zp + 2*x2*xp*z2_outer*zp + x1_outer*x1_outer*z2*z2 - 2*x1_outer*x1_outer*z2*zp + x1_outer*x1_outer*zp*zp - 2*x1_outer*x2_outer*z1*z2 + 2*x1_outer*x2_outer*z1*zp + 2*x1_outer*x2_outer*z2*zp - 2*x1_outer*x2_outer*zp*zp + 2*x1_outer*xp*z1*z2 - 4*x1_outer*xp*z1*z2_outer + 2*x1_outer*xp*z1*zp - 2*x1_outer*xp*z2*z2 + 2*x1_outer*xp*z2*z1_outer + 2*x1_outer*xp*z2*z2_outer - 2*x1_outer*xp*z2*zp - 2*x1_outer*xp*z1_outer*zp + 2*x1_outer*xp*z2_outer*zp + x2_outer*x2_outer*z1*z1 - 2*x2_outer*x2_outer*z1*zp + x2_outer*x2_outer*zp*zp - 2*x2_outer*xp*z1*z1 + 2*x2_outer*xp*z1*z2 + 2*x2_outer*xp*z1*z1_outer + 2*x2_outer*xp*z1*z2_outer - 2*x2_outer*xp*z1*zp - 4*x2_outer*xp*z2*z1_outer + 2*x2_outer*xp*z2*zp + 2*x2_outer*xp*z1_outer*zp - 2*x2_outer*xp*z2_outer*zp + xp*xp*z1*z1 - 2*xp*xp*z1*z2 - 2*xp*xp*z1*z1_outer + 2*xp*xp*z1*z2_outer + xp*xp*z2*z2 + 2*xp*xp*z2*z1_outer - 2*xp*xp*z2*z2_outer + xp*xp*z1_outer*z1_outer - 2*xp*xp*z1_outer*z2_outer + xp*xp*z2_outer*z2_outer) - z2_outer*sqrt(x1*x1*z2_outer*z2_outer - 2*x1*x1*z2_outer*zp + x1*x1*zp*zp - 2*x1*x2*z1_outer*z2_outer + 2*x1*x2*z1_outer*zp + 2*x1*x2*z2_outer*zp - 2*x1*x2*zp*zp - 2*x1*x1_outer*z2*z2_outer + 2*x1*x1_outer*z2*zp + 2*x1*x1_outer*z2_outer*zp - 2*x1*x1_outer*zp*zp - 2*x1*x2_outer*z1*z2_outer + 2*x1*x2_outer*z1*zp + 4*x1*x2_outer*z2*z1_outer - 4*x1*x2_outer*z2*zp - 4*x1*x2_outer*z1_outer*zp + 2*x1*x2_outer*z2_outer*zp + 2*x1*x2_outer*zp*zp + 2*x1*xp*z1*z2_outer - 2*x1*xp*z1*zp - 4*x1*xp*z2*z1_outer + 2*x1*xp*z2*z2_outer + 2*x1*xp*z2*zp + 2*x1*xp*z1_outer*z2_outer + 2*x1*xp*z1_outer*zp - 2*x1*xp*z2_outer*z2_outer - 2*x1*xp*z2_outer*zp + x2*x2*z1_outer*z1_outer - 2*x2*x2*z1_outer*zp + x2*x2*zp*zp + 4*x2*x1_outer*z1*z2_outer - 4*x2*x1_outer*z1*zp - 2*x2*x1_outer*z2*z1_outer + 2*x2*x1_outer*z2*zp + 2*x2*x1_outer*z1_outer*zp - 4*x2*x1_outer*z2_outer*zp + 2*x2*x1_outer*zp*zp - 2*x2*x2_outer*z1*z1_outer + 2*x2*x2_outer*z1*zp + 2*x2*x2_outer*z1_outer*zp - 2*x2*x2_outer*zp*zp + 2*x2*xp*z1*z1_outer - 4*x2*xp*z1*z2_outer + 2*x2*xp*z1*zp + 2*x2*xp*z2*z1_outer - 2*x2*xp*z2*zp - 2*x2*xp*z1_outer*z1_outer + 2*x2*xp*z1_outer*z2_outer - 2*x2*xp*z1_outer*zp + 2*x2*xp*z2_outer*zp + x1_outer*x1_outer*z2*z2 - 2*x1_outer*x1_outer*z2*zp + x1_outer*x1_outer*zp*zp - 2*x1_outer*x2_outer*z1*z2 + 2*x1_outer*x2_outer*z1*zp + 2*x1_outer*x2_outer*z2*zp - 2*x1_outer*x2_outer*zp*zp + 2*x1_outer*xp*z1*z2 - 4*x1_outer*xp*z1*z2_outer + 2*x1_outer*xp*z1*zp - 2*x1_outer*xp*z2*z2 + 2*x1_outer*xp*z2*z1_outer + 2*x1_outer*xp*z2*z2_outer - 2*x1_outer*xp*z2*zp - 2*x1_outer*xp*z1_outer*zp + 2*x1_outer*xp*z2_outer*zp + x2_outer*x2_outer*z1*z1 - 2*x2_outer*x2_outer*z1*zp + x2_outer*x2_outer*zp*zp - 2*x2_outer*xp*z1*z1 + 2*x2_outer*xp*z1*z2 + 2*x2_outer*xp*z1*z1_outer + 2*x2_outer*xp*z1*z2_outer - 2*x2_outer*xp*z1*zp - 4*x2_outer*xp*z2*z1_outer + 2*x2_outer*xp*z2*zp + 2*x2_outer*xp*z1_outer*zp - 2*x2_outer*xp*z2_outer*zp + xp*xp*z1*z1 - 2*xp*xp*z1*z2 - 2*xp*xp*z1*z1_outer + 2*xp*xp*z1*z2_outer + xp*xp*z2*z2 + 2*xp*xp*z2*z1_outer - 2*xp*xp*z2*z2_outer + xp*xp*z1_outer*z1_outer - 2*xp*xp*z1_outer*z2_outer + xp*xp*z2_outer*z2_outer) + x1*z2*z2_outer - 2*x2*z1*z2_outer + x2*z2*z1_outer + x2_outer*z1*z2 + x2*z1_outer*z2_outer + x1_outer*z2*z2_outer + x2_outer*z1*z2_outer - 2*x2_outer*z2*z1_outer + x1*z2*zp - xp*z1*z2 - x2*z2*zp - x1*z2_outer*zp - x1_outer*z2*zp + xp*z1*z2_outer + xp*z2*z1_outer + x2*z2_outer*zp + x2_outer*z2*zp - 2*xp*z2*z2_outer + x1_outer*z2_outer*zp - xp*z1_outer*z2_outer - x2_outer*z2_outer*zp)/(2*(x1*z2 - x2*z1 - x1*z2_outer + x2*z1_outer - x1_outer*z2 + x2_outer*z1 + x1_outer*z2_outer - x2_outer*z1_outer)))/(z2 - z2_outer);
				
				ratio_pt_wrt_t = 1 - ratio_pt_wrt_b;
				
				/*
				cout << "xp,yp,zp,ratio_pt_wrt_l,ratio_pt_wrt_r,ratio_pt_wrt_t,ratio_pt_wrt_b\n";
				
				cout << xp << '\n';
				cout << zp << '\n';
				cout << ratio_pt_wrt_l << '\n';
				cout << ratio_pt_wrt_r << '\n';
				cout << ratio_pt_wrt_t << '\n';
				cout << ratio_pt_wrt_b << '\n';
				*/
				
				
				
				
				x1_1_outer = ratio_pt_wrt_t*x1 + ratio_pt_wrt_b*x1_outer_cur;
				z1_1_outer = ratio_pt_wrt_t*z1 + ratio_pt_wrt_b*z1_outer_cur;
				
				x2_2_outer = ratio_pt_wrt_t*x2 + ratio_pt_wrt_b*x2_outer_cur;
				z2_2_outer = ratio_pt_wrt_t*z2 + ratio_pt_wrt_b*z2_outer_cur;
				
				x12 = ratio_pt_wrt_r*x1 + ratio_pt_wrt_l*x2;
				z12 = ratio_pt_wrt_r*z1 + ratio_pt_wrt_l*z2;
				
				x12_outer = ratio_pt_wrt_r*x1_outer_cur + ratio_pt_wrt_l*x2_outer_cur;
				z12_outer = ratio_pt_wrt_r*z1_outer_cur + ratio_pt_wrt_l*z2_outer_cur;
				
				dd = (x1_1_outer - x2_2_outer)*(z12 - z12_outer) - (z1_1_outer - z2_2_outer)*(x12 - x12_outer);
				
				xp_new = ((x1_1_outer*z2_2_outer - z1_1_outer*x2_2_outer)*(x12 - x12_outer) - (x1_1_outer - x2_2_outer)*(x12*z12_outer - z12*x12_outer))/dd;
				yp_new = yp;
				zp_new = ((x1_1_outer*z2_2_outer - z1_1_outer*x2_2_outer)*(z12 - z12_outer) - (z1_1_outer - z2_2_outer)*(x12*z12_outer - z12*x12_outer))/dd;
				
				//rotate /pitch - no need to
				
				/*xp_new2 = (xp_new - pitch_origin[0])*cos(orient_pitch_cur) - (yp_new - pitch_origin[1])*sin(orient_pitch_cur) + pitch_origin[0];
				yp_new2 = (xp_new - pitch_origin[0])*sin(orient_pitch_cur) + (yp_new - pitch_origin[1])*cos(orient_pitch_cur) + pitch_origin[1];
				zp_new2 = zp_new;*/
				
				xp_new2 = xp_new;
				yp_new2 = yp_new;
				zp_new2 = zp_new;
				
				sd[iter].component(vector::X) = xp_new2 - x_init;
					
				sd[iter].component(vector::Y) = yp_new2 - y_init;
					
				sd[iter].component(vector::Z) = zp_new2 - z_init;
				
				/*
				cout << "x1_1_outer,x2_2_outer,x12,x12_outer,sd x z\n";
				
				cout << x1_1_outer << '\n' << z1_1_outer << '\n';
				cout << x2_2_outer << '\n' << z2_2_outer << '\n';
				cout << x12 << '\n' << z12 << '\n';
				cout << x12_outer << '\n' << z12_outer << '\n';
				cout << sd[iter].component(vector::X) << '\n' << sd[iter].component(vector::Z) << '\n';
				
				getchar();
				*/
		
				
			//red circle region
				
			} else if (which_side_outer > 0) {
				
				//sd0Rel.component(vector::X) = (xp - x1_outer);
				//sd0Rel.component(vector::Y) = 0.;
				//sd0Rel.component(vector::Z) = (zp - z1_outer);
				
				//sd[iter] = (sd0Rel*(cos(small_angle) - 1) + (axisHat ^ sd0Rel*sin(small_angle)) + (axisHat & sd0Rel)*(1 - cos(small_angle))*axisHat) + outer1_delta;
				
				xp_new = x1_outer_cur - x1_outer + xp;
				yp_new = yp;
				zp_new = z1_outer_cur - z1_outer + zp;
				
				//rotate /pitch - no need to
				
				/*xp_new2 = (xp_new - pitch_origin[0])*cos(orient_pitch_cur) - (yp_new - pitch_origin[1])*sin(orient_pitch_cur) + pitch_origin[0];
				yp_new2 = (xp_new - pitch_origin[0])*sin(orient_pitch_cur) + (yp_new - pitch_origin[1])*cos(orient_pitch_cur) + pitch_origin[1];
				zp_new2 = zp_new;*/
				xp_new2 = xp_new;
				yp_new2 = yp_new;
				zp_new2 = zp_new;
				
				sd[iter].component(vector::X) = xp_new2 - x_init;
					
				sd[iter].component(vector::Y) = yp_new2 - y_init;
					
				sd[iter].component(vector::Z) = zp_new2 - z_init;
				
				
				/* outer1_delta.component(vector::X) = x1_outer_cur - x1_outer;
				outer1_delta.component(vector::Y) = 0.;
				outer1_delta.component(vector::Z) = z1_outer_cur - z1_outer;
				
				
				sd[iter] = outer1_delta; */
					
			} else {
				
				//print out pts to check its location, see whether suitable for orientation reference
				/* if (iter%100 == 0) {
					
					cout << "iter,xp,yp,zp\n";
					
					cout << iter << '\n';
					cout << xp << '\n';
					cout << yp << '\n';
					cout << zp << '\n';
					
				} */

				
				sd[iter] = {0., 0., 0.};
					
					
			}
			
			//pitching rotation
			
			/* scalar new_xp = (sd[iter].component(vector::Y) - pitch_origin.component(vector::Y))*sin(pitch_angle) + (sd[iter].component(vector::X) -pitch_origin.component(vector::X))*cos(pitch_angle) + pitch_origin.component(vector::X);
			scalar new_yp = (sd[iter].component(vector::Y) - pitch_origin.component(vector::Y))*cos(pitch_angle) - (sd[iter].component(vector::X) -pitch_origin.component(vector::X))*sin(pitch_angle) + pitch_origin.component(vector::Y);
			
			tmp.component(vector::X) = new_xp - sd[iter].component(vector::X);
			tmp.component(vector::Y) = new_yp - sd[iter].component(vector::Y);
			tmp.component(vector::Z) = 0.; */
			
			//sd0Rel = p0_[iter] - pitch_origin;
			//sd0Rel = sd[iter] - pitch_origin;
			
			//tmp = sd0Rel*(cos(pitch_angle) - 1) + (pitch_axisHat ^ sd0Rel*sin(pitch_angle)) + (pitch_axisHat & sd0Rel)*(1 - cos(pitch_angle))*pitch_axisHat;
			
			//sd[iter] = tmp;
			
			//sd[iter] = sd[iter] + tmp;
			
			/* sd[iter].component(vector::X) = 0.;
			sd[iter].component(vector::Y) = 0.;
			sd[iter].component(vector::Z) = 0.; */
				
			
		};
		
		saved_time_ = t.value();
		
		vectorField::operator=
		(
			
			sd

		);

		fixedValuePointPatchField<vector>::updateCoeffs();
		
	}
}



void prototype_retract_motion_gen::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("angle0", angle0_);
	os.writeEntry("x1", x1);
	os.writeEntry("z1", z1);
	os.writeEntry("x1_outer", x1_outer);
	os.writeEntry("z1_outer", z1_outer);
	os.writeEntry("x2", x2);
	os.writeEntry("z2", z2);
	os.writeEntry("x2_outer", x2_outer);
	os.writeEntry("z2_outer", z2_outer);
	os.writeEntry("amplitude", amplitude_);
	os.writeEntry("maximum_retract_angle", max_retract_angle_);
    os.writeEntry("omega", omega_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    prototype_retract_motion_gen
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
