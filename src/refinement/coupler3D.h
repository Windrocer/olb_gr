/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Zhishang Xu 
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef REFINEMENT_COUPLER_3D_H
#define REFINEMENT_COUPLER_3D_H

#include "grid3D.h"

namespace olb {

/// General coupling structre
template <typename T, typename DESCRIPTOR>
class Coupler3D {
protected:
	Grid3D<T,DESCRIPTOR>& _coarse;
	Grid3D<T,DESCRIPTOR>& _fine;

	const int _coarseSize_a;
	const int _coarseSize_b;
	const int _fineSize_a;
	const int _fineSize_b;

	// const bool _vertical ---> will not use in 3D 
	// instead the direction will be judged according to extend vector
	
	const Vector<T,3> _physOrigin;

	const Vector<int,4>& getFineLatticeR(int x1, int x2) const;
	const Vector<int,4>& getCoarseLatticeR(int x1, int x2) const;

	T getScalingFactor() const;
	T getInvScalingFactor() const;

private:
	std::vector<std::vector<Vector<int,4>>> _coarseLatticeR;
	std::vector<std::vector<Vector<int,4>>> _fineLatticeR;

public:
	Coupler3D(Grid3D<T,DESCRIPTOR>& coarse, Grid3D<T,DESCRIPTOR>& fine,
			  Vector<T,3> origin, Vector<T,3> extend);
	// Parameter overlap is for bulk initialization
	Coupler3D(Grid3D<T,DESCRIPTOR>& coarse, Grid3D<T,DESCRIPTOR>& fine,
			  Vector<T,3> origin, Vector<T,3> extend, int overlap);
};

/// FineCoupler3D is for coupling from coarse to fine. 
//  Interpolation is needed.
template <typename T, typename DESCRIPTOR>
class FineCoupler3D : public Coupler3D<T,DESCRIPTOR> {
private:
	std::vector<std::vector<Vector<T,1>>>				_c2f_rho;
	std::vector<std::vector<Vector<T,DESCRIPTOR::d>>>	_c2f_u;
	std::vector<std::vector<Vector<T,DESCRIPTOR::q>>>	_c2f_fneq;

public:
	FineCoupler3D(Grid3D<T,DESCRIPTOR>& coarse, Grid3D<T,DESCRIPTOR>& fine,
				  Vector<T,3> origin, Vector<T,3> extend);

	// self-functions
	void store();
	void interpolate();
	void couple();
};

/// BulkFineCoupler3D is specific for in-stream nodes coupling
template <typename T, typename DESCRIPTOR>
class BulkFineCoupler3D : public Coupler3D<T,DESCRIPTOR> {
private:
	std::vector<std::vector<Vector<T,1>>>				_c2f_rho;
	std::vector<std::vector<Vector<T,DESCRIPTOR::d>>>	_c2f_u;
	std::vector<std::vector<Vector<T,DESCRIPTOR::q>>>	_c2f_fneq;
	const int _overlap;

public:
	BulkFineCoupler3D(Grid3D<T,DESCRIPTOR>& coarse, Grid3D<T,DESCRIPTOR>& fine,
				  Vector<T,3> origin, Vector<T,3> extend, int overlap);

	// self-functions
	void store();
	void interpolate();
	void couple();
};

/// CoarseCoupler3D is for coupling from fine to coarse
template <typename T, typename DESCRIPTOR>
class CoarseCoupler3D : public Coupler3D<T,DESCRIPTOR> {
public:
	CoarseCoupler3D(Grid3D<T,DESCRIPTOR>& coarse, Grid3D<T,DESCRIPTOR>& fine,
					Vector<T,3> origin, Vector<T,3> extend);

	// from fine to coarse, only couple is needed --- no interpolation!
	void couple();
};

}
#endif
