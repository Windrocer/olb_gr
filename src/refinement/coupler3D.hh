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

#ifndef REFINEMENT_COUPLER_3D_HH
#define REFINEMENT_COUPLER_3D_HH

#include "coupler3D.h"

#include "dynamics/lbHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
Coupler3D<T,DESCRIPTOR>::Coupler3D(Grid3D<T,DESCRIPTOR>& coarse, 
		Grid3D<T,DESCRIPTOR>& fine, Vector<T,3> origin, Vector<T,3> extend):
	_coarse(coarse),
	_fine(fine),
	// if extend[0] is 0, then y-z, use extend[1]
	// lx / dx + 1 is discretized total points
	_coarseSize_a(std::floor((util::nearZero(extend[0]) ? extend[1] :
		extend[0]) / coarse.getConverter().getPhysDeltaX() + 0.5) + 1),
	// if extend[0] or extend[1] is 0, then must be x-z or y-z
	// if extend[0] or extend[1] is not 0, then must be x-y
	_coarseSize_b(std::floor((util::nearZero(extend[0] * extend[1]) ? 
			extend[2] : extend[1]) / coarse.getConverter().getPhysDeltaX() 
			+ 0.5) + 1),
	_fineSize_a(2*_coarseSize_a - 1),
	_fineSize_b(2*_coarseSize_b - 1),
	_physOrigin(_coarse.alignOriginToGrid(origin)),
	_coarseLatticeR(_coarseSize_a, 
			std::vector<Vector<int,4>>(_coarseSize_b)),
	_fineLatticeR(_fineSize_a, std::vector<Vector<int,4>>(_fineSize_b))
{
	OLB_ASSERT( util::nearZero(extend[0]) || util::nearZero(extend[1]) || 
			util::nearZero(extend[2]), "Coupling is only implemented alongside unit plane vectors.");

	// Get cuboid geometry to align the matrix position to lattice R
	const auto& coarseGeometry	= _coarse.getCuboidGeometry();
	const auto& fineGeometry	= _fine.getCuboidGeometry();
	// step size based on fine grid
	const T deltaX = _fine.getConverter().getPhysDeltaX();
	const Vector<T,3> stepPhysR_a = (util::nearZero(extend[0])) ? 
									Vector<T,3> {0, deltaX, 0} :
									Vector<T,3> {deltaX, 0, 0};
	const Vector<T,3> stepPhysR_b = (util::nearZero(extend[0] * 
				extend[1])) ? Vector<T,3> {0, 0, deltaX} : 
							Vector<T,3> {0, deltaX, 0};
	// fulfil latticeR matrix
	for (int i = 0; i < _fineSize_a; ++i) {
		for (int j = 0; j < _fineSize_b; ++j) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				coarseGeometry.getLatticeR(_physOrigin + i*stepPhysR_a + 
						j*stepPhysR_b, _coarseLatticeR[i/2][j/2]);
			}
			fineGeometry.getLatticeR(_physOrigin + i*stepPhysR_a + 
					j*stepPhysR_b, _fineLatticeR[i][j]);
		}
	}
}

template <typename T, typename DESCRIPTOR>
Coupler3D<T,DESCRIPTOR>::Coupler3D(Grid3D<T,DESCRIPTOR>& coarse, 
		Grid3D<T,DESCRIPTOR>& fine, Vector<T,3> origin, Vector<T,3> extend,
		int overlap):
	_coarse(coarse),
	_fine(fine),
	_coarseSize_a(std::floor((util::nearZero(extend[0]) ? extend[1] :
					extend[0]) / coarse.getConverter().getPhysDeltaX() + 0.5) + 1),
	_coarseSize_b(std::floor((util::nearZero(extend[0]*extend[1]) ? extend[2] :
					extend[1]) / coarse.getConverter().getPhysDeltaX() + 0.5) + 1),
	_fineSize_a(2*_coarseSize_a - 1),
	_fineSize_b(2*_coarseSize_b - 1),
	_physOrigin(_coarse.alignOriginToGrid(origin)),
	_coarseLatticeR(_coarseSize_a, std::vector<Vector<int,4>>(_coarseSize_b)),
	_fineLatticeR(_fineSize_a, std::vector<Vector<int,4>>(_fineSize_b))
{
	OLB_ASSERT( util::nearZero(extend[0]) || util::nearZero(extend[1]) || 
			util::nearZero(extend[2]), "Coupling is only implemented alongside unit plane vectors.");

	// Get cuboid geometry to align the matrix position to lattice R
	const auto& coarseGeometry	= _coarse.getCuboidGeometry();
	const auto& fineGeometry	= _fine.getCuboidGeometry();
	// step size based on fine grid
	const T deltaX = _fine.getConverter().getPhysDeltaX();
	const Vector<T,3> stepPhysR_a = (util::nearZero(extend[0])) ? 
									Vector<T,3> {0, deltaX, 0} :
									Vector<T,3> {deltaX, 0, 0};
	const Vector<T,3> stepPhysR_b = (util::nearZero(extend[0] * extend[1])) ?
									Vector<T,3> {0, 0, deltaX} :
									Vector<T,3> {0, deltaX, 0};
	// fulfil latticeR matrix
	for (int i = 0; i < _fineSize_a; ++i) {
		for (int j = 0; j < _fineSize_b; ++j) {
			if ((i % 2 == 0) && (j % 2 == 0)) {
				coarseGeometry.getLatticeR(_physOrigin + i*stepPhysR_a + j*stepPhysR_b,
									   _coarseLatticeR[i/2][j/2]);
			}
			if ((i >= overlap) && (i < _fineSize_a - overlap) &&
				(j >= overlap) && (j < _fineSize_b - overlap)) {
				fineGeometry.getLatticeR(_physOrigin + i*stepPhysR_a + j*stepPhysR_b,
										 _fineLatticeR[i][j]);
			}
			else {
				for (int m = 0; m < 4; ++m) {
					_fineLatticeR[i][j][m] = 0;
				}
			}
		}
	}
}

template <typename T, typename DESCRIPTOR>
T Coupler3D<T,DESCRIPTOR>::getScalingFactor() const
{
	const T coarseTau = _coarse.getConverter().getLatticeRelaxationTime();
	return (coarseTau - 0.25) / coarseTau;
}

template <typename T, typename DESCRIPTOR>
T Coupler3D<T,DESCRIPTOR>::getInvScalingFactor() const
{
	return 1./getScalingFactor();
}

template <typename T, typename DESCRIPTOR>
const Vector<int,4>& Coupler3D<T,DESCRIPTOR>::getFineLatticeR(int x1, int x2) const
{
	return _fineLatticeR[x1][x2];
}

template <typename T, typename DESCRIPTOR>
const Vector<int,4>& Coupler3D<T,DESCRIPTOR>::getCoarseLatticeR(int x1, int x2) const
{
	return _coarseLatticeR[x1][x2];
}

template <typename T, typename DESCRIPTOR>
FineCoupler3D<T,DESCRIPTOR>::FineCoupler3D(Grid3D<T,DESCRIPTOR>& coarse, 
		Grid3D<T,DESCRIPTOR>& fine, Vector<T,3> origin, Vector<T,3> extend):
	Coupler3D<T,DESCRIPTOR>(coarse, fine, origin, extend),
	_c2f_rho(this->_coarseSize_a, std::vector<Vector<T,1>>(this->_coarseSize_b)),
	_c2f_u(this->_coarseSize_a, std::vector<Vector<T,DESCRIPTOR::d>>(this->_coarseSize_b)),
//				Vector<T,DESCRIPTOR::d>(T{}))),
	_c2f_fneq(this->_coarseSize_a, std::vector<Vector<T,DESCRIPTOR::q>>(this->_coarseSize_b))
//				Vector<T,DESCRIPTOR::q>(T{})))
{
	OstreamManager clout(std::cout, "C2F");

	const auto& coarseOrigin = this->getCoarseLatticeR(0,0);
	const auto& fineOrigin = this->getFineLatticeR(0,0);

	clout << "coarse origin: " << coarseOrigin[0] << " " << coarseOrigin[1] << " " 
		  << coarseOrigin[2] << " " << coarseOrigin[3] << std::endl;
	clout << "coarse size:	" << this->_coarseSize_a << " * " << this->_coarseSize_b
		  << std::endl;
	clout << "fine origin: " << fineOrigin[0] << " " << fineOrigin[1] << " "
		  << fineOrigin[2] << " " << fineOrigin[3] << std::endl;
	clout << "fine size:	" << this->_fineSize_a << " * " << this->_fineSize_b
		  << std::endl;
}

template <typename T, typename DESCRIPTOR>
void FineCoupler3D<T,DESCRIPTOR>::store()
{
	auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
			const auto pos = this->getCoarseLatticeR(i, j);
			T rho{};
			T u[DESCRIPTOR::d] {};
			T fNeq[DESCRIPTOR::q] {};
			T pi[util::TensorVal<DESCRIPTOR>::n] {};
			Cell<T,DESCRIPTOR> coarseCell;
			coarseLattice.get(pos, coarseCell);
			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
//			coarseCell.computeRhoU(rho, u);
			lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);
//			coarseCell.computeAllMomenta(rho, u, pi);
//			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
//				fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
//			}
			
			_c2f_rho[i][j] = Vector<T,1>(rho);
//			for ( int id = 0; id < DESCRIPTOR::d; ++id ) {
//				if ( util::nearZero(u[id], 1.0e-15) ) {
//					u[id] = 0;
//				}
//			}
			_c2f_u[i][j] = Vector<T,DESCRIPTOR::d>(u);
//			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
//			lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);
			_c2f_fneq[i][j] = Vector<T,DESCRIPTOR::q>(fNeq);
		}
	}
}

template <typename T, unsigned N>
Vector<T,N> order2interpolation(const Vector<T,N>& f0, const Vector<T,N>& f1)
{
	return 0.5 * (f0 + f1);
}

template <typename T, unsigned N>
Vector<T,N> order2interpolation(const std::vector<std::vector<Vector<T,N>>>& data, int i, int j)
{
	return 0.25 * (data[i][j] + data[i+1][j] + data[i][j+1] + data[i+1][j+1]);
}

// order 3 interpolation in 1D --- will use order 4 instead
template <typename T, unsigned N>
Vector<T,N> order3interpolation(const std::vector<std::vector<Vector<T,N>>>& data, int i, int j,
		bool vertical, bool ascending)
{
	if (!vertical) {
		if (ascending) {
			return 3./8. * data[i][j] + 3./4. * data[i+1][j] - 1./8. * data[i+2][j];
		}
		else {
			return 3./8. * data[i+1][j] + 3./4. * data[i][j] - 1./8. * data[i-1][j];
		}
	}
	else {
		if (ascending) {
			return 3./8. * data[i][j] + 3./4. * data[i][j+1] - 1./8. * data[i][j+2];
		}
		else {
			return 3./8. * data[i][j+1] + 3./4. * data[i][j] - 1./8. * data[i][j-1];
		}
	}
}

// order 4 interpolation in 1D --- for points lying between two coarse nodes
// centre points
template <typename T, unsigned N>
Vector<T,N> order4interpolation1DCentre(const std::vector<std::vector<Vector<T,N>>>& data, int i, 
		int j, bool vertical)
{
	if (!vertical) {
		return 9./16. * (data[i][j] + data[i+1][j]) - 1./16. * (data[i-1][j] + data[i+2][j]);
	}
	else {
		return 9./16. * (data[i][j] + data[i][j+1]) - 1./16. * (data[i][j-1] + data[i][j+2]);
	}
}

// edge points
template <typename T, unsigned N>
Vector<T,N> order4interpolation1DEdge(const std::vector<std::vector<Vector<T,N>>>& data, int i, 
		int j, bool vertical, bool ascending)
{
	std::vector<T> weight {5./16., 15./16., -5./16., 1./16.};
	Vector<T,N> output {0.};
	if (!vertical) {
		if (ascending) {
			for (int k = 0; k < 4; ++k) {
				output += weight[k] * data[i+k][j];
			}
		}
		else {
			for (int k = 0; k < 4; ++k) {
				output += weight[k] * data[i+1-k][j];
			}
		}
	}
	else {
		if (ascending) {
			for (int k = 0; k < 4; ++k) {
				output += weight[k] * data[i][j+k];
			}
		}
		else {
			for (int k = 0; k < 4; ++k) {
				output += weight[k] * data[i][j+1-k];
			}
		}
	}
	return output;
}

// only for periodic boundary in z-direction
// this is for z-edge points 
template <typename T, unsigned N>
Vector<T,N> order4interpolation1DEdgePeriodic(const std::vector<std::vector<Vector<T,N>>>& data, int i, bool bottom)
{
	int zLength = data[0].size();
	if (bottom) {
		return 9./16. * (data[i][0] + data[i][1]) - 1./16. * (data[i][2] + data[i][zLength-1]);
	}
	else {
		return 9./16. * (data[i][zLength-1] + data[i][zLength-2]) - 1./16. * (data[i][zLength-3] + data[i][0]);
	}
}
// order 4 interpolation in 2D --- for points lying between four coarse nodes
// centre points
template <typename T, unsigned N>
Vector<T,N> order4interpolation2DCentre(const std::vector<std::vector<Vector<T,N>>>& data, 
		int i, int j)
{
	std::vector<T> weight {-1./16., 9./16., 9./16., -1./16.};
	Vector<T,N> output {0.};
	for (int m = 0; m < 4; ++m) {
		for (int n = 0; n < 4; ++n) {
			output += weight[m] * weight[n] * data[i+m-1][j+n-1];
		}
	}
	return output;
}

// edge points
// has to separate implementations according to positions
// vertical and ascending --- coupled to form position control parameters
// check notebook for more information
template <typename T, unsigned N>
Vector<T,N> order4interpolation2DEdge(const std::vector<std::vector<Vector<T,N>>>& data,
		int i, int j, bool vertical, bool ascending)
{
	std::vector<T> weightEdge {5./16., 15./16., -5./16., 1./16.};
	std::vector<T> weightCentre {-1./16., 9./16., 9./16., -1./16.};
	Vector<T,N> output {0.};
	if (!vertical) {
		// bottom -- 0,0
		if (!ascending) {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weightCentre[m] * weightEdge[n] * data[i+m-1][j+n];
				}
			}
		}
		// top -- 0,1
		else {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weightCentre[m] * weightEdge[n] * data[i+m-1][j+1-n];
				}
			}
		}
	}
	else {
		// left -- 1,0
		if (!ascending) {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weightEdge[m] * weightCentre[n] * data[i+m][j+n-1];
				}
			}
		}
		// right -- 1,1
		else {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weightEdge[m] * weightCentre[n] * data[i+1-m][j+n-1];
				}
			}
		}
	}
	return output;
}

// only for periodic boundary in z-direction
// this is for z-edge_centre points 
template <typename T, unsigned N>
Vector<T,N> order4interpolation2DEdgePeriodic(const std::vector<std::vector<Vector<T,N>>>& data, int i, bool bottom)
{
	int zLength = data[0].size();
	std::vector<T> weight {-1./16., 9./16., 9./16., -1./16.};
	Vector<T,N> output {0.};
	if (bottom) {
		for (int m = 0; m < 4; ++m) {
			output += weight[m] * weight[0] * data[i+m-1][zLength-1];
			output += weight[m] * weight[1] * data[i+m-1][0];
			output += weight[m] * weight[2] * data[i+m-1][1];
			output += weight[m] * weight[3] * data[i+m-1][2];
		}
	}
	else if (!bottom) {
		for (int m = 0; m < 4; ++m) {
			output += weight[m] * weight[0] * data[i+m-1][0];
			output += weight[m] * weight[1] * data[i+m-1][zLength-1];
			output += weight[m] * weight[2] * data[i+m-1][zLength-2];
			output += weight[m] * weight[3] * data[i+m-1][zLength-3];
		}
	}
	return output;
}

// corner points
template <typename T, unsigned N>
Vector<T,N> order4interpolation2DCorner(const std::vector<std::vector<Vector<T,N>>>& data,
		int i, int j, bool vertical, bool ascending)
{
	std::vector<T> weight {5./16., 15./16., -5./16., 1./16.};
	Vector<T,N> output {0.};
	if (!vertical) {
		// left bottom corner -- 0,0
		if (!ascending) {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weight[m] * weight[n] * data[i+m][j+n];
				}
			}
		}
		// left top corner -- 0,1
		else {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weight[m] * weight[n] * data[i+m][j+1-n];
				}
			}
		}
	}
	else {
		// right bottom corner -- 1,0
		if (!ascending) {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weight[m] * weight[n] * data[i+1-m][j+n];
				}
			}
		}
		// right top corner -- 1,1
		else {
			for (int m = 0; m < 4; ++m) {
				for (int n = 0; n < 4; ++n) {
					output += weight[m] * weight[n] * data[i+1-m][j+1-n];
				}
			}
		}
	}
	return output;
}

// only for periodic boundary in z-direction
// this is for z-corner_centre points 
template <typename T, unsigned N>
Vector<T,N> order4interpolation2DCornerPeriodic(const std::vector<std::vector<Vector<T,N>>>& data, int i, bool bottom)
{
	int zLength = data[0].size();
	std::vector<T> weightz {-1./16., 9./16., 9./16., -1./16.};
	std::vector<T> weight  {5./16., 15./16., -5./16., 1./16.};
	Vector<T,N> output {0.};
	if (bottom) {
		// left bottom corner
		if (i == 0) {
			for (int m = 0; m < 4; ++m) {
				output += weight[m] * weightz[0] * data[m][zLength-1];
				output += weight[m] * weightz[1] * data[m][0];
				output += weight[m] * weightz[2] * data[m][1];
				output += weight[m] * weightz[3] * data[m][2];
			}
		}
		// right bottom corner
		else {
			for (int m = 0; m < 4; ++m) {
				output += weight[m] * weightz[0] * data[i+1-m][zLength-1];
				output += weight[m] * weightz[1] * data[i+1-m][0];
				output += weight[m] * weightz[2] * data[i+1-m][1];
				output += weight[m] * weightz[3] * data[i+1-m][2];
			}
		}
	}
	else if (!bottom) {
		// left top corner
		if (i == 0) {
			for (int m = 0; m < 4; ++m) {
				output += weight[m] * weightz[0] * data[m][0];
				output += weight[m] * weightz[1] * data[m][zLength-1];
				output += weight[m] * weightz[2] * data[m][zLength-2];
				output += weight[m] * weightz[3] * data[m][zLength-3];
			}
		}
		// right top corner
		else {
			for (int m = 0; m < 4; ++m) {
				output += weight[m] * weightz[0] * data[i+1-m][0];
				output += weight[m] * weightz[1] * data[i+1-m][zLength-1];
				output += weight[m] * weightz[2] * data[i+1-m][zLength-2];
				output += weight[m] * weightz[3] * data[i+1-m][zLength-3];
			}
		}
	}
	return output;
}

// interpolate between two fine timestep to update coarse node values 
// --- only 2nd order accuracy
template <typename T, typename DESCRIPTOR>
void FineCoupler3D<T,DESCRIPTOR>::interpolate()
{
	auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
			Cell<T,DESCRIPTOR> coarseCell;
			coarseLattice.get(this->getCoarseLatticeR(i,j), coarseCell);

			T rho{};
			T u[DESCRIPTOR::d] {};
			T fNeq[DESCRIPTOR::q] {};
			T pi[util::TensorVal<DESCRIPTOR>::n] {};
			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
//			coarseCell.computeRhoU(rho, u);
//			coarseCell.computeAllMomenta(rho, u, pi);
//			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
//				fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
//			}

			_c2f_rho[i][j] = order2interpolation(Vector<T,1>(rho), _c2f_rho[i][j]);
//			for ( int id = 0; id < DESCRIPTOR::d; ++id ) {
//				if ( util::nearZero(u[id], 1.0e-15) ) {
//					u[id] = 0;
//				}
//			}
			_c2f_u[i][j] = order2interpolation(Vector<T,DESCRIPTOR::d>(u), _c2f_u[i][j]);

//			T fNeq[DESCRIPTOR::q] {};
//			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
			lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);

			_c2f_fneq[i][j] = order2interpolation(Vector<T,DESCRIPTOR::q>(fNeq), 
					_c2f_fneq[i][j]);
		}
	}
}

// Coupling step --- from coarse to fine, including interpolation
template <typename T, typename DESCRIPTOR>
void FineCoupler3D<T,DESCRIPTOR>::couple()
{
//	const auto& coarseLattice = this->_coarse.getSuperLattice();
	auto& fineLattice = this->_fine.getSuperLattice();

	// coupling of two coincide points 
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
//			const auto& coarsePos = this->getCoarseLatticeR(i, j);
			const auto& finePos = this->getFineLatticeR(2*i, 2*j);

			// coupling of feq
//			T rho{};
//			T u[DESCRIPTOR::d] {};
//			T fEq[DESCRIPTOR::q] {};
//			Cell<T,DESCRIPTOR> coarseCell;
//			coarseLattice.get(coarsePos, coarseCell);
//			coarseCell.computeRhoU(rho, u);
//			const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

			// Wait, shouldn't we use _c2f_rho and _c2f_u?
			T fEq[DESCRIPTOR::q] {};
			const T uSqr = _c2f_u[i][j][0]*_c2f_u[i][j][0]
						 + _c2f_u[i][j][1]*_c2f_u[i][j][1]
						 + _c2f_u[i][j][2]*_c2f_u[i][j][2];
//			lbHelpers<T,DESCRIPTOR>::computeFeq(coarseCell, fEq);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
//				fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
				fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, _c2f_rho[i][j][0], _c2f_u[i][j].data, uSqr);
			}

			// couling of fneq and assign the value back
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = fEq[iPop] + this->getScalingFactor() * _c2f_fneq[i][j][iPop];
			}
			fineLattice.set(finePos, fineCell);
		}
	}

	// coupling of 1d centre points --- y aligned
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a-1; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
			// left edge
			if ( i == 0 ) {
				// interpolation 
				const auto rho	= order4interpolation1DEdge(_c2f_rho, i, j, false, true);
				const auto u	= order4interpolation1DEdge(_c2f_u, i, j, false, true);
				const auto fneq	= order4interpolation1DEdge(_c2f_fneq, i, j, false, true);
				// get position and calculate new populations 
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i+1, 2*j);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}

			// right edge
			else if ( i == this->_coarseSize_a-2) {
				// interpolation
				const auto rho	= order4interpolation1DEdge(_c2f_rho, i, j, false, false);
				const auto u	= order4interpolation1DEdge(_c2f_u, i, j, false, false);
				const auto fneq	= order4interpolation1DEdge(_c2f_fneq, i, j, false, false);
				// get position and calculate new populations 
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i+1, 2*j);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}

			// inner points
			else {
				// interpolation
				const auto rho	= order4interpolation1DCentre(_c2f_rho, i, j, false);
				const auto u	= order4interpolation1DCentre(_c2f_u, i, j, false);
				const auto fneq	= order4interpolation1DCentre(_c2f_fneq, i, j, false);
				// get position and calculate new populations 
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i+1, 2*j);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}
		}
	}

	// coupling of 1d inner points --- x aligned
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b-1; ++j) {
			// bottom edge
			if ( j == 0 ) {
				// interpolation
				const auto rho	= order4interpolation1DEdge(_c2f_rho, i, j, true, true);
				const auto u	= order4interpolation1DEdge(_c2f_u, i, j, true, true);
				const auto fneq	= order4interpolation1DEdge(_c2f_fneq, i, j, true, true);

				// use periodic interpolation instead
//				const auto rho	= order4interpolation1DEdgePeriodic(_c2f_rho, i, true);
//				const auto u	= order4interpolation1DEdgePeriodic(_c2f_u, i, true);
//				const auto fneq	= order4interpolation1DEdgePeriodic(_c2f_fneq, i, true);

				// get position and calculate new populations
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i, 2*j+1);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}

			// top edge
			else if ( j == this->_coarseSize_b-2 ) {
				// interpolation
				const auto rho	= order4interpolation1DEdge(_c2f_rho, i, j, true, false);
				const auto u	= order4interpolation1DEdge(_c2f_u, i, j, true, false);
				const auto fneq	= order4interpolation1DEdge(_c2f_fneq, i, j, true, false);

				// use periodic interpolation instead
//				const auto rho	= order4interpolation1DEdgePeriodic(_c2f_rho, i, false);
//				const auto u	= order4interpolation1DEdgePeriodic(_c2f_u, i, false);
//				const auto fneq	= order4interpolation1DEdgePeriodic(_c2f_fneq, i, false);

				// get position and calculate new populations
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i, 2*j+1);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}

			// inner points
			else {
				// interpolation
				const auto rho	= order4interpolation1DCentre(_c2f_rho, i, j, true);
				const auto u	= order4interpolation1DCentre(_c2f_u, i, j, true);
				const auto fneq	= order4interpolation1DCentre(_c2f_fneq, i, j, true);
				// get position and calculate new populations
				const T uSqr = u*u;
				const auto finePos = this->getFineLatticeR(2*i, 2*j+1);
				Cell<T,DESCRIPTOR> fineCell;
				fineLattice.get(finePos, fineCell);
				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
					fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
							u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
				}
				// update populations
				fineLattice.set(finePos, fineCell);
			}
		}
	}

	//// coupling of 2d inner points
	// inner points
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 1; i < this->_coarseSize_a-2; ++i) {
		for (int j = 1; j < this->_coarseSize_b-2; ++j) {
			// interpolation
			const auto rho	= order4interpolation2DCentre(_c2f_rho, i, j);
			const auto u	= order4interpolation2DCentre(_c2f_u, i, j);
			const auto fneq	= order4interpolation2DCentre(_c2f_fneq, i, j);
			// get position and calculate new populations
			const T uSqr = u*u;
			const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
						u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
			}
			// update populations
			fineLattice.set(finePos, fineCell);
		}
	}
	// bottom edge
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 1; i < this->_coarseSize_a-2; ++i) {
		int j {0};
		// interpolation
		const auto rho	= order4interpolation2DEdge(_c2f_rho, i, j, false, false);
		const auto u	= order4interpolation2DEdge(_c2f_u, i, j, false, false);
		const auto fneq	= order4interpolation2DEdge(_c2f_fneq, i, j, false, false);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DEdgePeriodic(_c2f_rho, i, true);
//		const auto u	= order4interpolation2DEdgePeriodic(_c2f_u, i, true);
//		const auto fneq	= order4interpolation2DEdgePeriodic(_c2f_fneq, i, true);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}
	// top edge
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 1; i < this->_coarseSize_a-2; ++i) {
		int j {this->_coarseSize_b-2};
		// interpolation
		const auto rho	= order4interpolation2DEdge(_c2f_rho, i, j, false, true);
		const auto u	= order4interpolation2DEdge(_c2f_u, i, j, false, true);
		const auto fneq	= order4interpolation2DEdge(_c2f_fneq, i, j, false, true);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DEdgePeriodic(_c2f_rho, i, false);
//		const auto u	= order4interpolation2DEdgePeriodic(_c2f_u, i, false);
//		const auto fneq	= order4interpolation2DEdgePeriodic(_c2f_fneq, i, false);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}
	// left edge
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int j = 1; j < this->_coarseSize_b-2; ++j) {
		int i {0};
		// interpolation
		const auto rho	= order4interpolation2DEdge(_c2f_rho, i, j, true, false);
		const auto u	= order4interpolation2DEdge(_c2f_u, i, j, true, false);
		const auto fneq	= order4interpolation2DEdge(_c2f_fneq, i, j, true, false);
		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}
	// right edge
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int j = 1; j < this->_coarseSize_b-2; ++j) {
		int i {this->_coarseSize_a-2};
		// interpolation
		const auto rho	= order4interpolation2DEdge(_c2f_rho, i, j, true, true);
		const auto u	= order4interpolation2DEdge(_c2f_u, i, j, true, true);
		const auto fneq	= order4interpolation2DEdge(_c2f_fneq, i, j, true, true);
		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}

	// left bottom corner
	{
		int i {0};
		int j {0};
		// interpolation
		const auto rho	= order4interpolation2DCorner(_c2f_rho, i, j, false, false);
		const auto u	= order4interpolation2DCorner(_c2f_u, i, j, false, false);
		const auto fneq	= order4interpolation2DCorner(_c2f_fneq, i, j, false, false);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DCornerPeriodic(_c2f_rho, i, true);
//		const auto u	= order4interpolation2DCornerPeriodic(_c2f_u, i, true);
//		const auto fneq	= order4interpolation2DCornerPeriodic(_c2f_fneq, i, true);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}

	// right bottom corner
	{
		int i {this->_coarseSize_a-2};
		int j {0};
		// interpolation
		const auto rho	= order4interpolation2DCorner(_c2f_rho, i, j, true, false);
		const auto u	= order4interpolation2DCorner(_c2f_u, i, j, true, false);
		const auto fneq	= order4interpolation2DCorner(_c2f_fneq, i, j, true, false);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DCornerPeriodic(_c2f_rho, i, true);
//		const auto u	= order4interpolation2DCornerPeriodic(_c2f_u, i, true);
//		const auto fneq	= order4interpolation2DCornerPeriodic(_c2f_fneq, i, true);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}

	// left top corner
	{
		int i {0};
		int j {this->_coarseSize_b-2};
		// interpolation
		const auto rho	= order4interpolation2DCorner(_c2f_rho, i, j, false, true);
		const auto u	= order4interpolation2DCorner(_c2f_u, i, j, false, true);
		const auto fneq	= order4interpolation2DCorner(_c2f_fneq, i, j, false, true);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DCornerPeriodic(_c2f_rho, i, false);
//		const auto u	= order4interpolation2DCornerPeriodic(_c2f_u, i, false);
//		const auto fneq	= order4interpolation2DCornerPeriodic(_c2f_fneq, i, false);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}

	// right top corner
	{
		int i {this->_coarseSize_a-2};
		int j {this->_coarseSize_b-2};
		// interpolation
		const auto rho	= order4interpolation2DCorner(_c2f_rho, i, j, true, true);
		const auto u	= order4interpolation2DCorner(_c2f_u, i, j, true, true);
		const auto fneq	= order4interpolation2DCorner(_c2f_fneq, i, j, true, true);

		// use periodic interpolation
//		const auto rho	= order4interpolation2DCornerPeriodic(_c2f_rho, i, false);
//		const auto u	= order4interpolation2DCornerPeriodic(_c2f_u, i, false);
//		const auto fneq	= order4interpolation2DCornerPeriodic(_c2f_fneq, i, false);

		// get position and calculate new populations
		const T uSqr = u*u;
		const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
		Cell<T,DESCRIPTOR> fineCell;
		fineLattice.get(finePos, fineCell);
		for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
			fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
					u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
		}
		// update populations
		fineLattice.set(finePos, fineCell);
	}
}

// Overlap indicates the number of extra ghost fine lattice cells 
template <typename T, typename DESCRIPTOR>
BulkFineCoupler3D<T,DESCRIPTOR>::BulkFineCoupler3D(Grid3D<T,DESCRIPTOR>& coarse, 
		Grid3D<T,DESCRIPTOR>& fine, Vector<T,3> origin, Vector<T,3> extend, int overlap):
	Coupler3D<T,DESCRIPTOR>(coarse, fine, origin, extend, overlap),
	_c2f_rho(this->_coarseSize_a, std::vector<Vector<T,1>>(this->_coarseSize_b)),
	_c2f_u(this->_coarseSize_a, std::vector<Vector<T,DESCRIPTOR::d>>(this->_coarseSize_b)),
//				Vector<T,DESCRIPTOR::d>(T{}))),
	_c2f_fneq(this->_coarseSize_a, std::vector<Vector<T,DESCRIPTOR::q>>(this->_coarseSize_b)),
//				Vector<T,DESCRIPTOR::q>(T{})))
	_overlap(overlap)
{
	OstreamManager clout(std::cout, "C2F");

	const auto& coarseOrigin = this->getCoarseLatticeR(0,0);
	const auto& fineOrigin = this->getFineLatticeR(overlap,overlap);

	clout << "coarse origin: " << coarseOrigin[0] << " " << coarseOrigin[1] << " " 
		  << coarseOrigin[2] << " " << coarseOrigin[3] << std::endl;
	clout << "coarse size:	" << this->_coarseSize_a << " * " << this->_coarseSize_b
		  << std::endl;
	clout << "fine origin: " << fineOrigin[0] << " " << fineOrigin[1] << " "
		  << fineOrigin[2] << " " << fineOrigin[3] << std::endl;
	clout << "fine size:	" << this->_fineSize_a - overlap * 2 
		  << " * " << this->_fineSize_b - overlap * 2
		  << std::endl;
}

// Same as FineCoupler3D
template <typename T, typename DESCRIPTOR>
void BulkFineCoupler3D<T,DESCRIPTOR>::store()
{
	auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
//			if ( (i==0) || (j==0) || (i==this->_coarseSize_a - 1)
//					|| (j==this->_coarseSize_b - 1) ) {
			const auto pos = this->getCoarseLatticeR(i, j);
			T rho{};
			T u[DESCRIPTOR::d] {};
			T fNeq[DESCRIPTOR::q] {};
			T pi[util::TensorVal<DESCRIPTOR>::n] {};
			Cell<T,DESCRIPTOR> coarseCell;
			coarseLattice.get(pos, coarseCell);
			coarseCell.computeAllMomenta(rho, u, pi);
//			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
//			lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
			}

			_c2f_rho[i][j] = Vector<T,1>(rho);
			_c2f_u[i][j] = Vector<T,DESCRIPTOR::d>(u);
			_c2f_fneq[i][j] = Vector<T,DESCRIPTOR::q>(fNeq);
//			}
		}
	}
}

// Same as FineCoupler3D
template <typename T, typename DESCRIPTOR>
void BulkFineCoupler3D<T,DESCRIPTOR>::interpolate()
{
	auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
			Cell<T,DESCRIPTOR> coarseCell;
			coarseLattice.get(this->getCoarseLatticeR(i,j), coarseCell);

			T rho{};
			T u[DESCRIPTOR::d] {};
			T fNeq[DESCRIPTOR::q] {};
			T pi[util::TensorVal<DESCRIPTOR>::n] {};
//			lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
			coarseCell.computeAllMomenta(rho, u, pi);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
			}

			_c2f_rho[i][j] = order2interpolation(Vector<T,1>(rho), _c2f_rho[i][j]);
			_c2f_u[i][j] = order2interpolation(Vector<T,DESCRIPTOR::d>(u), _c2f_u[i][j]);

//			T fNeq[DESCRIPTOR::q] {};
//			lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);

			_c2f_fneq[i][j] = order2interpolation(Vector<T,DESCRIPTOR::q>(fNeq), 
					_c2f_fneq[i][j]);
		}
	}
}

// Coupling step --- no edge points now in this bulkCoupling class
template <typename T, typename DESCRIPTOR>
void BulkFineCoupler3D<T,DESCRIPTOR>::couple()
{
//	const auto& coarseLattice = this->_coarse.getSuperLattice();
	auto& fineLattice = this->_fine.getSuperLattice();
	const int coarseOverlap = this->_overlap / 2;

	// coupling of two coincide points 
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = coarseOverlap; i < this->_coarseSize_a - coarseOverlap; ++i) {
		for (int j = coarseOverlap; j < this->_coarseSize_b - coarseOverlap; ++j) {
//			const auto& coarsePos = this->getCoarseLatticeR(i, j);
			const auto& finePos = this->getFineLatticeR(2*i, 2*j);

			// coupling of feq
//			T fEq[DESCRIPTOR::q] {};
//			Cell<T,DESCRIPTOR> coarseCell;
//			coarseLattice.get(coarsePos, coarseCell);
//			lbHelpers<T,DESCRIPTOR>::computeFeq(coarseCell, fEq);
			T fEq[DESCRIPTOR::q] {};
			const T uSqr = _c2f_u[i][j][0]*_c2f_u[i][j][0]
						 + _c2f_u[i][j][1]*_c2f_u[i][j][1]
						 + _c2f_u[i][j][2]*_c2f_u[i][j][2];
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, _c2f_rho[i][j][0], _c2f_u[i][j].data, uSqr);
			}

			// couling of fneq and assign the value back
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = fEq[iPop] + this->getScalingFactor() * _c2f_fneq[i][j][iPop];
			}
			fineLattice.set(finePos, fineCell);
		}
	}

	// coupling of 1d centre points --- y aligned
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = coarseOverlap; i < this->_coarseSize_a-1-coarseOverlap; ++i) {
		for (int j = coarseOverlap; j < this->_coarseSize_b-coarseOverlap; ++j) {
			// interpolation
			const auto rho	= order4interpolation1DCentre(_c2f_rho, i, j, false);
			const auto u	= order4interpolation1DCentre(_c2f_u, i, j, false);
			const auto fneq	= order4interpolation1DCentre(_c2f_fneq, i, j, false);
			// get position and calculate new populations 
			const T uSqr = u*u;
			const auto finePos = this->getFineLatticeR(2*i+1, 2*j);
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
						u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
			}
			// update populations
			fineLattice.set(finePos, fineCell);
		}
	}

	// coupling of 1d inner points --- x aligned
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = coarseOverlap; i < this->_coarseSize_a-coarseOverlap; ++i) {
		for (int j = coarseOverlap; j < this->_coarseSize_b-1-coarseOverlap; ++j) {
			// interpolation
			const auto rho	= order4interpolation1DCentre(_c2f_rho, i, j, true);
			const auto u	= order4interpolation1DCentre(_c2f_u, i, j, true);
			const auto fneq	= order4interpolation1DCentre(_c2f_fneq, i, j, true);
			// get position and calculate new populations
			const T uSqr = u*u;
			const auto finePos = this->getFineLatticeR(2*i, 2*j+1);
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], 
						u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
			}
			// update populations
			fineLattice.set(finePos, fineCell);
		}
	}

	//// coupling of 2d inner points
	// inner points
#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = coarseOverlap; i < this->_coarseSize_a-1-coarseOverlap; ++i) {
		for (int j = coarseOverlap; j < this->_coarseSize_b-1-coarseOverlap; ++j) {
			// interpolation
			const auto rho	= order4interpolation2DCentre(_c2f_rho, i, j);
			const auto u	= order4interpolation2DCentre(_c2f_u, i, j);
			const auto fneq	= order4interpolation2DCentre(_c2f_fneq, i, j);
			// get position and calculate new populations
			const T uSqr = u*u;
			const auto finePos = this->getFineLatticeR(2*i+1, 2*j+1);
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0],
						u.data, uSqr) + this->getScalingFactor() * fneq[iPop];
			}
			// update populations
			fineLattice.set(finePos, fineCell);
		}
	}
}

template <typename T, typename DESCRIPTOR>
CoarseCoupler3D<T,DESCRIPTOR>::CoarseCoupler3D(Grid3D<T,DESCRIPTOR>& coarse,
		Grid3D<T,DESCRIPTOR>& fine, Vector<T,3> origin, Vector<T,3> extend) : 
	Coupler3D<T,DESCRIPTOR>(coarse, fine, origin, extend)
{
	OstreamManager clout(std::cout, "F2C");
	
	const auto& coarseOrigin	= this->getCoarseLatticeR(0, 0);
	const auto& fineOrigin		= this->getFineLatticeR(0, 0);

	clout << "coarse origin: " << coarseOrigin[0] << " " << coarseOrigin[1] << " "
		  << coarseOrigin[2] << " " << coarseOrigin[3] << std::endl;
	clout << "coarse size:	" << this->_coarseSize_a << " * " 
		  << this->_coarseSize_b << std::endl;
	clout << "fine origin: " << fineOrigin[0] << " " << fineOrigin[1] << " " 
		  << fineOrigin[2] << " " << fineOrigin[3] << std::endl;
	clout << "fine size:	" << this->_fineSize_a << " * "
		  << this->_fineSize_b << std::endl;
}

// this fneq is calculated by averaging all neighbourhood fneq of 
// the fine grids to get a restricted fneq value for the coarse grid
//
// Don't worry! There's overlap ghost points so it's okay!
// The coupling from fine to coarse nodes are done on at least i=3
// So it's okay
template <typename T, typename DESCRIPTOR>
void computeRestrictedFneq(const SuperLattice3D<T,DESCRIPTOR>& lattice,
						   Vector<int,4> latticeR,
						   T restrictedFneq[DESCRIPTOR::q])
{
	// sum all neighbourhood's fneq up respectively
	for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
		const auto neighbour = latticeR + Vector<int,4> {0, 
			descriptors::c<DESCRIPTOR>(iPop,0), descriptors::c<DESCRIPTOR>(iPop,1), 
			descriptors::c<DESCRIPTOR>(iPop,2)};
		Cell<T,DESCRIPTOR> cell;
		lattice.get(neighbour, cell);

		T fNeq[DESCRIPTOR::q] {};
		T rho {};
		T u[3] {};
//		T pi[util::TensorVal<DESCRIPTOR>::n] {};
		lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq);
//		cell.computeAllMomenta(rho, u, pi);

		for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
			restrictedFneq[jPop] += fNeq[jPop];
//			restrictedFneq[jPop] += firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(jPop, pi);
		}
	}

	// average 
	for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
		restrictedFneq[iPop] /= DESCRIPTOR::q;
	}
}

// couple from fine to coarse --- use restricted fNeq instead
template <typename T, typename DESCRIPTOR>
void CoarseCoupler3D<T,DESCRIPTOR>::couple() 
{
	const auto& fineLattice = this->_fine.getSuperLattice();
	auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
	#pragma omp parallel for
#endif
	for (int i = 0; i < this->_coarseSize_a; ++i) {
		for (int j = 0; j < this->_coarseSize_b; ++j) {
			const auto& finePos		= this->getFineLatticeR(2*i, 2*j);
			const auto& coarsePos	= this->getCoarseLatticeR(i, j);

			T rho {};
			T u[3] {};
			T fEq[DESCRIPTOR::q] {};
			T fNeq[DESCRIPTOR::q] {};
			T pi[util::TensorVal<DESCRIPTOR>::n] {};
			Cell<T,DESCRIPTOR> fineCell;
			fineLattice.get(finePos, fineCell);
			fineCell.computeRhoU(rho, u);
//			fineCell.computeAllMomenta(rho, u, pi);
			const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
//			const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
			for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
				fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
			}
//			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
//				fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
//			}

//			T fEq[DESCRIPTOR::q] {};
//			Cell<T,DESCRIPTOR> fineCell;
//			fineLattice.get(finePos, fineCell);
//			lbHelpers<T,DESCRIPTOR>::computeFeq(fineCell, fEq);

//			T fNeq[DESCRIPTOR::q] {};
			// Use restricted fneq
			computeRestrictedFneq(fineLattice, finePos, fNeq);
			// Alternative method for directly calculating fneq
//			T rho{};
//			T u[3] {};

//			if ( (i == 0) || (i == this->_coarseSize_a-1) 
//					|| (j == 0) || (j == this->_coarseSize_b-1) ) {
//			if ( (i == 0) || (i == this->_coarseSize_a-1) ) {
//				lbHelpers<T,DESCRIPTOR>::computeRhoU(fineCell, rho, u);
//				lbHelpers<T,DESCRIPTOR>::computeFneq(fineCell, fNeq, rho, u);
//				for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
//					fNeq[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
//				}
//			}
//			else {
//				computeRestrictedFneq(fineLattice, finePos, fNeq);
//			}
			
			Cell<T,DESCRIPTOR> coarseCell;
			coarseLattice.get(coarsePos, coarseCell);

			for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
				coarseCell[iPop] = fEq[iPop] + this->getInvScalingFactor() * fNeq[iPop];
			}

			coarseLattice.set(coarsePos, coarseCell);
		}
	}
}

}
#endif
