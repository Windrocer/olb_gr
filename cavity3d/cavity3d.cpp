/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>

// Parameters for the simulation setup
const int N			= 100;			// resolution of the model
const T Re			= 3200.;			// Reynolds number
const T maxPhysT	= 100.;		// max. simulation time in s, SI unit
const T D			= 1;			// Characteristic physical lenght
const T tau			= 0.501875;			// relaxation time 

const T lx			= 1.;
const T ly			= 1.;
const T lz			= 1.;


// Characteristics needed by Grid3D
const Characteristics<T> PhysCharacteristics(
		D,
		1.,
		1.*D/Re,
		1.0);

// Assign material numbers for boundary conditions
void prepareGeometry(Grid3D<T,DESCRIPTOR>& grid, 
					 Vector<T,3> origin, 
					 Vector<T,3> extend) {

	OstreamManager clout(std::cout, "prepareGeometry");
	clout << "Prepare Geometry ..." << std::endl;

	auto& converter		= grid.getConverter();
	auto& sGeometry		= grid.getSuperGeometry();
	const T deltaX		= converter.getPhysDeltaX();

	sGeometry.rename(0, 1);

	// x left wall
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4.,
									  origin[1] - deltaX*3./4.,
									  origin[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {deltaX*3./2., 
									  extend[1] + deltaX*3./2., 
									  extend[2] + deltaX*3./2.}; 

		IndicatorCuboid3D<T> leftWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, leftWall);
	}
	// x right wall
	{
		const Vector<T,3> wallOrigin {extend[0] - deltaX*3./4., 
									  origin[1] - deltaX*3./4.,
									  origin[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {deltaX*3./2., 
									  extend[1] + deltaX*3./2., 
									  extend[2] + deltaX*3./2.}; 

		IndicatorCuboid3D<T> rightWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, rightWall);
	}

	// y front wall
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4., 
									  origin[1] - deltaX*3./4.,
									  origin[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  deltaX*3./2.,
									  extend[2] + deltaX*3./2.};

		IndicatorCuboid3D<T> frontWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, frontWall);
	}
	// y back wall
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4., 
									  extend[1] - deltaX*3./4.,
									  origin[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  deltaX*3./2.,
									  extend[2] + deltaX*3./2.};

		IndicatorCuboid3D<T> backWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, backWall);
	}

	// z bottom wall
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4.,
									  origin[1] - deltaX*3./4.,
									  origin[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  extend[1] + deltaX*3./2.,
									  deltaX*3./2.};

		IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, bottomWall);
	}
	// z top boundary
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4.,
									  origin[1] - deltaX*3./4.,
									  extend[2] - deltaX*3./4.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  extend[1] + deltaX*3./2.,
									  deltaX*3./2.};

		IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 3, topWall);
	}

	// Clean
	sGeometry.clean();
	sGeometry.innerClean();
	sGeometry.checkForErrors();

	clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
					 Vector<T,3> domainOrigin, Vector<T,3> domainExtend) {
		
	OstreamManager clout(std::cout, "setupRefinement");
	clout << "Setup Refinement ..." << std::endl;

	// First level refinement 
	auto& boundFineGrid = coarseGrid.refine(domainOrigin, domainExtend,
			false, false, false, false);
	prepareGeometry(boundFineGrid, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> voidOrigin = {0.1, 0.1, 0.1};
		const Vector<T,3> voidExtend = {0.8, 0.8, 0.7};

		// add fine coupler --- from coarse to fine --- inner layer
		const Vector<T,3> extendX = {voidExtend[0], 0, 0};
		const Vector<T,3> extendY = {0, voidExtend[1], 0};
		const Vector<T,3> extendZ = {0, 0, voidExtend[2]};

		const Vector<T,3> extendXY = {voidExtend[0], voidExtend[1], 0};
		const Vector<T,3> extendYZ = {0, voidExtend[1], voidExtend[2]};
		const Vector<T,3> extendXZ = {voidExtend[0], 0, voidExtend[2]};

		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin, extendXY);
		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin, extendYZ);
		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin, extendXZ);

		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin + extendX, extendYZ);
		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin + extendY, extendXZ);
		coarseGrid.addFineCoupling(boundFineGrid, voidOrigin + extendZ, extendXY);

		// add coarse coupler --- from fine to coarse --- outer layer
		const Vector<T,3> outerOrigin = voidOrigin - coarseDeltaX;
		const Vector<T,3> outerExtend = voidExtend + 2*coarseDeltaX;

		const Vector<T,3> outerExtendX = {outerExtend[0], 0, 0};
		const Vector<T,3> outerExtendY = {0, outerExtend[1], 0};
		const Vector<T,3> outerExtendZ = {0, 0, outerExtend[2]};

		const Vector<T,3> outerExtendXY = {outerExtend[0], outerExtend[1], 0};
		const Vector<T,3> outerExtendYZ = {0, outerExtend[1], outerExtend[2]};
		const Vector<T,3> outerExtendXZ = {outerExtend[0], 0, outerExtend[2]};

		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin, outerExtendXY);
		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin, outerExtendYZ);
		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin, outerExtendXZ);

		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin + outerExtendX,
				outerExtendYZ);
		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin + outerExtendY,
				outerExtendXZ);
		coarseGrid.addCoarseCoupling(boundFineGrid, outerOrigin + outerExtendZ,
				outerExtendXY);

		// reset overlapped region on fineGrid
		const T fineDeltaX = boundFineGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> refinedFineOrigin = voidOrigin + fineDeltaX;
		const Vector<T,3> refinedFineExtend = voidExtend - 2*fineDeltaX;
		IndicatorCuboid3D<T> refinedFine(refinedFineExtend, refinedFineOrigin);
		boundFineGrid.getSuperGeometry().reset(refinedFine);

		// reset overlapped region on coarseGrid
		const Vector<T,3> refinedCoarseExtendXY = {domainExtend[0],
			domainExtend[1], outerOrigin[2] - coarseDeltaX};
		const Vector<T,3> refinedCoarseExtendYZ = {outerOrigin[0] - 
			coarseDeltaX, domainExtend[1], domainExtend[2]};
		const Vector<T,3> refinedCoarseExtendXZ = {domainExtend[0],
			outerOrigin[1] - coarseDeltaX, domainExtend[2]};

		IndicatorCuboid3D<T> refinedCoarseXY(refinedCoarseExtendXY,
				domainOrigin);
		coarseGrid.getSuperGeometry().reset(refinedCoarseXY);
		IndicatorCuboid3D<T> refinedCoarseYZ(refinedCoarseExtendYZ,
				domainOrigin);
		coarseGrid.getSuperGeometry().reset(refinedCoarseYZ);
		IndicatorCuboid3D<T> refinedCoarseXZ(refinedCoarseExtendXZ,
				domainOrigin);
		coarseGrid.getSuperGeometry().reset(refinedCoarseXZ);

		const Vector<T,3> refinedCoarseOriginX = {outerOrigin[0] + 
			outerExtend[0] + coarseDeltaX, domainOrigin[1], domainOrigin[2]};
		const Vector<T,3> refinedCoarseOriginY = {domainOrigin[0],
			outerOrigin[1] + outerExtend[1] + coarseDeltaX, domainOrigin[2]};
		const Vector<T,3> refinedCoarseOriginZ = {domainOrigin[0],
			domainOrigin[1], outerOrigin[2] + outerExtend[2] + coarseDeltaX};

		const Vector<T,3> refinedCoarseExtendXY2 = {domainExtend[0],
			domainExtend[1], domainExtend[2] - refinedCoarseOriginZ[2]};
		const Vector<T,3> refinedCoarseExtendYZ2 = {domainExtend[0] -
			refinedCoarseOriginX[0], domainExtend[1], domainExtend[2]};
		const Vector<T,3> refinedCoarseExtendXZ2 = {domainExtend[0],
			domainExtend[1] - refinedCoarseOriginY[1], domainExtend[2]};

		IndicatorCuboid3D<T> refinedCoarseXY2(refinedCoarseExtendXY2,
				refinedCoarseOriginZ);
		coarseGrid.getSuperGeometry().reset(refinedCoarseXY2);
		IndicatorCuboid3D<T> refinedcoarseYZ2(refinedCoarseExtendYZ2,
				refinedCoarseOriginX);
		coarseGrid.getSuperGeometry().reset(refinedcoarseYZ2);
		IndicatorCuboid3D<T> refinedCoarseXZ2(refinedCoarseExtendXZ2,
				refinedCoarseOriginY);
		coarseGrid.getSuperGeometry().reset(refinedCoarseXZ2);
	}

	// Second level refinement 
	const Vector<T,3> fineOrigin2 = {0.0, 0.0, 0.9};
	const Vector<T,3> fineExtend2 = {1.0, 1.0, 0.1};
	auto& boundFineGrid2 = boundFineGrid.refine(fineOrigin2, fineExtend2,
			false, false, false, false);
	prepareGeometry(boundFineGrid2, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = boundFineGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> voidOrigin = {0.0, 0.0, 0.9};
		const Vector<T,3> voidExtend = {1.0, 1.0, 0.1};

		// add fine coupler --- from coarse to fine --- inner layer
//		const Vector<T,3> extendX = {voidExtend[0], 0, 0};
//		const Vector<T,3> extendY = {0, voidExtend[1], 0};
//		const Vector<T,3> extendZ = {0, 0, voidExtend[2]};

		const Vector<T,3> extendXY = {voidExtend[0], voidExtend[1], 0};
//		const Vector<T,3> extendYZ = {0, voidExtend[1], voidExtend[2]};
//		const Vector<T,3> extendXZ = {voidExtend[0], 0, voidExtend[2]};

		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin, extendXY);
//		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin, extendYZ);
//		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin, extendXZ);
//
//		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin + extendX, extendYZ);
//		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin + extendY, extendXZ);
//		boundFineGrid.addFineCoupling(boundFineGrid2, voidOrigin + extendZ, extendXY);

		// add coarse coupler --- from fine to coarse --- outer layer
		const Vector<T,3> outerOrigin = voidOrigin 
			+ Vector<T,3> {0., 0., coarseDeltaX};
//		const Vector<T,3> outerExtend = voidExtend + 2*coarseDeltaX;
//
//		const Vector<T,3> outerExtendX = {outerExtend[0], 0, 0};
//		const Vector<T,3> outerExtendY = {0, outerExtend[1], 0};
//		const Vector<T,3> outerExtendZ = {0, 0, outerExtend[2]};
//
//		const Vector<T,3> outerExtendXY = {outerExtend[0], outerExtend[1], 0};
//		const Vector<T,3> outerExtendYZ = {0, outerExtend[1], outerExtend[2]};
//		const Vector<T,3> outerExtendXZ = {outerExtend[0], 0, outerExtend[2]};

		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin, extendXY);
//		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin, outerExtendYZ);
//		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin, outerExtendXZ);
//
//		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin + outerExtendX,
//				outerExtendYZ);
//		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin + outerExtendY,
//				outerExtendXZ);
//		boundFineGrid.addCoarseCoupling(boundFineGrid2, outerOrigin + outerExtendZ,
//				outerExtendXY);

//		// reset overlapped region on fine grid
//		const T fineDeltaX = boundFineGrid2.getConverter().getPhysDeltaX();
//		const Vector<T,3> refinedFineOrigin = voidOrigin + fineDeltaX;
//		const Vector<T,3> refinedFineExtend = voidExtend - 2*fineDeltaX;
//		IndicatorCuboid3D<T> refinedFine(refinedFineExtend, refinedFineOrigin);
//		boundFineGrid2.getSuperGeometry().reset(refinedFine);

		// reset overlapped region on coarse grid
		const Vector<T,3> refinedOrigin = voidOrigin 
			+ Vector<T,3> {0., 0., 2*coarseDeltaX};
		const Vector<T,3> refinedExtend = voidExtend
			- Vector<T,3> {0., 0., 2*coarseDeltaX};
		IndicatorCuboid3D<T> refinedCoarse(refinedExtend, refinedOrigin);
		boundFineGrid.getSuperGeometry().reset(refinedCoarse);
//		const Vector<T,3> refinedCoarseExtendXY = {domainExtend[0],
//			domainExtend[1], outerOrigin[2] - coarseDeltaX};
//		const Vector<T,3> refinedCoarseExtendYZ = {outerOrigin[0] - 
//			coarseDeltaX, domainExtend[1], domainExtend[2]};
//		const Vector<T,3> refinedCoarseExtendXZ = {domainExtend[0],
//			outerOrigin[1] - coarseDeltaX, domainExtend[2]};
//
//		IndicatorCuboid3D<T> refinedCoarseXY(refinedCoarseExtendXY,
//				domainOrigin);
//		boundFineGrid.getSuperGeometry().reset(refinedCoarseXY);
//		IndicatorCuboid3D<T> refinedCoarseYZ(refinedCoarseExtendYZ,
//				domainOrigin);
//		boundFineGrid.getSuperGeometry().reset(refinedCoarseYZ);
//		IndicatorCuboid3D<T> refinedCoarseXZ(refinedCoarseExtendXZ,
//				domainOrigin);
//		boundFineGrid.getSuperGeometry().reset(refinedCoarseXZ);
//
//		const Vector<T,3> refinedCoarseOriginX = {outerOrigin[0] + 
//			outerExtend[0] + coarseDeltaX, domainOrigin[1], domainOrigin[2]};
//		const Vector<T,3> refinedCoarseOriginY = {domainOrigin[0],
//			outerOrigin[1] + outerExtend[1] + coarseDeltaX, domainOrigin[2]};
//		const Vector<T,3> refinedCoarseOriginZ = {domainOrigin[0],
//			domainOrigin[1], outerOrigin[2] + outerExtend[2] + coarseDeltaX};
//
//		const Vector<T,3> refinedCoarseExtendXY2 = {domainExtend[0],
//			domainExtend[1], domainExtend[2] - refinedCoarseOriginZ[2]};
//		const Vector<T,3> refinedCoarseExtendYZ2 = {domainExtend[0] -
//			refinedCoarseOriginX[0], domainExtend[1], domainExtend[2]};
//		const Vector<T,3> refinedCoarseExtendXZ2 = {domainExtend[0],
//			domainExtend[1] - refinedCoarseOriginY[1], domainExtend[2]};
//
//		IndicatorCuboid3D<T> refinedCoarseXY2(refinedCoarseExtendXY2,
//				refinedCoarseOriginZ);
//		boundFineGrid.getSuperGeometry().reset(refinedCoarseXY2);
//		IndicatorCuboid3D<T> refinedcoarseYZ2(refinedCoarseExtendYZ2,
//				refinedCoarseOriginX);
//		boundFineGrid.getSuperGeometry().reset(refinedcoarseYZ2);
//		IndicatorCuboid3D<T> refinedCoarseXZ2(refinedCoarseExtendXZ2,
//				refinedCoarseOriginY);
//		boundFineGrid.getSuperGeometry().reset(refinedCoarseXZ2);
	}

//	// Third level refinement 
//	auto& boundFineGrid3 = boundFineGrid2.refine(domainOrigin, domainExtend,
//			false, false, false, false);
//	prepareGeometry(boundFineGrid3, domainOrigin, domainExtend);
//	{
//		const T coarseDeltaX = boundFineGrid2.getConverter().getPhysDeltaX();
//		const Vector<T,3> voidOrigin = {0.05, 0.05, 0.05};
//		const Vector<T,3> voidExtend = {0.9, 0.9, 0.8};
//
//		// add fine coupler --- from coarse to fine --- inner layer
//		const Vector<T,3> extendX = {voidExtend[0], 0, 0};
//		const Vector<T,3> extendY = {0, voidExtend[1], 0};
//		const Vector<T,3> extendZ = {0, 0, voidExtend[2]};
//
//		const Vector<T,3> extendXY = {voidExtend[0], voidExtend[1], 0};
//		const Vector<T,3> extendYZ = {0, voidExtend[1], voidExtend[2]};
//		const Vector<T,3> extendXZ = {voidExtend[0], 0, voidExtend[2]};
//
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin, extendXY);
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin, extendYZ);
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin, extendXZ);
//
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin + extendX, extendYZ);
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin + extendY, extendXZ);
//		boundFineGrid2.addFineCoupling(boundFineGrid3, voidOrigin + extendZ, extendXY);
//
//		// add coarse coupler --- from fine to coarse --- outer layer
//		const Vector<T,3> outerOrigin = voidOrigin - coarseDeltaX;
//		const Vector<T,3> outerExtend = voidExtend + 2*coarseDeltaX;
//
//		const Vector<T,3> outerExtendX = {outerExtend[0], 0, 0};
//		const Vector<T,3> outerExtendY = {0, outerExtend[1], 0};
//		const Vector<T,3> outerExtendZ = {0, 0, outerExtend[2]};
//
//		const Vector<T,3> outerExtendXY = {outerExtend[0], outerExtend[1], 0};
//		const Vector<T,3> outerExtendYZ = {0, outerExtend[1], outerExtend[2]};
//		const Vector<T,3> outerExtendXZ = {outerExtend[0], 0, outerExtend[2]};
//
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin, outerExtendXY);
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin, outerExtendYZ);
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin, outerExtendXZ);
//
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin + outerExtendX,
//				outerExtendYZ);
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin + outerExtendY,
//				outerExtendXZ);
//		boundFineGrid2.addCoarseCoupling(boundFineGrid3, outerOrigin + outerExtendZ,
//				outerExtendXY);
//
//		// reset overlapped region on fine grid
//		const T fineDeltaX = boundFineGrid3.getConverter().getPhysDeltaX();
//		const Vector<T,3> refinedFineOrigin = voidOrigin + fineDeltaX;
//		const Vector<T,3> refinedFineExtend = voidExtend - 2*fineDeltaX;
//		IndicatorCuboid3D<T> refinedFine(refinedFineExtend, refinedFineOrigin);
//		boundFineGrid3.getSuperGeometry().reset(refinedFine);
//
//		// reset overlapped region on coarse grid
//		const Vector<T,3> refinedCoarseExtendXY = {domainExtend[0],
//			domainExtend[1], outerOrigin[2] - coarseDeltaX};
//		const Vector<T,3> refinedCoarseExtendYZ = {outerOrigin[0] - 
//			coarseDeltaX, domainExtend[1], domainExtend[2]};
//		const Vector<T,3> refinedCoarseExtendXZ = {domainExtend[0],
//			outerOrigin[1] - coarseDeltaX, domainExtend[2]};
//
//		IndicatorCuboid3D<T> refinedCoarseXY(refinedCoarseExtendXY,
//				domainOrigin);
//		boundFineGrid2.getSuperGeometry().reset(refinedCoarseXY);
//
//		IndicatorCuboid3D<T> refinedCoarseYZ(refinedCoarseExtendYZ,
//				domainOrigin);
//		boundFineGrid2.getSuperGeometry().reset(refinedCoarseYZ);
//
//		IndicatorCuboid3D<T> refinedCoarseXZ(refinedCoarseExtendXZ,
//				domainOrigin);
//		boundFineGrid2.getSuperGeometry().reset(refinedCoarseXZ);
//
//		const Vector<T,3> refinedCoarseOriginX = {outerOrigin[0] + 
//			outerExtend[0] + coarseDeltaX, domainOrigin[1], domainOrigin[2]};
//		const Vector<T,3> refinedCoarseOriginY = {domainOrigin[0],
//			outerOrigin[1] + outerExtend[1] + coarseDeltaX, domainOrigin[2]};
//		const Vector<T,3> refinedCoarseOriginZ = {domainOrigin[0],
//			domainOrigin[1], outerOrigin[2] + outerExtend[2] + coarseDeltaX};
//
//		const Vector<T,3> refinedCoarseExtendXY2 = {domainExtend[0],
//			domainExtend[1], domainExtend[2] - refinedCoarseOriginZ[2]};
//		const Vector<T,3> refinedCoarseExtendYZ2 = {domainExtend[0] -
//			refinedCoarseOriginX[0], domainExtend[1], domainExtend[2]};
//		const Vector<T,3> refinedCoarseExtendXZ2 = {domainExtend[0],
//			domainExtend[1] - refinedCoarseOriginY[1], domainExtend[2]};
//
//		IndicatorCuboid3D<T> refinedCoarseXY2(refinedCoarseExtendXY2,
//				refinedCoarseOriginZ);
//		boundFineGrid2.getSuperGeometry().reset(refinedCoarseXY2);
//		IndicatorCuboid3D<T> refinedcoarseYZ2(refinedCoarseExtendYZ2,
//				refinedCoarseOriginX);
//		boundFineGrid2.getSuperGeometry().reset(refinedcoarseYZ2);
//		IndicatorCuboid3D<T> refinedCoarseXZ2(refinedCoarseExtendXZ2,
//				refinedCoarseOriginY);
//		boundFineGrid2.getSuperGeometry().reset(refinedCoarseXZ2);
//	}

	clout << "Setup Refinement ... OK" << std::endl;
}

// Create lattice structures
void prepareLattice(Grid3D<T,DESCRIPTOR>& grid) {

	OstreamManager clout(std::cout, "prepareLattice");
	clout << "Prepare lattice ..." << std::endl;

	auto& converter = grid.getConverter();
	auto& sGeometry = grid.getSuperGeometry();
	auto& sLattice  = grid.getSuperLattice();
//	const T deltaX  = converter.getPhysDeltaX();
	const T omega	= converter.getLatticeRelaxationFrequency();

	// Initialize dynamics
	Dynamics<T,DESCRIPTOR>& bulkDynamics = grid.addDynamics(
			std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
				new BGKdynamics<T,DESCRIPTOR>(
					omega, instances::getBulkMomenta<T,DESCRIPTOR>())));

	// Initialize boundary condition types
	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc = 
		grid.getOnLatticeBoundaryCondition();
	createLocalBoundaryCondition3D<T,DESCRIPTOR>(bc);
	
	// Define dynamics
	sLattice.defineDynamics(sGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>());
	sLattice.defineDynamics(sGeometry, 2, &instances::getBounceBack<T,DESCRIPTOR>());
//	sLattice.defineDynamics(sGeometry, 2, &bulkDynamics);

	auto bulkIndicator = sGeometry.getMaterialIndicator({1, 3});
	sLattice.defineDynamics(bulkIndicator, &bulkDynamics);

	// Define boundary conditions
	bc.addVelocityBoundary(sGeometry, 3, omega);

	// Initial conditions
	AnalyticalConst3D<T,T> rhoF {1.};
	Vector<T,3> velocityV {0.};
	Vector<T,3> velocityBC {converter.getCharLatticeVelocity(), 0, 0};
	AnalyticalConst3D<T,T> uF(velocityV);
	AnalyticalConst3D<T,T> uBC(velocityBC);

	// Initialize all values of distribution functions to their local equilibrium
	sLattice.defineRhoU( sGeometry, 1, rhoF, uF );
	sLattice.iniEquilibrium( sGeometry, 1, rhoF, uF );
	sLattice.defineRhoU( sGeometry, 3, rhoF, uBC );
	sLattice.iniEquilibrium( sGeometry, 3, rhoF, uBC );

	// Make the lattice ready for simulation
	sLattice.initialize();
	
	clout << "Prepare Lattice ... OK" << std::endl;
}

//// Set boundary values for start-scale inlet velocity
//void setBoundaryValues(Grid3D<T,DESCRIPTOR>& grid, int iT) {
//
//	OstreamManager clout(std::cout, "setBoundaryValues");
//	
//	auto& converter	= grid.getConverter();
//	auto& sGeometry	= grid.getSuperGeometry();
//	auto& sLattice	= grid.getSuperLattice();
//	
//	// No of time steps for smooth start-up
////	int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
//	int iTmaxStart = 1000;
//	int iTupdate = 10;
//	
//	if (iT%iTupdate == 0 && iT <= iTmaxStart) {
//		// Smooth start curve, sinus
////		SinusStartScale<T,int> StartScale(iTmaxStart, T(1));
//		
//		// Smooth start curve, polynomial
//		PolynomialStartScale<T,int> StartScale(iTmaxStart, T( 1 ));
//		
//		// Creates and sets the Poiseuille inflow profile using functors
//		int iTvec[1] = {iT};
//		T frac[1] = {};
//		StartScale( frac,iTvec );
//		std::vector<T> maxVelocity(3, 0);
//		maxVelocity[0] = 1.5*frac[0]*converter.getCharLatticeVelocity();
//		
//		T distance2Wall = converter.getConversionFactorLength()/2.;
////		RectanglePoiseuille3D<T> poiseuilleU(sGeometry, 3, maxVelocity,
////				distance2Wall, distance2Wall, distance2Wall);
//		Rectangle1DPoiseuille3D<T> poiseuilleU( sGeometry, 3, maxVelocity, 
//				distance2Wall, 0 );
//		sLattice.defineU(sGeometry, 3, poiseuilleU);
//	
//		clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
//	}
//}

// Output to vtk files 
void getResults(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT) {

	auto& converter = grid.getConverter();
	auto& sLattice  = grid.getSuperLattice();
	auto& sGeometry = grid.getSuperGeometry();
	
	SuperVTMwriter3D<T> vtmWriter(prefix);
	SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
	SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
	SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, sGeometry);
	SuperLatticeKnudsen3D<T,DESCRIPTOR> knudsen(sLattice);
	SuperLatticeRefinementMetricKnudsen3D<T,DESCRIPTOR> quality(sLattice, converter);
	vtmWriter.addFunctor(geometry);
	vtmWriter.addFunctor(velocity);
	vtmWriter.addFunctor(pressure);
	vtmWriter.addFunctor(knudsen);
	vtmWriter.addFunctor(quality);
	
	if (iT==0) {
	  vtmWriter.createMasterFile();
	}
	
	vtmWriter.write(iT);
}

//// Measure physical forces and pressure
//void takeMeasurements(Grid3D<T,DESCRIPTOR>& grid, int iT, bool print=true) {
//	auto& sLattice  = grid.getSuperLattice();
//	auto& sGeometry = grid.getSuperGeometry();
//	auto& converter = grid.getConverter();
//
//	SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
//	AnalyticalFfromSuperF3D<T> intpolatePressure(pressure, true);
//	SuperLatticePhysDrag3D<T,DESCRIPTOR> dragF(sLattice, sGeometry, 5, converter);
//
//	const T point1[3] {cylinderCenterX - D/2., cylinderCenterY, lz/2.};
//	const T point2[3] {cylinderCenterX + D/2., cylinderCenterY, lz/2.};
//
//	T pressureInFrontOfCylinder, pressureBehindCylinder;
//	intpolatePressure(&pressureInFrontOfCylinder, point1);
//	intpolatePressure(&pressureBehindCylinder,    point2);
//	const T pressureDrop = pressureInFrontOfCylinder - pressureBehindCylinder;
//	
//	const int input[4] {};
//	T drag[dragF.getTargetDim()];
//	dragF(drag, input);
//
//	if (print) {
//		OstreamManager clout(std::cout, "measurement");
//		clout << "pressureDrop=" << pressureDrop
//		      << "; drag=" << drag[0]
//		      << "; lift=" << drag[1]
//		      << endl;
//	}
//}
//
//// Capture the pressure around the cylinder --- middle z
//void capturePressure(Grid3D<T,DESCRIPTOR>& coarseGrid, int iT) {
//	auto& sLattice = coarseGrid.getSuperLattice();
//	auto& converter = coarseGrid.getConverter();
//	const T dR = converter.getPhysDeltaX();
//	SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
//	AnalyticalFfromSuperF3D<T> interpolatePressure(pressure, true);
//	const T radiusCylinder = D/2;
//	
//	ofstream myfile;
//	std::string filename {"tmp/pressure_" + to_string(iT) + ".csv"};
//	myfile.open(filename, fstream::app);
//	
//	for(int i = 0; i < 180; ++i) {
//		const T degree = i * 2. * 3.14159 / 180.0;
//		const T point[3] {cylinderCenterX 
//						  + (radiusCylinder + 1.5 * dR) * std::cos(degree),
//						  cylinderCenterY 
//						  + (radiusCylinder + 1.5 * dR) * std::sin(degree),
//						  lz/2.};
//		T pressureAtPoint {};
//		interpolatePressure(&pressureAtPoint, point);
//		
//		myfile << pressureAtPoint << std::endl;
//	}
//
//   myfile.close();
//}

int main(int argc, char* argv[]) {
	// OLB initialization
	olbInit(&argc, &argv);
	singleton::directories().setOutputDir("./tmp/");
	OstreamManager clout(std::cout, "main");
	
	// An overall indicator for all domain parts
	Vector<T,3> origin {0.};
	Vector<T,3> extend {lx, ly, lz};
	IndicatorCuboid3D<T> coarseDomain(extend, origin);

	// Construct a background coarse grid
	Grid3D<T,DESCRIPTOR> coarseGrid(
			coarseDomain,
			RelaxationTime<T>(tau),
			N,
			PhysCharacteristics,
			false, false, false);

	const Vector<T,3> domainOrigin = 
		coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
	const Vector<T,3> domainExtend = 
		coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

	// Prepare geometry ... 
	prepareGeometry(coarseGrid, domainOrigin, domainExtend);

	setupRefinement(coarseGrid, domainOrigin, domainExtend);

//	prepareLattice(coarseGrid);
	coarseGrid.forEachGrid(prepareLattice);

	clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() << std::endl;
	clout << "starting simulation..." << endl;

	Timer<T> timer( 
			coarseGrid.getConverter().getLatticeTime(maxPhysT),
			coarseGrid.getSuperGeometry().getStatistics().getNvoxel() );
	timer.start();

	// Convergence tracer every physical second
	util::ValueTracer<T> converge(coarseGrid.getConverter().getLatticeTime(0.1), 1e-8);

	const int statIter	= 10;
	const int vtkIter	= 1000;
	const int checkIter	= 50000;
	const int recordIter	= 50000;
    
	// matrix for recording sum of velocity and pressure of the coarse grid
	T sectionUx[N+1][N+1] {0.};
	T sectionUy[N+1][N+1] {0.};
	T sectionUz[N+1][N+1] {0.};
	T sectionP[N+1][N+1] {0.};
	// matrix for recording sum of velocity and pressure of the first fine grid
	T fine1sectionUx[2*N+1][2*N+1] {0.};
	T fine1sectionUy[2*N+1][2*N+1] {0.};
	T fine1sectionUz[2*N+1][2*N+1] {0.};
	T fine1sectionP[2*N+1][2*N+1] {0.};
	// matrix for recording sum of velocity and pressure of the second fine grid
	const int f2lz = (int)(0.1*4*N + 0.5);
	T fine2sectionUx[4*N+1][f2lz+1] {0.};
	T fine2sectionUy[4*N+1][f2lz+1] {0.};
	T fine2sectionUz[4*N+1][f2lz+1] {0.};
	T fine2sectionP[4*N+1][f2lz+1] {0.};

//	for (int iT = 0; iT <= coarseGrid.getConverter().getLatticeTime(maxPhysT); ++iT) {
 	for (int iT = 0; iT <= recordIter; ++iT) {

		// Load last checkpoint at startup
		if (iT == 0) {
			coarseGrid.forEachGrid("cavity3D", [&](Grid3D<T,DESCRIPTOR>& grid, 
						const std::string& id)
					{grid.getSuperLattice().load(id+".checkpoint");
					 clout << "Checkpoint loaded." << std::endl;});
		}
 
		if ( iT % vtkIter == 0 ) {
			coarseGrid.forEachGrid("cavity3D", [&](Grid3D<T,DESCRIPTOR>& grid,
						const std::string& id)
					{getResults(grid, id, iT); });
		}

		// coarse grid
		// extract point data recursively and added to global array
		{
			auto& sLattice = coarseGrid.getSuperLattice();
			auto& converter = coarseGrid.getConverter();
			for (int i = 0; i <= N; ++i) {
				for (int j = 0; j <= N; ++j) {
					const int point[4] {0, i, N/2, j};
					SuperLatticePhysVelocity3D<T,DESCRIPTOR> u(sLattice, converter);
					SuperLatticePhysPressure3D<T,DESCRIPTOR> p(sLattice, converter);
					T currentU[3] {};
					T currentP {};
					u(&currentU[0], &point[0]);
					p(&currentP, &point[0]);
					// sum up and divide them in the end
					sectionUx[i][j] += currentU[0];
					sectionUy[i][j] += currentU[1];
					sectionUz[i][j] += currentU[2];
					sectionP[i][j] += currentP;
				}
			}
		}

		// output matrix to a csv file
		if (iT == recordIter) {
			ofstream fileUx, fileUy, fileUz, fileP;
			fileUx.open("ux.csv", fstream::app);
			fileUy.open("uy.csv", fstream::app);
			fileUz.open("uz.csv", fstream::app);
			fileP.open("p.csv", fstream::app);
			for (int i = 0; i <= N; ++i) {
				for (int j = 0; j <= N; ++j) {
					fileUx << sectionUx[i][j]/recordIter << ",";
					fileUy << sectionUy[i][j]/recordIter << ",";
					fileUz << sectionUz[i][j]/recordIter << ",";
					fileP << sectionP[i][j]/recordIter << ",";
				}
				fileUx << std::endl;
				fileUy << std::endl;
				fileUz << std::endl;
				fileP << std::endl;
			}
			fileUx.close();
			fileUy.close();
			fileUz.close();
			fileP.close();
		}

		// fine 1 grid
		// extract point data recursively and added to global array
		{
			// locate the fine 1 grid using the centre node
			const Vector<T,3> cavityCentre {0.5, 0.5, 0.5};
			Grid3D<T,DESCRIPTOR>& fine1Grid = coarseGrid.locate(cavityCentre);

			auto& sLattice = fine1Grid.getSuperLattice();
			auto& converter = fine1Grid.getConverter();
			for (int i = 0; i <= 2*N; ++i) {
				for (int j = 0; j <= 2*N; ++j) {
					const int point[4] {0, i, 2*N/2, j};
					SuperLatticePhysVelocity3D<T,DESCRIPTOR> u(sLattice, converter);
					SuperLatticePhysPressure3D<T,DESCRIPTOR> p(sLattice, converter);
					T currentU[3] {};
					T currentP {};
					u(&currentU[0], &point[0]);
					p(&currentP, &point[0]);
					// sum up and divide them in the end
					fine1sectionUx[i][j] += currentU[0];
					fine1sectionUy[i][j] += currentU[1];
					fine1sectionUz[i][j] += currentU[2];
					fine1sectionP[i][j] += currentP;
				}
			}
		}

		// output matrix to a csv file
		if (iT == recordIter) {
			ofstream fine1fileUx, fine1fileUy, fine1fileUz, fine1fileP;
			fine1fileUx.open("f1ux.csv", fstream::app);
			fine1fileUy.open("f1uy.csv", fstream::app);
			fine1fileUz.open("f1uz.csv", fstream::app);
			fine1fileP.open("f1p.csv", fstream::app);
			for (int i = 0; i <= 2*N; ++i) {
				for (int j = 0; j <= 2*N; ++j) {
					fine1fileUx << fine1sectionUx[i][j]/recordIter << ",";
					fine1fileUy << fine1sectionUy[i][j]/recordIter << ",";
					fine1fileUz << fine1sectionUz[i][j]/recordIter << ",";
					fine1fileP << fine1sectionP[i][j]/recordIter << ",";
				}
				fine1fileUx << std::endl;
				fine1fileUy << std::endl;
				fine1fileUz << std::endl;
				fine1fileP << std::endl;
			}
			fine1fileUx.close();
			fine1fileUy.close();
			fine1fileUz.close();
			fine1fileP.close();
		}

		// fine 2 grid
		// extract point data recursively and added to global array
		{
			// locate the fine 1 grid using the centre node
			const Vector<T,3> cavityCentre {0.5, 0.5, 0.95};
			Grid3D<T,DESCRIPTOR>& fine2Grid = coarseGrid.locate(cavityCentre);

			auto& sLattice = fine2Grid.getSuperLattice();
			auto& converter = fine2Grid.getConverter();
			for (int i = 0; i <= 4*N; ++i) {
				for (int j = 0; j <= f2lz; ++j) {
					const int point[4] {0, i, 4*N/2, j + 4*N - f2lz};
					SuperLatticePhysVelocity3D<T,DESCRIPTOR> u(sLattice, converter);
					SuperLatticePhysPressure3D<T,DESCRIPTOR> p(sLattice, converter);
					T currentU[3] {};
					T currentP {};
					u(&currentU[0], &point[0]);
					p(&currentP, &point[0]);
					// sum up and divide them in the end
					fine2sectionUx[i][j] += currentU[0];
					fine2sectionUy[i][j] += currentU[1];
					fine2sectionUz[i][j] += currentU[2];
					fine2sectionP[i][j] += currentP;
				}
			}
		}

		// output matrix to a csv file
		if (iT == recordIter) {
			ofstream fine2fileUx, fine2fileUy, fine2fileUz, fine2fileP;
			fine2fileUx.open("f2ux.csv", fstream::app);
			fine2fileUy.open("f2uy.csv", fstream::app);
			fine2fileUz.open("f2uz.csv", fstream::app);
			fine2fileP.open("f2p.csv", fstream::app);
			for (int i = 0; i <= 4*N; ++i) {
				for (int j = 0; j <= f2lz; ++j) {
					fine2fileUx << fine2sectionUx[i][j]/recordIter << ",";
					fine2fileUy << fine2sectionUy[i][j]/recordIter << ",";
					fine2fileUz << fine2sectionUz[i][j]/recordIter << ",";
					fine2fileP << fine2sectionP[i][j]/recordIter << ",";
				}
				fine2fileUx << std::endl;
				fine2fileUy << std::endl;
				fine2fileUz << std::endl;
				fine2fileP << std::endl;
			}
			fine2fileUx.close();
			fine2fileUy.close();
			fine2fileUz.close();
			fine2fileP.close();
		}

//		setBoundaryValues(coarseGrid, iT);

//		coarseSLattice.collideAndStream();
		coarseGrid.collideAndStream();

		// Convergence check
		converge.takeValue(
				coarseGrid.getSuperLattice().getStatistics().getAverageEnergy(), true);

		// Update state
		if (iT % statIter == 0) {
			timer.update(iT);
			timer.printStep();

//			takeMeasurements(cylinderGrid, iT);
		}

//		if ( iT % recordIter == 0 ) {
//			capturePressure(coarseGrid, iT);
//		}
 
		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0) ) {
			if (iT % (2 * checkIter) == 0) {
				coarseGrid.forEachGrid("cavity3D_even", [&](Grid3D<T,DESCRIPTOR>& grid,
							const std::string& id)
						{grid.getSuperLattice().save(id+".checkpoint");
						 clout << "Even checkpoint saved." << std::endl;});
			}
			else {
				coarseGrid.forEachGrid("cavity3D_odd", [&](Grid3D<T,DESCRIPTOR>& grid,
							const std::string& id)
						{grid.getSuperLattice().save(id+".checkpoint");
						 clout << "Odd checkpoint saved." << std::endl;});
			}
		}

		if ((iT > 100) && converge.hasConverged()) {
			clout << "Simulation converged." << endl;
			break;
		}
	}

	timer.stop();
	timer.printSummary();

	// Assign material numbers 
	return 0;
}
