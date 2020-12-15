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

#ifndef REFINEMENT_GRID_3D_H
#define REFINEMENT_GRID_3D_H

#include <memory>
#include <vector>
#include <string>
#include <functional>

#include "core/unitConverter.h"
#include "core/superLattice3D.h"
#include "geometry/superGeometry3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/superGeometry3D.h"
#include "communication/loadBalancer.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "utilities/functorPtr.h"
#include "utilities/namedType.h"

namespace olb {

// Declaration as they will used in Grid2D class (in this file)
template <typename T, typename DESCRIPTOR> class FineCoupler3D;
template <typename T, typename DESCRIPTOR> class BulkFineCoupler3D;
template <typename T, typename DESCRIPTOR> class CoarseCoupler3D;
template <typename T, typename DESCRIPTOR> class RefiningGrid3D;

template <typename T, typename DESCRIPTOR> class sOnLatticeBoundaryCondition3D;
template <typename T, typename DESCRIPTOR> class sOffLatticeBoundaryCondition3D;

template <typename T>
using RelaxationTime = utilities::NamedType<T, struct NamedRelaxationTime>;

template <typename T>
using LatticeVelocity = utilities::NamedType<T, struct NamedLatticeVelocity>;

// Flow characteristic parameters --- used to construct the converter
template <typename T>
struct Characteristics {
	Characteristics(T l, T u, T nu, T rho):
		length(l),
		velocity(u),
		viscosity(nu),
		density(rho) {}

	// definition based on Re --- characterized parameter
	Characteristics(T Re):
		Characteristics(1.0, 1.0, 1.0/Re, 1.0) {}

	const T length;
	const T velocity;
	const T viscosity;
	const T density;
};

template <typename T, typename DESCRIPTOR>
class Grid3D{
protected:
	FunctorPtr<IndicatorF3D<T>> _domainF;
	const Characteristics<T> _characteristics;
	Vector<bool, 3> _periodicityOn;

	// components to define a superlattice --- group to be Grid3D lattice structure
	std::unique_ptr<UnitConverter<T,DESCRIPTOR>>	_converter;
	std::unique_ptr<CuboidGeometry3D<T>>			_cuboids;
	std::unique_ptr<LoadBalancer<T>>				_balancer;
	std::unique_ptr<SuperGeometry3D<T>>				_geometry;
	std::unique_ptr<SuperLattice3D<T,DESCRIPTOR>>	_lattice;

	// list of dynamics and BCs of all coarse and fine grids
	std::vector<std::unique_ptr<Dynamics<T,DESCRIPTOR>>>	_dynamics;
	std::vector<std::unique_ptr<sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>>>
		_onLatticeBoundaryConditions;
	std::vector<std::unique_ptr<sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>>>
		_offLatticeBoundaryConditions;

	// list of derived fine grid and their couplers
	std::vector<std::unique_ptr<RefiningGrid3D<T,DESCRIPTOR>>>		_fineGrids;
	std::vector<std::unique_ptr<FineCoupler3D<T,DESCRIPTOR>>>		_fineCouplers;
	std::vector<std::unique_ptr<BulkFineCoupler3D<T,DESCRIPTOR>>>	_bulkFineCouplers;
	std::vector<std::unique_ptr<CoarseCoupler3D<T,DESCRIPTOR>>>		_coarseCouplers;

public:
	Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
		   RelaxationTime<T> tau,
		   int resolution,
		   Characteristics<T> characteristics);
	Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
		   LatticeVelocity<T> latticeVelocity,
		   int resolution,
		   Characteristics<T> characteristics);
	Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF, 
		   RelaxationTime<T> tau,
		   int resolution,
		   Characteristics<T> characteristics,
		   bool periodicityX,
		   bool periodicityY,
		   bool periodicityZ);

	// constructor based on Re
	Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF, 
		   RelaxationTime<T> tau, 
		   int resolution, 
		   T re);
	Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF, 
		   LatticeVelocity<T> uMax, 
		   int resolution, 
		   T re);

	// Get-functions
	Characteristics<T> getCharacteristics() const;
	UnitConverter<T,DESCRIPTOR>&	getConverter();
	CuboidGeometry3D<T>&			getCuboidGeometry();
	LoadBalancer<T>&				getLoadBalancer();
	SuperGeometry3D<T>&				getSuperGeometry();
	SuperLattice3D<T,DESCRIPTOR>&	getSuperLattice();

	Dynamics<T,DESCRIPTOR>& addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>&& dynamics);
	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& getOnLatticeBoundaryCondition();
	sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& getOffLatticeBoundaryCondition();

	// core functions
	void collideAndStream();
	
	FineCoupler3D<T,DESCRIPTOR>& addFineCoupling(
			Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend);
	BulkFineCoupler3D<T,DESCRIPTOR>& addBulkFineCoupling(
			Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend,
			int overlap);
	CoarseCoupler3D<T,DESCRIPTOR>& addCoarseCoupling(
			Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend);

	// functions to align physical position to exact node 
	Vector<T,3> alignOriginToGrid(Vector<T,3> physR) const;
	Vector<T,3> alignExtendToGrid(Vector<T,3> physR) const;

	// create refining grid based on origin and extend
//	RefiningGrid3D<T,DESCRIPTOR>& refine(Vector<T,3> origin, Vector<T,3> extend, 
//			bool addCouplers=true);
	RefiningGrid3D<T,DESCRIPTOR>& refine(Vector<T,3> origin, Vector<T,3> extend, 
			bool periodicityX, bool periodicityY, bool periodicityZ, 
			bool addCouplers=true);

	// process function f on every grid including the coarse one and fine ones
	void forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&)>&& f);
	void forEachGrid(const std::string& id, std::function<void(Grid3D<T,DESCRIPTOR>&, 
				const std::string&)>&& f);
	void forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&, int)>&& f, int iT);
	void forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&, IndicatorF3D<T>&,
				STLreader<T>&)>&& f, IndicatorF3D<T>& indicator, 
				STLreader<T>& stlReader);
	void forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&, STLreader<T>&)>&& f,
			STLreader<T>& stlReader);
	void forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&, STLreader<T>&,
				Timer<T>&, int)>&& f, STLreader<T>& stlReader, Timer<T>& timer, int iT);


	/// Returns the finest grid representing a physical position
	/**
	* Only works if pos is actually contained in a node of the refinement tree.
	**/
	Grid3D<T,DESCRIPTOR>& locate(Vector<T,3> pos);

	std::size_t getActiveVoxelN() const;
};

template <typename T, typename DESCRIPTOR>
class RefiningGrid3D : public Grid3D<T,DESCRIPTOR> {
private:
	const Vector<T,3> _origin;
	const Vector<T,3> _extend;
	Vector<bool,3> _periodicityOn;

	Grid3D<T,DESCRIPTOR>& _parentGrid;

public:
	RefiningGrid3D(Grid3D<T,DESCRIPTOR>& parentGrid, Vector<T,3> origin, Vector<T,3> extend);
	RefiningGrid3D(Grid3D<T,DESCRIPTOR>& parentGrid, Vector<T,3> origin, Vector<T,3> extend,
				   bool periodicityX, bool periodicityY, bool periodicityZ);

	Vector<T,3> getOrigin() const;
	Vector<T,3> getExtend() const;

	/// Indicates the subdomain of the coarse grid rendered moot by refinement
	std::unique_ptr<IndicatorF3D<T>> getRefinedOverlap() const;
};

}

#endif
