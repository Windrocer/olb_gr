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

#ifndef REFINEMENT_GRID_3D_HH
#define REFINEMENT_GRID_3D_HH

#include "utilities/vectorHelpers.h"
#include "boundary/superBoundaryCondition3D.h"
#include "boundary/superOffBoundaryCondition3D.h"
#include "communication/heuristicLoadBalancer.h"

#include "grid3D.h"
//#include "coupler3D.hh"

namespace olb {

// definition of Grid3D class based on relaxation time
template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>::Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
							 RelaxationTime<T> tau,
							 int resolution,
							 Characteristics<T> characteristics):
	_domainF(std::move(domainF)),
	_characteristics(characteristics),
	_periodicityOn(false),
	_converter(new UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>(
				resolution,
				tau,
				characteristics.length,
				characteristics.velocity,
				characteristics.viscosity,
				characteristics.density)),
	_cuboids(new CuboidGeometry3D<T>(
				*domainF,
				_converter->getConversionFactorLength(),
#ifdef PARALLEL_MODE_MPI
				singleton::mpi().getSize()
#else
				1
#endif
				)),
	_balancer(new HeuristicLoadBalancer<T>(
				*_cuboids)),
	_geometry(new SuperGeometry3D<T>(
				*_cuboids,
				*_balancer,
				2)),
	_lattice(new SuperLattice3D<T,DESCRIPTOR>(
				*_geometry))
{
	_converter->print();
}

// definition of Grid3D class based on lattice velocity
template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>::Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
							 LatticeVelocity<T> latticeVelocity,
							 int resolution,
							 Characteristics<T> characteristics):
	_domainF(std::move(domainF)),
	_characteristics(characteristics),
	_periodicityOn(false),
	_converter(new UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>(
				resolution,
				latticeVelocity,
				characteristics.length,
				characteristics.velocity,
				characteristics.viscosity,
				characteristics.density)),
	_cuboids(new CuboidGeometry3D<T>(
				*domainF,
				_converter->getConversionFactorLength(),
#ifdef PARALLEL_MODE_MPI
				singleton::mpi().getSize()
#else
				1
#endif
				)),
	_balancer(new HeuristicLoadBalancer<T>(
				*_cuboids)),
	_geometry(new SuperGeometry3D<T>(
				*_cuboids,
				*_balancer,
				2)),
	_lattice(new SuperLattice3D<T,DESCRIPTOR>(
				*_geometry))
{
	_converter->print();
}

template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>::Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
							 RelaxationTime<T> tau,
							 int resolution,
							 Characteristics<T> characteristics,
							 bool periodicityX,
							 bool periodicityY,
							 bool periodicityZ):
	_domainF(std::move(domainF)),
	_characteristics(characteristics),
	_periodicityOn(periodicityX, periodicityY, periodicityZ),
	_converter(new UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>(
				resolution,
				tau,
				characteristics.length,
				characteristics.velocity,
				characteristics.viscosity,
				characteristics.density)),
	_cuboids(new CuboidGeometry3D<T>(
				*domainF,
				_converter->getConversionFactorLength(),
				periodicityX,
				periodicityY,
				periodicityZ,
#ifdef PARALLEL_MODE_MPI
				singleton::mpi().getSize()
#else
				1
#endif
				)),
	_balancer(new HeuristicLoadBalancer<T>(
				*_cuboids)),
	_geometry(new SuperGeometry3D<T>(
				*_cuboids,
				*_balancer,
				2)),
	_lattice(new SuperLattice3D<T,DESCRIPTOR>(
				*_geometry))
{
	_converter->print();
}

template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>::Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
							 RelaxationTime<T> tau,
							 int resolution,
							 T re):
	Grid3D(std::forward<decltype(domainF)>(domainF),
		   tau,
		   resolution,
		   Characteristics<T>(re))
{ }

template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>::Grid3D(FunctorPtr<IndicatorF3D<T>>&& domainF,
							 LatticeVelocity<T> latticeVelocity,
							 int resolution,
							 T re):
	Grid3D(std::forward<decltype(domainF)>(domainF),
		   latticeVelocity,
		   resolution,
		   Characteristics<T>(re))
{ }

template <typename T, typename DESCRIPTOR>
Characteristics<T> Grid3D<T,DESCRIPTOR>::getCharacteristics() const
{
	return _characteristics;
}

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::getConverter()
{
	return *_converter;
}

template <typename T, typename DESCRIPTOR>
CuboidGeometry3D<T>& Grid3D<T,DESCRIPTOR>::getCuboidGeometry()
{
	return *_cuboids;
}

template <typename T, typename DESCRIPTOR>
LoadBalancer<T>& Grid3D<T,DESCRIPTOR>::getLoadBalancer()
{
	return *_balancer;
}

template <typename T, typename DESCRIPTOR>
SuperGeometry3D<T>& Grid3D<T,DESCRIPTOR>::getSuperGeometry()
{
	return *_geometry;
}

template <typename T, typename DESCRIPTOR>
SuperLattice3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::getSuperLattice()
{
	return *_lattice;
}

template <typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::addDynamics(
		std::unique_ptr<Dynamics<T,DESCRIPTOR>>&& dynamics)
{
	Dynamics<T,DESCRIPTOR>& ref = *dynamics;
	_dynamics.emplace_back(std::move(dynamics));
	return ref;
}

template <typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::getOnLatticeBoundaryCondition()
{
	_onLatticeBoundaryConditions.emplace_back(
			new sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>(getSuperLattice()));
	return *_onLatticeBoundaryConditions.back();
}

template <typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::getOffLatticeBoundaryCondition()
{
	_offLatticeBoundaryConditions.emplace_back(
			new sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>(getSuperLattice()));
	return *_offLatticeBoundaryConditions.back();
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::collideAndStream()
{
	for ( auto& fineCoupler : _fineCouplers ) {
		fineCoupler->store();
	}

	for ( auto& bulkFineCoupler : _bulkFineCouplers ) {
		bulkFineCoupler->store();
	}

	this->getSuperLattice().collideAndStream();

	for ( auto& fineGrid : _fineGrids ) {
		fineGrid->collideAndStream();
	}

	// interpolation and first coupling
	for ( auto& fineCoupler : _fineCouplers ) {
		fineCoupler->interpolate();
		fineCoupler->couple();
	}

	for ( auto& bulkFineCoupler : _bulkFineCouplers ) {
		bulkFineCoupler->interpolate();
		bulkFineCoupler->couple();
	}

	for ( auto& fineGrid : _fineGrids ) {
		fineGrid->getSuperLattice().communicate();
	}

	for ( auto& fineGrid : _fineGrids ) {
		fineGrid->collideAndStream();
	}

	// use new coarse data and second coupling
	for ( auto& fineCoupler : _fineCouplers ) {
		fineCoupler->store();
		fineCoupler->couple();
	}

	for ( auto& bulkFineCoupler : _bulkFineCouplers ) {
		bulkFineCoupler->store();
		bulkFineCoupler->couple();
	}

	for ( auto& fineGrid : _fineGrids ) {
		fineGrid->getSuperLattice().communicate();
	}

	for ( auto& coarseCoupler : _coarseCouplers ) {
		coarseCoupler->couple();
	}

	this->getSuperLattice().communicate();
}

template <typename T, typename DESCRIPTOR>
FineCoupler3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::addFineCoupling(
		Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend)
{
	_fineCouplers.emplace_back(
			new FineCoupler3D<T,DESCRIPTOR>(
				*this, fineGrid, origin, extend));
	return *_fineCouplers.back();
}

template <typename T, typename DESCRIPTOR>
BulkFineCoupler3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::addBulkFineCoupling(
		Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend,
		int overlap)
{
	_bulkFineCouplers.emplace_back(
			new BulkFineCoupler3D<T,DESCRIPTOR>(
				*this, fineGrid, origin, extend, overlap));
	return *_bulkFineCouplers.back();
}

template <typename T, typename DESCRIPTOR>
CoarseCoupler3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::addCoarseCoupling(
		Grid3D<T,DESCRIPTOR>& fineGrid, Vector<T,3> origin, Vector<T,3> extend)
{
	_coarseCouplers.emplace_back(
			new CoarseCoupler3D<T,DESCRIPTOR>(
				*this, fineGrid, origin, extend));
	return *_coarseCouplers.back();
}

template <typename T, typename DESCRIPTOR>
Vector<T,3> Grid3D<T,DESCRIPTOR>::alignOriginToGrid(Vector<T,3> physR) const
{
	Vector<int,4> latticeR{};
	_cuboids->getLatticeR(physR, latticeR);
//	T alignedPhysR[3] {};
//	_cuboids->getPhysR(alignedPhysR, latticeR); 
//	return alignedPhysR;
	return _cuboids->getPhysR(latticeR.toStdVector());
}

template <typename T, typename DESCRIPTOR>
Vector<T,3> Grid3D<T,DESCRIPTOR>::alignExtendToGrid(Vector<T,3> extend) const
{
	const T deltaX = _converter->getPhysDeltaX();
	return {
		static_cast<int>(std::round(extend[0] / deltaX)) * deltaX,
		static_cast<int>(std::round(extend[1] / deltaX)) * deltaX,
		static_cast<int>(std::round(extend[2] / deltaX)) * deltaX
	};
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&)>&& f)
{
	f(*this);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f));
	}
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(
		const std::string& id,
		std::function<void(Grid3D<T,DESCRIPTOR>&, const std::string&)>&& f)
{
	f(*this, id);
	for(std::size_t i = 0; i < _fineGrids.size(); ++i) {
		_fineGrids[i]->forEachGrid(id + "_" + std::to_string(i),
									std::forward<decltype(f)>(f));
	}
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(
		std::function<void(Grid3D<T,DESCRIPTOR>&, int)>&& f, int iT)
{
	f(*this, iT);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f), iT);
	}
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(
		std::function<void(Grid3D<T,DESCRIPTOR>&, IndicatorF3D<T>&, 
			STLreader<T>&)>&& f, IndicatorF3D<T>& indicator,
		STLreader<T>& stlReader)
{
	f(*this, indicator, stlReader);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f), indicator, stlReader);
	}
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(
		std::function<void(Grid3D<T,DESCRIPTOR>&, STLreader<T>&)>&& f,
		STLreader<T>& stlReader)
{
	f(*this, stlReader);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f), stlReader);
	}
}

template <typename T, typename DESCRIPTOR>
void Grid3D<T,DESCRIPTOR>::forEachGrid(
		std::function<void(Grid3D<T,DESCRIPTOR>&, STLreader<T>&, Timer<T>&, int)>&& f,
		STLreader<T>& stlReader, Timer<T>& timer, int iT)
{
	f(*this, stlReader, timer, iT);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f), stlReader, timer, iT);
	}
}

template <typename T, typename DESCRIPTOR>
Grid3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::locate(Vector<T,3> pos)
{
	int iC;
	for (auto& grid : _fineGrids) {
		if (grid->getCuboidGeometry().getC(pos, iC)) {
			return grid->locate(pos);
		}
	}
	return *this;
}

template <typename T, typename DESCRIPTOR>
std::size_t Grid3D<T,DESCRIPTOR>::getActiveVoxelN() const
{
	std::size_t n = _geometry->getStatistics().getNvoxel();
	for (const auto& grid : _fineGrids) {
		n += grid->getActiveVoxelN();
	}
	return n;
}

template <typename T, typename DESCRIPTOR>
RefiningGrid3D<T,DESCRIPTOR>& Grid3D<T,DESCRIPTOR>::refine(
		Vector<T,3> wantedOrigin, Vector<T,3> wantedExtend, 
		bool periodicityX, bool periodicityY, bool periodicityZ, 
		bool addCouplers)
{
	// addCouplers indicator the refined grid lay inside the stream
	// therefore use BulkFineCoupler instead
	if (addCouplers) {
		auto& fineGrid = refine(wantedOrigin, wantedExtend, 
				periodicityX, periodicityY, periodicityZ, false);

		const T coarseDeltaX = getConverter().getPhysDeltaX();

		// Coarse to fine --- choose bulk or with boundary
		const Vector<T,3> origin = fineGrid.getOrigin();
		const Vector<T,3> extend = fineGrid.getExtend();

		const Vector<T,3> extendX = {extend[0], 0, 0};
		const Vector<T,3> extendY = {0, extend[1], 0}; 
		const Vector<T,3> extendZ = {0, 0, extend[2]};

		const Vector<T,3> extendXY = {extend[0], extend[1], 0};
		const Vector<T,3> extendYZ = {0, extend[1], extend[2]};
		const Vector<T,3> extendXZ = {extend[0], 0, extend[2]};

		addFineCoupling(fineGrid, origin, extendXY);
		addFineCoupling(fineGrid, origin, extendYZ);
		addFineCoupling(fineGrid, origin, extendXZ);

		addFineCoupling(fineGrid, origin + extendZ, extendXY);
		addFineCoupling(fineGrid, origin + extendX, extendYZ);
		addFineCoupling(fineGrid, origin + extendY, extendXZ);

//		// Bulk fine coupler exclusive
//		const Vector<T,3> origin = fineGrid.getOrigin() - coarseDeltaX;
//		const Vector<T,3> extend = fineGrid.getExtend() + coarseDeltaX * 2;
//
//		const Vector<T,3> originX = {origin[0] + coarseDeltaX, 
//									 origin[1], origin[2]};
//		const Vector<T,3> originY = {origin[0], origin[1] + coarseDeltaX,
//									 origin[2]};
//		const Vector<T,3> originZ = {origin[0], origin[1],
//									 origin[2] + coarseDeltaX};
//
//		const Vector<T,3> extendX = {extend[0] - coarseDeltaX, 0, 0};
//		const Vector<T,3> extendY = {0, extend[1] - coarseDeltaX, 0}; 
//		const Vector<T,3> extendZ = {0, 0, extend[2] - coarseDeltaX};
//
//		const Vector<T,3> extendXY = {extend[0], extend[1], 0};
//		const Vector<T,3> extendYZ = {0, extend[1], extend[2]};
//		const Vector<T,3> extendXZ = {extend[0], 0, extend[2]};
//
//
//		addBulkFineCoupling(fineGrid, originZ, extendXY, 2);
//		addBulkFineCoupling(fineGrid, originX, extendYZ, 2);
//		addBulkFineCoupling(fineGrid, originY, extendXZ, 2);
//
//		addBulkFineCoupling(fineGrid, origin + extendZ, extendXY, 2);
//		addBulkFineCoupling(fineGrid, origin + extendX, extendYZ, 2);
//		addBulkFineCoupling(fineGrid, origin + extendY, extendXZ, 2);

		// Fine to coarse --- not related to interpolation 
		const Vector<T,3> innerOrigin = origin + 1*coarseDeltaX;
		const Vector<T,3> innerExtend = extend - 2*coarseDeltaX;

		const Vector<T,3> innerExtendX = {innerExtend[0], 0, 0};
		const Vector<T,3> innerExtendY = {0, innerExtend[1], 0};
		const Vector<T,3> innerExtendZ = {0, 0, innerExtend[2]};

		const Vector<T,3> innerExtendXY = {innerExtend[0], innerExtend[1], 0};
		const Vector<T,3> innerExtendYZ = {0, innerExtend[1], innerExtend[2]};
		const Vector<T,3> innerExtendXZ = {innerExtend[0], 0, innerExtend[2]};

		addCoarseCoupling(fineGrid, innerOrigin, innerExtendXY);
		addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);
		addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);

		addCoarseCoupling(fineGrid, innerOrigin + innerExtendZ, innerExtendXY);
		addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendYZ);
		addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendXZ);

		auto refinedOverlap = fineGrid.getRefinedOverlap();
		_geometry->reset(*refinedOverlap);

		return fineGrid;
	}
	else {
		Vector<T,3> origin = alignOriginToGrid(wantedOrigin);
		Vector<T,3> extend = alignExtendToGrid(wantedExtend);

		const T coarseDeltaX = getConverter().getPhysDeltaX();
		const Vector<bool,3> periodicity {periodicityX, periodicityY, periodicityZ};

		for (int i = 0; i < 3; ++i) {
			if (periodicity[i] == true) {
				origin[i] -= 0.5 * coarseDeltaX;
				extend[i] += 0.5 * coarseDeltaX;
			}
		}

//		_fineGrids.emplace_back(
//				new RefiningGrid3D<T,DESCRIPTOR>(*this, origin, extend));
		_fineGrids.emplace_back(
				new RefiningGrid3D<T,DESCRIPTOR>(*this, origin, extend, 
					periodicityX, periodicityY, periodicityZ));
		RefiningGrid3D<T,DESCRIPTOR>& fineGrid = *_fineGrids.back();
//		auto refinedOverlap = fineGrid.getRefinedOverlap();
//		_geometry->reset(*refinedOverlap);

		return fineGrid;
	}
}

template <typename T, typename DESCRIPTOR>
RefiningGrid3D<T,DESCRIPTOR>::RefiningGrid3D(
		Grid3D<T,DESCRIPTOR>& parentGrid, Vector<T,3> origin, Vector<T,3> extend) :
	Grid3D<T,DESCRIPTOR>(
			std::unique_ptr<IndicatorF3D<T>>(new IndicatorCuboid3D<T>(extend, origin)),
			RelaxationTime<T>(2*parentGrid.getConverter().getLatticeRelaxationTime() - 0.5),
			2*parentGrid.getConverter().getResolution(),
			parentGrid.getCharacteristics() ),
	_origin(origin),
	_extend(extend),
	_periodicityOn(false),
	_parentGrid(parentGrid) { }

template <typename T, typename DESCRIPTOR>
RefiningGrid3D<T,DESCRIPTOR>::RefiningGrid3D(
		Grid3D<T,DESCRIPTOR>& parentGrid, Vector<T,3> origin, Vector<T,3> extend,
		bool periodicityX, bool periodicityY, bool periodicityZ):
	Grid3D<T,DESCRIPTOR>(
			std::unique_ptr<IndicatorF3D<T>>(new IndicatorCuboid3D<T>(extend, origin)),
			RelaxationTime<T>(2*parentGrid.getConverter().getLatticeRelaxationTime() - 0.5),
			2*parentGrid.getConverter().getResolution(),
			parentGrid.getCharacteristics(), 
			periodicityX,
			periodicityY,
			periodicityZ),
	_origin(origin),
	_extend(extend),
	_periodicityOn(periodicityX, periodicityY, periodicityZ),
	_parentGrid(parentGrid) { }

template <typename T, typename DESCRIPTOR>
Vector<T,3> RefiningGrid3D<T,DESCRIPTOR>::getOrigin() const
{
	return _origin;
}

template <typename T, typename DESCRIPTOR>
Vector<T,3> RefiningGrid3D<T,DESCRIPTOR>::getExtend() const
{
	return _extend;
}

template <typename T, typename DESCRIPTOR>
std::unique_ptr<IndicatorF3D<T>> RefiningGrid3D<T,DESCRIPTOR>::getRefinedOverlap() const
{
	const T coarseDeltaX = _parentGrid.getConverter().getPhysDeltaX();

	return std::unique_ptr<IndicatorF3D<T>>(
			new IndicatorCuboid3D<T>(_extend - 4*coarseDeltaX, _origin + 2*coarseDeltaX));
}


}

#endif
