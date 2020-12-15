/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerländer
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

#ifndef REFINEMENT_GRID_2D_HH
#define REFINEMENT_GRID_2D_HH

#include "grid2D.h"
#include "coupler2D.hh"

#include "utilities/vectorHelpers.h"
#include "boundary/superBoundaryCondition2D.h"
#include "boundary/superOffBoundaryCondition2D.h"
#include "communication/heuristicLoadBalancer.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
Grid2D<T,DESCRIPTOR>::Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF,
                             RelaxationTime<T> tau,
                             int resolution,
                             Characteristics<T> characteristics):
  _domainF(std::move(domainF)),
  _characteristics(characteristics),
  _converter(new UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>(
               resolution, // resolution: number of voxels per charPhysL
               tau,        // latticeRelaxationTime: relaxation time, has to be greater than 0.5!
               characteristics.length,     // charPhysLength: reference length of simulation geometry
               characteristics.velocity,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
               characteristics.viscosity,  // physViscosity: physical kinematic viscosity in __m^2 / s__
               characteristics.density)),  // physDensity: physical density in __kg / m^3__
  _cuboids(new CuboidGeometry2D<T>(
             *_domainF,
             _converter->getConversionFactorLength(),
#ifdef PARALLEL_MODE_MPI
             singleton::mpi().getSize()
#else
             1
#endif
           )),
  _balancer(new HeuristicLoadBalancer<T>(
              *_cuboids)),
  _geometry(new SuperGeometry2D<T>(
              *_cuboids,
              *_balancer,
              2)),
  _lattice(new SuperLattice2D<T,DESCRIPTOR>(
             *_geometry))
{
  _converter->print();
}

template <typename T, typename DESCRIPTOR>
Grid2D<T,DESCRIPTOR>::Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF,
                             LatticeVelocity<T> latticeVelocity,
                             int resolution,
                             Characteristics<T> characteristics):
  _domainF(std::move(domainF)),
  _characteristics(characteristics),
  _converter(new UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>(
               resolution,      // resolution: number of voxels per charPhysL
               latticeVelocity, // charLatticeVelocity
               characteristics.length,     // charPhysLength: reference length of simulation geometry
               characteristics.velocity,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
               characteristics.viscosity,  // physViscosity: physical kinematic viscosity in __m^2 / s__
               characteristics.density)),  // physDensity: physical density in __kg / m^3__
  _cuboids(new CuboidGeometry2D<T>(
             *_domainF,
             _converter->getConversionFactorLength(),
#ifdef PARALLEL_MODE_MPI
             singleton::mpi().getSize()
#else
             1
#endif
           )),
  _balancer(new HeuristicLoadBalancer<T>(
              *_cuboids)),
  _geometry(new SuperGeometry2D<T>(
              *_cuboids,
              *_balancer,
              2)),
  _lattice(new SuperLattice2D<T,DESCRIPTOR>(
             *_geometry))
{
  _converter->print();
}

template <typename T, typename DESCRIPTOR>
Grid2D<T,DESCRIPTOR>::Grid2D(
  FunctorPtr<IndicatorF2D<T>>&& domainF,
  RelaxationTime<T> tau,
  int resolution,
  int re):
  Grid2D(std::forward<decltype(domainF)>(domainF),
         tau,
         resolution,
         Characteristics<T>(re))
{ }

template <typename T, typename DESCRIPTOR>
Grid2D<T,DESCRIPTOR>::Grid2D(
  FunctorPtr<IndicatorF2D<T>>&& domainF,
  LatticeVelocity<T> latticeVelocity,
  int resolution,
  int re):
  Grid2D(std::forward<decltype(domainF)>(domainF),
         latticeVelocity,
         resolution,
         Characteristics<T>(re))
{ }

template <typename T, typename DESCRIPTOR>
Characteristics<T> Grid2D<T,DESCRIPTOR>::getCharacteristics() const
{
  return _characteristics;
}

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::getConverter()
{
  return *_converter;
}

template <typename T, typename DESCRIPTOR>
CuboidGeometry2D<T>& Grid2D<T,DESCRIPTOR>::getCuboidGeometry()
{
  return *_cuboids;
}

template <typename T, typename DESCRIPTOR>
LoadBalancer<T>& Grid2D<T,DESCRIPTOR>::getLoadBalancer()
{
  return *_balancer;
}

template <typename T, typename DESCRIPTOR>
SuperGeometry2D<T>& Grid2D<T,DESCRIPTOR>::getSuperGeometry()
{
  return *_geometry;
}

template <typename T, typename DESCRIPTOR>
SuperLattice2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::getSuperLattice()
{
  return *_lattice;
}

template <typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::addDynamics(
  std::unique_ptr<Dynamics<T,DESCRIPTOR>>&& dynamics)
{
  Dynamics<T,DESCRIPTOR>& ref = *dynamics;
  _dynamics.emplace_back(std::move(dynamics));
  return ref;
}

template <typename T, typename DESCRIPTOR>
Spongedynamics<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::addSpongeDynamics(
		std::unique_ptr<Spongedynamics<T,DESCRIPTOR>>&& spongedynamics)
{
	Spongedynamics<T,DESCRIPTOR>& ref = *spongedynamics;
	_spongedynamics.emplace_back(std::move(spongedynamics));
	return ref;
}

template <typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::getOnLatticeBoundaryCondition()
{
  _onLatticeBoundaryConditions.emplace_back(
    new sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>(getSuperLattice()));
  return *_onLatticeBoundaryConditions.back();
}

template <typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::getOffLatticeBoundaryCondition()
{
  _offLatticeBoundaryConditions.emplace_back(
    new sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>(getSuperLattice()));
  return *_offLatticeBoundaryConditions.back();
}

template <typename T, typename DESCRIPTOR>
void Grid2D<T,DESCRIPTOR>::collideAndStream()
{
  for ( auto& fineCoupler : _fineCouplers ) {
    fineCoupler->store();
  }

  this->getSuperLattice().collideAndStream();

  for ( auto& fineGrid : _fineGrids ) {
    fineGrid->collideAndStream();
  }

  for ( auto& fineCoupler : _fineCouplers ) {
    fineCoupler->interpolate();
    fineCoupler->couple();
  }

  for ( auto& fineGrid : _fineGrids ) {
    fineGrid->collideAndStream();
  }

  for ( auto& fineCoupler : _fineCouplers ) {
    fineCoupler->store();
    fineCoupler->couple();
  }

  for ( auto& coarseCoupler : _coarseCouplers ) {
    coarseCoupler->couple();
  }
}

template <typename T, typename DESCRIPTOR>
FineCoupler2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::addFineCoupling(
  Grid2D<T,DESCRIPTOR>& fineGrid, Vector<T,2> origin, Vector<T,2> extend)
{
  _fineCouplers.emplace_back(
    new FineCoupler2D<T,DESCRIPTOR>(
      *this, fineGrid, origin, extend));
  return *_fineCouplers.back();
}

template <typename T, typename DESCRIPTOR>
CoarseCoupler2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::addCoarseCoupling(
  Grid2D<T,DESCRIPTOR>& fineGrid, Vector<T,2> origin, Vector<T,2> extend)
{
  _coarseCouplers.emplace_back(
    new CoarseCoupler2D<T,DESCRIPTOR>(
      *this, fineGrid, origin, extend));
  return *_coarseCouplers.back();
}

template <typename T, typename DESCRIPTOR>
Vector<T,2> Grid2D<T,DESCRIPTOR>::alignOriginToGrid(Vector<T,2> physR) const
{
  Vector<int,3> latticeR{};
  _cuboids->getLatticeR(physR, latticeR);
  return _cuboids->getPhysR(latticeR.toStdVector());
}

template <typename T, typename DESCRIPTOR>
Vector<T,2> Grid2D<T,DESCRIPTOR>::alignExtendToGrid(Vector<T,2> extend) const
{
  const T deltaX = _converter->getPhysDeltaX();
  return {
    static_cast<int>(std::round(extend[0] / deltaX)) * deltaX,
    static_cast<int>(std::round(extend[1] / deltaX)) * deltaX
  };
}

template <typename T, typename DESCRIPTOR>
void Grid2D<T,DESCRIPTOR>::forEachGrid(std::function<void(Grid2D<T,DESCRIPTOR>&)>&& f)
{
  f(*this);
  for (auto& grid : _fineGrids) {
    grid->forEachGrid(std::forward<decltype(f)>(f));
  }
}

template <typename T, typename DESCRIPTOR>
void Grid2D<T,DESCRIPTOR>::forEachGrid(
  const std::string& id,
  std::function<void(Grid2D<T,DESCRIPTOR>&,const std::string&)>&& f)
{
  f(*this, id);
  for (std::size_t i = 0; i < _fineGrids.size(); ++i) {
    _fineGrids[i]->forEachGrid(id + "_" + std::to_string(i),
                               std::forward<decltype(f)>(f));

  }
}

template <typename T, typename DESCRIPTOR>
void Grid2D<T,DESCRIPTOR>::forEachGrid(std::function<void(Grid2D<T,DESCRIPTOR>&, int)>&& f, 
		int iT)
{
	f(*this, iT);
	for (auto& grid : _fineGrids) {
		grid->forEachGrid(std::forward<decltype(f)>(f), iT);
	}
}

template <typename T, typename DESCRIPTOR>
Grid2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::locate(Vector<T,2> pos)
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
std::size_t Grid2D<T,DESCRIPTOR>::getActiveVoxelN() const
{
  std::size_t n = _geometry->getStatistics().getNvoxel();
  for (const auto& grid : _fineGrids) {
    n += grid->getActiveVoxelN();
  }
  return n;
}

template <typename T, typename DESCRIPTOR>
RefiningGrid2D<T,DESCRIPTOR>& Grid2D<T,DESCRIPTOR>::refine(
  Vector<T,2> wantedOrigin, Vector<T,2> wantedExtend, bool addCouplers)
{
  if (addCouplers) {
    auto& fineGrid = refine(wantedOrigin, wantedExtend, false);

    const Vector<T,2> origin = fineGrid.getOrigin();
    const Vector<T,2> extend = fineGrid.getExtend();

    const Vector<T,2> extendX = {extend[0],0};
    const Vector<T,2> extendY = {0,extend[1]};

    addFineCoupling(fineGrid, origin,           extendY);
    addFineCoupling(fineGrid, origin + extendX, extendY);
    addFineCoupling(fineGrid, origin + extendY, extendX);
    addFineCoupling(fineGrid, origin,           extendX);

    const T coarseDeltaX = getConverter().getPhysDeltaX();
    const Vector<T,2> innerOrigin  = origin + coarseDeltaX;
    const Vector<T,2> innerExtendX = extendX - Vector<T,2> {2*coarseDeltaX,0};
    const Vector<T,2> innerExtendY = extendY - Vector<T,2> {0,2*coarseDeltaX};

    addCoarseCoupling(fineGrid, innerOrigin,                innerExtendY);
    addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendY);
    addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendX);
    addCoarseCoupling(fineGrid, innerOrigin,                innerExtendX);

    return fineGrid;
  }
  else {
    const Vector<T,2> origin = alignOriginToGrid(wantedOrigin);
    const Vector<T,2> extend = alignExtendToGrid(wantedExtend);

    _fineGrids.emplace_back(
      new RefiningGrid2D<T,DESCRIPTOR>(*this, origin, extend));
    RefiningGrid2D<T,DESCRIPTOR>& fineGrid = *_fineGrids.back();

    auto refinedOverlap = fineGrid.getRefinedOverlap();
    _geometry->reset(*refinedOverlap);

    return fineGrid;
  }
}


template <typename T, typename DESCRIPTOR>
RefiningGrid2D<T,DESCRIPTOR>::RefiningGrid2D(
  Grid2D<T,DESCRIPTOR>& parentGrid, Vector<T,2> origin, Vector<T,2> extend):
  Grid2D<T,DESCRIPTOR>(
    std::unique_ptr<IndicatorF2D<T>>(new IndicatorCuboid2D<T>(extend, origin)),
    RelaxationTime<T>(2*parentGrid.getConverter().getLatticeRelaxationTime() - 0.5),
    2*parentGrid.getConverter().getResolution(),
    parentGrid.getCharacteristics()),
  _origin(origin),
  _extend(extend),
  _parentGrid(parentGrid) { }

template <typename T, typename DESCRIPTOR>
Vector<T,2> RefiningGrid2D<T,DESCRIPTOR>::getOrigin() const
{
  return _origin;
}

template <typename T, typename DESCRIPTOR>
Vector<T,2> RefiningGrid2D<T,DESCRIPTOR>::getExtend() const
{
  return _extend;
}

template <typename T, typename DESCRIPTOR>
std::unique_ptr<IndicatorF2D<T>> RefiningGrid2D<T,DESCRIPTOR>::getRefinedOverlap() const
{
  const T coarseDeltaX = _parentGrid.getConverter().getPhysDeltaX();

  return std::unique_ptr<IndicatorF2D<T>>(
           new IndicatorCuboid2D<T>(_extend - 4*coarseDeltaX, _origin + 2*coarseDeltaX));
}

}

#endif
