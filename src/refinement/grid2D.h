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

#ifndef REFINEMENT_GRID_2D_H
#define REFINEMENT_GRID_2D_H

#include <memory>
#include <vector>
#include <string>
#include <functional>

#include "core/unitConverter.h"
#include "core/superLattice2D.h"
#include "geometry/superGeometry2D.h"
#include "geometry/cuboidGeometry2D.h"
#include "geometry/superGeometry2D.h"
#include "communication/loadBalancer.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "utilities/functorPtr.h"
#include "utilities/namedType.h"

namespace olb {


template <typename T, typename DESCRIPTOR> class FineCoupler2D;
template <typename T, typename DESCRIPTOR> class CoarseCoupler2D;
template <typename T, typename DESCRIPTOR> class RefiningGrid2D;

template <typename T, typename DESCRIPTOR> class sOnLatticeBoundaryCondition2D;
template <typename T, typename DESCRIPTOR> class sOffLatticeBoundaryCondition2D;

template <typename T>
using RelaxationTime = utilities::NamedType<T,struct NamedRelaxationTime>;

template <typename T>
using LatticeVelocity = utilities::NamedType<T,struct NamedLatticeVelocity>;

template <typename T>
struct Characteristics {
  Characteristics(T l, T u, T nu, T rho):
    length(l),
    velocity(u),
    viscosity(nu),
    density(rho) { }

  Characteristics(int Re):
    Characteristics(1.0, 1.0, 1.0/Re, 1.0) { }

  const T length;
  const T velocity;
  const T viscosity;
  const T density;
};


template <typename T, typename DESCRIPTOR>
class Grid2D {
protected:
  FunctorPtr<IndicatorF2D<T>> _domainF;
  const Characteristics<T> _characteristics;

  std::unique_ptr<UnitConverter<T,DESCRIPTOR>>  _converter;
  std::unique_ptr<CuboidGeometry2D<T>>          _cuboids;
  std::unique_ptr<LoadBalancer<T>>              _balancer;
  std::unique_ptr<SuperGeometry2D<T>>           _geometry;
  std::unique_ptr<SuperLattice2D<T,DESCRIPTOR>> _lattice;

  std::vector<std::unique_ptr<Dynamics<T,DESCRIPTOR>>>                       _dynamics;
  std::vector<std::unique_ptr<Spongedynamics<T,DESCRIPTOR>>>				_spongedynamics;
  std::vector<std::unique_ptr<sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>>>  _onLatticeBoundaryConditions;
  std::vector<std::unique_ptr<sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>>> _offLatticeBoundaryConditions;

  std::vector<std::unique_ptr<RefiningGrid2D<T,DESCRIPTOR>>> _fineGrids;

  std::vector<std::unique_ptr<FineCoupler2D<T,DESCRIPTOR>>>   _fineCouplers;
  std::vector<std::unique_ptr<CoarseCoupler2D<T,DESCRIPTOR>>> _coarseCouplers;

public:
  Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF,
         RelaxationTime<T> tau,
         int resolution,
         Characteristics<T> characteristics);
  Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF,
         LatticeVelocity<T> latticeVelocity,
         int resolution,
         Characteristics<T> characteristics);

  Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF, RelaxationTime<T> tau,   int resolution, int re);
  Grid2D(FunctorPtr<IndicatorF2D<T>>&& domainF, LatticeVelocity<T> uMax, int resolution, int re);

  Characteristics<T> getCharacteristics() const;

  UnitConverter<T,DESCRIPTOR>&  getConverter();
  CuboidGeometry2D<T>&          getCuboidGeometry();
  LoadBalancer<T>&              getLoadBalancer();
  SuperGeometry2D<T>&           getSuperGeometry();
  SuperLattice2D<T,DESCRIPTOR>& getSuperLattice();

  Dynamics<T,DESCRIPTOR>& addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>&& dynamics);
  Spongedynamics<T,DESCRIPTOR>& addSpongeDynamics(
		  std::unique_ptr<Spongedynamics<T,DESCRIPTOR>>&& spongedynamics);
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>&  getOnLatticeBoundaryCondition();
  sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& getOffLatticeBoundaryCondition();

  void collideAndStream();

  FineCoupler2D<T,DESCRIPTOR>& addFineCoupling(
    Grid2D<T,DESCRIPTOR>& fineGrid, Vector<T,2> origin, Vector<T,2> extend);
  CoarseCoupler2D<T,DESCRIPTOR>& addCoarseCoupling(
    Grid2D<T,DESCRIPTOR>& fineGrid, Vector<T,2> origin, Vector<T,2> extend);

  Vector<T,2> alignOriginToGrid(Vector<T,2> physR) const;
  Vector<T,2> alignExtendToGrid(Vector<T,2> physR) const;

  RefiningGrid2D<T,DESCRIPTOR>& refine(Vector<T,2> origin, Vector<T,2> extend, bool addCouplers=true);

  void forEachGrid(std::function<void(Grid2D<T,DESCRIPTOR>&)>&& f);
  void forEachGrid(const std::string& id, std::function<void(Grid2D<T,DESCRIPTOR>&,const std::string&)>&& f);
  void forEachGrid(std::function<void(Grid2D<T,DESCRIPTOR>&, int)>&& f, int iT);

  /// Returns the finest grid representing a physical position
  /**
   * Only works if pos is actually contained in a node of the refinement tree.
   **/
  Grid2D<T,DESCRIPTOR>& locate(Vector<T,2> pos);

  std::size_t getActiveVoxelN() const;

};

template <typename T, typename DESCRIPTOR>
class RefiningGrid2D : public Grid2D<T,DESCRIPTOR> {
private:
  const Vector<T,2> _origin;
  const Vector<T,2> _extend;

  Grid2D<T,DESCRIPTOR>& _parentGrid;

public:
  RefiningGrid2D(Grid2D<T,DESCRIPTOR>& parentGrid, Vector<T,2> origin, Vector<T,2> extend);

  Vector<T,2> getOrigin() const;
  Vector<T,2> getExtend() const;

  /// Indicates the subdomain of the coarse grid rendered moot by refinement
  std::unique_ptr<IndicatorF2D<T>> getRefinedOverlap() const;

};


}

#endif
