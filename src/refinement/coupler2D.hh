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

#ifndef REFINEMENT_COUPLER_2D_HH
#define REFINEMENT_COUPLER_2D_HH

#include "coupler2D.h"

#include "dynamics/lbHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
T Coupler2D<T,DESCRIPTOR>::getScalingFactor() const
{
  const T coarseTau = _coarse.getConverter().getLatticeRelaxationTime();
  return (coarseTau - 0.25) / coarseTau;
}

template <typename T, typename DESCRIPTOR>
T Coupler2D<T,DESCRIPTOR>::getInvScalingFactor() const
{
  return 1./getScalingFactor();
}

template <typename T, typename DESCRIPTOR>
const Vector<int,3>& Coupler2D<T,DESCRIPTOR>::getFineLatticeR(int y) const
{
  return _fineLatticeR[y];
}

template <typename T, typename DESCRIPTOR>
const Vector<int,3>& Coupler2D<T,DESCRIPTOR>::getCoarseLatticeR(int y) const
{
  return _coarseLatticeR[y];
}

template <typename T, typename DESCRIPTOR>
Coupler2D<T,DESCRIPTOR>::Coupler2D(Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
                                   Vector<T,2> origin, Vector<T,2> extend):
  _coarse(coarse),
  _fine(fine),
  _coarseSize(floor(extend.norm() / coarse.getConverter().getPhysDeltaX() + 0.5)+1),
  _fineSize(2*_coarseSize-1),
  _vertical(util::nearZero(extend[0])),
  _physOrigin(_coarse.alignOriginToGrid(origin)),
  _coarseLatticeR(_coarseSize),
  _fineLatticeR(_fineSize)
{
  OLB_ASSERT(util::nearZero(extend[0]) || util::nearZero(extend[1]), "Coupling is only implemented alongside unit vectors");

  const auto& coarseGeometry = _coarse.getCuboidGeometry();
  const auto& fineGeometry   = _fine.getCuboidGeometry();

  const T deltaX = _fine.getConverter().getPhysDeltaX();
  const Vector<T,2> stepPhysR = _vertical ? Vector<T,2> {0, deltaX} :
                                Vector<T,2> {deltaX, 0};

  for (int i=0; i < _fineSize; ++i) {
    if (i % 2 == 0) {
      coarseGeometry.getLatticeR(_physOrigin + i*stepPhysR, _coarseLatticeR[i/2]);
    }

    fineGeometry.getLatticeR(_physOrigin + i*stepPhysR, _fineLatticeR[i]);
  }
}


template <typename T, typename DESCRIPTOR>
FineCoupler2D<T,DESCRIPTOR>::FineCoupler2D(Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
    Vector<T,2> origin, Vector<T,2> extend):
  Coupler2D<T,DESCRIPTOR>(coarse, fine, origin, extend),
  _c2f_rho(this->_coarseSize),
  _c2f_u(this->_coarseSize, Vector<T,2>(T{})),
  _c2f_fneq(this->_coarseSize, Vector<T,DESCRIPTOR::q>(T{}))
{
  OstreamManager clout(std::cout,"C2F");

  const auto& coarseOrigin = this->getCoarseLatticeR(0);
  const auto& fineOrigin   = this->getFineLatticeR(0);

  clout << "coarse origin: " << coarseOrigin[0] << " " << coarseOrigin[1] << " " << coarseOrigin[2] << std::endl;
  clout << "fine origin:   " << fineOrigin[0]   << " " << fineOrigin[1]   << " " << fineOrigin[2]   << std::endl;
  clout << "fine size:     " << this->_fineSize << std::endl;
}

template <typename T, typename DESCRIPTOR>
void FineCoupler2D<T,DESCRIPTOR>::store()
{
  auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (int y=0; y < this->_coarseSize; ++y) {
    const auto pos = this->getCoarseLatticeR(y);
    T rho{};
    T u[2] {};
    T fNeq[DESCRIPTOR::q] {};
    Cell<T,DESCRIPTOR> coarseCell;
    coarseLattice.get(pos, coarseCell);
    lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);
    lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);

    _c2f_rho[y]  = Vector<T,1>(rho);
    _c2f_u[y]    = Vector<T,2>(u);
    _c2f_fneq[y] = Vector<T,DESCRIPTOR::q>(fNeq);
  }
}

template <typename T, unsigned N>
Vector<T,N> order2interpolation(const Vector<T,N>& f0, const Vector<T,N>& f1)
{
  return 0.5 * (f0 + f1);
}

template <typename T, unsigned N>
Vector<T,N> order2interpolation(const std::vector<Vector<T,N>>& data, int y)
{
  return 0.5 * (data[y] + data[y+1]);
}

template <typename T, unsigned N>
Vector<T,N> order3interpolation(const std::vector<Vector<T,N>>& data, int y, bool ascending)
{
  if (ascending) {
    return 3./8. * data[y] + 3./4. * data[y+1] - 1./8. * data[y+2];
  }
  else {
    return 3./8. * data[y] + 3./4. * data[y-1] - 1./8. * data[y-2];
  }
}

template <typename T, unsigned N>
Vector<T,N> order4interpolation(const std::vector<Vector<T,N>>& data, int y)
{
  return 9./16. * (data[y] + data[y+1]) - 1./16. * (data[y-1] + data[y+2]);
}


template <typename T, typename DESCRIPTOR>
void FineCoupler2D<T,DESCRIPTOR>::interpolate()
{
  auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (int y=0; y < this->_coarseSize; ++y) {
    Cell<T,DESCRIPTOR> coarseCell;
    coarseLattice.get(this->getCoarseLatticeR(y), coarseCell);

    T rho{};
    T u[2] {};
    lbHelpers<T,DESCRIPTOR>::computeRhoU(coarseCell, rho, u);

    _c2f_rho[y] = order2interpolation(Vector<T,1>(rho), _c2f_rho[y]);
    _c2f_u[y]   = order2interpolation(Vector<T,2>(u), _c2f_u[y]);

    T fNeq[DESCRIPTOR::q] {};
    lbHelpers<T,DESCRIPTOR>::computeFneq(coarseCell, fNeq, rho, u);

    _c2f_fneq[y] = order2interpolation(Vector<T,DESCRIPTOR::q>(fNeq), _c2f_fneq[y]);
  }
}

template <typename T, typename DESCRIPTOR>
void FineCoupler2D<T,DESCRIPTOR>::couple()
{
  const auto& coarseLattice = this->_coarse.getSuperLattice();
  auto& fineLattice   = this->_fine.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (int y=0; y < this->_coarseSize; ++y) {
    const auto& coarsePos = this->getCoarseLatticeR(y);
    const auto& finePos   = this->getFineLatticeR(2*y);

    T fEq[DESCRIPTOR::q] {};
    Cell<T,DESCRIPTOR> coarseCell;
    coarseLattice.get(coarsePos, coarseCell);
    lbHelpers<T,DESCRIPTOR>::computeFeq(coarseCell, fEq);

    Cell<T,DESCRIPTOR> cell;
    fineLattice.get(finePos, cell);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop] + this->getScalingFactor() * _c2f_fneq[y][iPop];
    }
    fineLattice.set(finePos, cell);
  }

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (int y=1; y < this->_coarseSize-2; ++y) {
    const auto rho  = order4interpolation(_c2f_rho,  y);
    const auto u    = order4interpolation(_c2f_u,    y);
    const auto fneq = order4interpolation(_c2f_fneq, y);

    const T uSqr = u*u;

    const auto finePos = this->getFineLatticeR(1+2*y);
    Cell<T,DESCRIPTOR> fineCell;
    fineLattice.get(finePos, fineCell);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], u.data, uSqr)
                       + this->getScalingFactor() * fneq[iPop];
    }

    fineLattice.set(finePos, fineCell);
  }

  {
    const auto rho  = order3interpolation(_c2f_rho,  0, true);
    const auto u    = order3interpolation(_c2f_u,    0, true);
    const auto fneq = order3interpolation(_c2f_fneq, 0, true);

    const T uSqr = u*u;

    const auto& finePos = this->getFineLatticeR(1);
    Cell<T,DESCRIPTOR> fineCell;
    fineLattice.get(finePos, fineCell);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], u.data, uSqr)
                       + this->getScalingFactor() * fneq[iPop];
    }

    fineLattice.set(finePos, fineCell);
  }

  {
    const auto rho  = order3interpolation(_c2f_rho,  this->_coarseSize-1, false);
    const auto u    = order3interpolation(_c2f_u,    this->_coarseSize-1, false);
    const auto fneq = order3interpolation(_c2f_fneq, this->_coarseSize-1, false);

    const T uSqr = u*u;

    const auto& finePos = this->getFineLatticeR(this->_fineSize-2);
    Cell<T,DESCRIPTOR> fineCell;
    fineLattice.get(finePos, fineCell);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fineCell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho[0], u.data, uSqr)
                       + this->getScalingFactor() * fneq[iPop];
    }

    fineLattice.set(finePos, fineCell);
  }
}


template <typename T, typename DESCRIPTOR>
void computeRestrictedFneq(const SuperLattice2D<T,DESCRIPTOR>& lattice,
                           Vector<int,3> latticeR,
                           T restrictedFneq[DESCRIPTOR::q])
{
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    const auto neighbor = latticeR + Vector<int,3> {0, descriptors::c<DESCRIPTOR>(iPop,0), descriptors::c<DESCRIPTOR>(iPop,1)};
    Cell<T,DESCRIPTOR> cell;
    lattice.get(neighbor, cell);

    T fNeq[DESCRIPTOR::q] {};
    lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq);

    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      restrictedFneq[jPop] += fNeq[jPop];
    }
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    restrictedFneq[iPop] /= DESCRIPTOR::q;
  }
}

template <typename T, typename DESCRIPTOR>
CoarseCoupler2D<T,DESCRIPTOR>::CoarseCoupler2D(
  Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
  Vector<T,2> origin, Vector<T,2> extend):
  Coupler2D<T,DESCRIPTOR>(coarse, fine, origin, extend)
{
  OstreamManager clout(std::cout,"F2C");

  const auto& coarseOrigin = this->getCoarseLatticeR(0);
  const auto& fineOrigin   = this->getFineLatticeR(0);

  clout << "coarse origin: " << coarseOrigin[0] << " " << coarseOrigin[1] << " " << coarseOrigin[2] << std::endl;
  clout << "fine origin:   " << fineOrigin[0]   << " " << fineOrigin[1]   << " " << fineOrigin[2]   << std::endl;
  clout << "coarse size:   " << this->_coarseSize << std::endl;
}

template <typename T, typename DESCRIPTOR>
void CoarseCoupler2D<T,DESCRIPTOR>::couple()
{
  const auto& fineLattice = this->_fine.getSuperLattice();
  auto& coarseLattice = this->_coarse.getSuperLattice();

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (int y=0; y < this->_coarseSize; ++y) {
    const auto& finePos   = this->getFineLatticeR(2*y);
    const auto& coarsePos = this->getCoarseLatticeR(y);

    T fEq[DESCRIPTOR::q] {};
    Cell<T,DESCRIPTOR> fineCell;
    fineLattice.get(finePos, fineCell);
    lbHelpers<T,DESCRIPTOR>::computeFeq(fineCell, fEq);

    T fNeq[DESCRIPTOR::q] {};
    computeRestrictedFneq(fineLattice, finePos, fNeq);

    Cell<T,DESCRIPTOR> coarseCell;
    coarseLattice.get(coarsePos, coarseCell);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      coarseCell[iPop] = fEq[iPop] + this->getInvScalingFactor() * fNeq[iPop];
    }

    coarseLattice.set(coarsePos, coarseCell);
  }
}


}

#endif
