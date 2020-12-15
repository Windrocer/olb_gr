/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * Template specializations for some computationally intensive LB
 * functions of the header file lbHelpers.h, for the D2Q9 grid.
 */

#ifndef LB_HELPERS_D2Q9_H
#define LB_HELPERS_D2Q9_H

namespace olb {

// Efficient specialization for D2Q9 base lattice
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D2Q9<FIELDS...> > {
  using SpecializedCellBase = CellBase<T,descriptors::D2Q9<FIELDS...>>;

  using SpecializedDescriptor = descriptors::D2Q9<FIELDS...>;

  static T equilibrium(int iPop, T rho, const T u[2], T uSqr)
  {
    typedef descriptors::D2Q9<> L;
    T c_u = descriptors::c<L>(iPop,0)*u[0] + descriptors::c<L>(iPop,1)*u[1];
    return rho * descriptors::t<T,L>(iPop) * (
             1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr )
           - descriptors::t<T,L>(iPop);
  }

  static T equilibriumFirstOrder(int iPop, T rho, const T u[2])
   {
     typedef descriptors::D2Q9<> L;
     T c_u = descriptors::c<L>(iPop,0) * u[0] + descriptors::c<L>(iPop,1) * u[1];

     return rho * descriptors::t<T,L>(iPop) * ( ( T )1 + c_u * descriptors::invCs2<T,L>() ) - descriptors::t<T,L>(iPop);
   }

  static T incEquilibrium(int iPop, const T j[2], const T jSqr, const T pressure)
  {
    typedef descriptors::D2Q9<> L;
    T c_j = descriptors::c<L>(iPop,0)*j[0] + descriptors::c<L>(iPop,1)*j[1];
    return descriptors::t<T,L>(iPop) * (
             3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr )
           - descriptors::t<T,L>(iPop);
  }

  static void computeFneq(SpecializedCellBase const& cell, T fNeq[9], T rho, const T u[2])
  {
    const T uSqr = u[0]*u[0] + u[1]*u[1];
    for (int iPop=0; iPop < 9; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  static T bgkCollision(SpecializedCellBase& cell, T const& rho, const T u[2], T const& omega)
  {
    T uxSqr = u[0]*u[0];
    T uySqr = u[1]*u[1];

    T ux_ = (T)3 * u[0];
    T uy_ = (T)3 * u[1];

    T uxSqr_ = (T)3 * uxSqr;
    T uySqr_ = (T)3 * uySqr;
    T uxSqr__ = (T)3/(T)2 * uxSqr;
    T uySqr__ = (T)3/(T)2 * uySqr;
    T uSqr_ = uxSqr__ + uySqr__;

    T uxPySqr_ = (T)9/(T)2 * (u[0]+u[1])*(u[0]+u[1]);
    T uxMySqr_ = (T)9/(T)2 * (u[0]-u[1])*(u[0]-u[1]);

    T rho_ = (T)4/(T)9 * rho;
    T cf_  = (T)4/(T)9 * (rho-(T)1);

    cell[0] *= (T)1-omega;
    cell[0] += omega*(cf_ + rho_*(- uxSqr__ - uySqr__));

    rho_ = (T)1/(T)9 * rho;
    cf_  = (T)1/(T)9 * (rho-(T)1);

    cell[6] *= (T)1-omega;
    cell[6] += omega*(cf_ + rho_*(ux_ + uxSqr_ - uySqr__));
    cell[8] *= (T)1-omega;
    cell[8] += omega*(cf_ + rho_*(uy_ + uySqr_ - uxSqr__));
    cell[2] *= (T)1-omega;
    cell[2] += omega*(cf_ + rho_*(-ux_ + uxSqr_ - uySqr__));
    cell[4] *= (T)1-omega;
    cell[4] += omega*(cf_ + rho_*(-uy_ + uySqr_ - uxSqr__));

    rho_ = (T)1/(T)36 * rho;
    cf_  = (T)1/(T)36 * (rho-(T)1);

    cell[7] *= (T)1-omega;
    cell[7] += omega*(cf_ + rho_*(ux_ + uy_ + uxPySqr_ - uSqr_));
    cell[1] *= (T)1-omega;
    cell[1] += omega*(cf_ + rho_*(-ux_ + uy_ + uxMySqr_ - uSqr_));
    cell[3] *= (T)1-omega;
    cell[3] += omega*(cf_ + rho_*(-ux_ - uy_ + uxPySqr_ - uSqr_));
    cell[5] *= (T)1-omega;
    cell[5] += omega*(cf_ + rho_*(ux_ - uy_ + uxMySqr_ - uSqr_));

    return uxSqr + uySqr;
  }

  static T spongeCollision(SpecializedCellBase& cell, T const& rho, const T u[2],
		  T const& omega, T const& delta, T const& rho_target, const T u_target[2])
  {
	  T sigma_target = (T)0.3*delta*delta; // here delta is delta/total_distance
	  T feq, f_target;

	  // constants of BGK D2Q9
	  T uxSqr = u[0]*u[0];
	  T uySqr = u[1]*u[1];
	  T ux_ = (T)3 * u[0]; // ux/cs2
	  T uy_ = (T)3 * u[1]; // uy/cs2
	  T uxSqr_ = (T)3 * uxSqr; // uxSqr/cs2
	  T uySqr_ = (T)3 * uySqr; // uySqr/cs2
	  T uxSqr__ = (T)3/(T)2 * uxSqr; // uxSqr_/2
	  T uySqr__ = (T)3/(T)2 * uySqr; // uySqr_/2
	  T uSqr_ = uxSqr__ + uySqr__; //u2/(2cs2)
	  T uxPySqr_ = (T)9/(T)2 * (u[0]+u[1])*(u[0]+u[1]); // (uxSqr+uySqr+2uxuy)/(2cs4)
	  T uxMySqr_ = (T)9/(T)2 * (u[0]-u[1])*(u[0]-u[1]); // (uxSqr+uySqr-2uxuy)/(2cs4)

	  // constants for target correction
	  T uxSqr_target = u_target[0]*u_target[0];
	  T uySqr_target = u_target[1]*u_target[1];
	  T ux_target_ = (T)3 * u_target[0];
	  T uy_target_ = (T)3 * u_target[1];
	  T uxSqr_target_ = (T)3 * uxSqr_target;
	  T uySqr_target_ = (T)3 * uySqr_target;
	  T uxSqr_target__ = (T)3/(T)2 * uxSqr_target;
	  T uySqr_target__ = (T)3/(T)2 * uySqr_target;
	  T uSqr_target_ = uxSqr_target__ + uySqr_target__;
	  T uxPySqr_target_ = (T)9/(T)2 * (u_target[0]+u_target[1])*(u_target[0]+u_target[1]);
	  T uxMySqr_target_ = (T)9/(T)2 * (u_target[0]-u_target[1])*(u_target[0]-u_target[1]);

	  // f0
	  T rho_ = (T)4/(T)9 * rho; //t0*rho
	  T cf_ = (T)4/(T)9 * (rho-(T)1); //t0*rho_bar
	  cell[0] *= (T)1-omega;
	  cell[0] += omega*(cf_ + rho_*(- uxSqr__ - uySqr__));

	  // f0 updated
	  T cf_target_ = (T)4/(T)9 * (rho_target-(T)1);
	  cell[0] += omega*rho*u[0]*u[0]*u[1]*u[1];
	  feq = cf_ + rho_*(- uxSqr__ - uySqr__) + rho*u[0]*u[0]*u[1]*u[1];
	  f_target = cf_target_ + rho_*(- uxSqr_target__ - uySqr_target__) 
		  + rho*u_target[0]*u_target[0]*u_target[1]*u_target[1];
	  cell[0] += -sigma_target*(feq - f_target); 

	  // f2,4,6,8
	  rho_ = (T)1/(T)9 * rho;
	  cf_ = (T)1/(T)9 * (rho-(T)1);
	  cf_target_ = (T)1/(T)9 * (rho_target-(T)1);

	  // f6
	  cell[6] *= (T)1-omega;
	  cell[6] += omega*(cf_ + rho_*(ux_ + uxSqr_ - uySqr__));

	  cell[6] += -omega*(T)0.5*rho*u[0]*u[1]*u[1]*(u[0]+(T)1);
	  feq = cf_ + rho_*(ux_ + uxSqr_ - uySqr__) - (T)0.5*rho*u[0]*u[1]*u[1]*(u[0]+(T)1);
	  f_target = cf_target_ + rho_*(ux_target_ + uxSqr_target_ - uySqr_target__) 
		  - (T)0.5*rho*u_target[0]*u_target[1]*u_target[1]*(u_target[0]+(T)1);
	  cell[6] += -sigma_target*(feq - f_target);

	  // f8
	  cell[8] *= (T)1-omega;
	  cell[8] += omega*(cf_ + rho_*(uy_ + uySqr_ - uxSqr__));

	  cell[8] += -omega*(T)0.5*rho*u[1]*u[0]*u[0]*(u[1]+(T)1);
	  feq = cf_ + rho_*(uy_ + uySqr_ - uxSqr__) - (T)0.5*rho*u[1]*u[0]*u[0]*(u[1]+(T)1);
	  f_target = cf_target_ + rho_*(uy_target_ + uySqr_target_ - uxSqr_target__) 
		  - (T)0.5*rho*u_target[1]*u_target[0]*u_target[0]*(u_target[1]+(T)1);
	  cell[8] += -sigma_target*(feq - f_target);

	  // f2
	  cell[2] *= (T)1-omega;
	  cell[2] += omega*(cf_ + rho_*(-ux_ + uxSqr_ - uySqr__));

	  cell[2] += -omega*(T)0.5*rho*u[0]*u[1]*u[1]*(u[0]-(T)1);
	  feq = cf_ + rho_*(-ux_ + uxSqr_ - uySqr__) - (T)0.5*rho*u[0]*u[1]*u[1]*(u[0]-(T)1);
	  f_target = cf_target_ + rho_*(-ux_target_ + uxSqr_target_ - uySqr_target__) 
		  - (T)0.5*rho*u_target[0]*u_target[1]*u_target[1]*(u_target[0]-(T)1);
	  cell[2] += -sigma_target*(feq - f_target);

	  // f4
	  cell[4] *= (T)1-omega;
	  cell[4] += omega*(cf_ + rho_*(-uy_ + uySqr_ - uxSqr__));

	  cell[4] += -omega*(T)0.5*rho*u[1]*u[0]*u[0]*(u[1]-(T)1);
	  feq = cf_ + rho_*(-uy_ + uySqr_ - uxSqr__) - (T)0.5*rho*u[1]*u[0]*u[0]*(u[1]-(T)1);
	  f_target = cf_target_ + rho_*(-uy_target_ + uySqr_target_ - uxSqr_target__) 
		  - (T)0.5*rho*u_target[1]*u_target[0]*u_target[0]*(u_target[1]-(T)1);
	  cell[4] += -sigma_target*(feq - f_target);


	  // f1,3,5,7
	  rho_ = (T)1/(T)36 * rho;
	  cf_ = (T)1/(T)36 * (rho-(T)1);
	  cf_target_ = (T)1/(T)36 * (rho_target-(T)1);

	  // f7
	  cell[7] *= (T)1-omega;
	  cell[7] += omega*(cf_ + rho_*(ux_ + uy_ + uxPySqr_ - uSqr_));

	  cell[7] += omega*(T)0.25*rho*u[0]*u[1]*(u[0]*u[1]+u[0]+u[1]);
	  feq = cf_ + rho_*(ux_ + uy_ + uxPySqr_ - uSqr_) 
		  + (T)0.25*rho*u[0]*u[1]*(u[0]*u[1]+u[0]+u[1]);
	  f_target = cf_target_ + rho_*(ux_target_ + uy_target_ + uxPySqr_target_ - uSqr_target_) 
		  + (T)0.25*rho*u_target[0]*u_target[1]*(u_target[0]*u_target[1]+u_target[0]+u_target[1]);
	  cell[7] += -sigma_target*(feq - f_target);

	  // f1
	  cell[1] *= (T)1-omega;
	  cell[1] += omega*(cf_ + rho_*(-ux_ + uy_ + uxMySqr_ - uSqr_));

	  cell[1] += omega*(T)0.25*rho*u[0]*u[1]*(u[0]*u[1]+u[0]-u[1]);
	  feq = cf_ + rho_*(-ux_ + uy_ + uxMySqr_ - uSqr_) 
		  + (T)0.25*rho*u[0]*u[1]*(u[0]*u[1]+u[0]-u[1]);
	  f_target = cf_target_ + rho_*(-ux_target_ + uy_target_ + uxMySqr_target_ - uSqr_target_) 
		  + (T)0.25*rho*u_target[0]*u_target[1]*(u_target[0]*u_target[1]+u_target[0]-u_target[1]);
	  cell[1] += -sigma_target*(feq - f_target);

	  // f3
	  cell[3] *= (T)1-omega;
	  cell[3] += omega*(cf_ + rho_*(-ux_ - uy_ + uxPySqr_ - uSqr_));

	  cell[3] += omega*(T)0.25*rho*u[0]*u[1]*(u[0]*u[1]-u[0]-u[1]);
	  feq = cf_ + rho_*(-ux_ - uy_ + uxPySqr_ - uSqr_) 
		  + (T)0.25*rho*u[0]*u[1]*(u[0]*u[1]-u[0]-u[1]);
	  f_target = cf_target_ + rho_*(-ux_target_ - uy_target_ + uxPySqr_target_ - uSqr_target_) 
		  + (T)0.25*rho*u_target[0]*u_target[1]*(u_target[0]*u_target[1]-u_target[0]-u_target[1]);
	  cell[3] += -sigma_target*(feq - f_target);

	  // f5
	  cell[5] *= (T)1-omega;
	  cell[5] += omega*(cf_ + rho_*(ux_ - uy_ + uxMySqr_ - uSqr_));

	  cell[5] += omega*(T)0.25*rho*u[0]*u[1]*(u[0]*u[1]-u[0]+u[1]);
	  feq = cf_ + rho_*(ux_ - uy_ + uxMySqr_ - uSqr_) 
		  + (T)0.25*rho*u[0]*u[1]*(u[0]*u[1]-u[0]+u[1]);
	  f_target = cf_target_ + rho_*(ux_target_ - uy_target_ + uxMySqr_target_ - uSqr_target_) 
		  + (T)0.25*rho*u_target[0]*u_target[1]*(u_target[0]*u_target[1]-u_target[0]+u_target[1]);
	  cell[5] += -sigma_target*(feq - f_target);

	  return uxSqr + uySqr;

  }

  static T incBgkCollision(SpecializedCellBase& cell, T pressure, const T j[2], T omega)
  {
    const T jSqr = util::normSqr<T,descriptors::D2Q9<>::d>(j);
    for (int iPop=0; iPop < descriptors::D2Q9<>::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbHelpers<T,SpecializedDescriptor>::incEquilibrium (
                      iPop, j, jSqr, pressure );
    }
    return jSqr;
  }

  static T constRhoBgkCollision(SpecializedCellBase& cell, T rho, const T u[2], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,descriptors::D2Q9<>::d>(u);
    for (int iPop=0; iPop < descriptors::D2Q9<>::q; ++iPop) {
      T feq = lbHelpers<T,SpecializedDescriptor>::equilibrium(iPop, rho, u, uSqr );
      cell[iPop] = ratioRho*(feq+descriptors::t<T,SpecializedDescriptor>(iPop))
                   -descriptors::t<T,SpecializedDescriptor>(iPop) +
                   ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }


  static void partial_rho(SpecializedCellBase const& cell,
                          T& lineX_P1, T& lineX_0, T& lineX_M1, T& lineY_P1, T& lineY_M1)
  {
    lineX_P1  = cell[5] + cell[6] + cell[7];
    lineX_0   = cell[0] + cell[4] + cell[8];
    lineX_M1  = cell[1] + cell[2] + cell[3];

    lineY_P1  = cell[7] + cell[8] + cell[1];
    lineY_M1  = cell[3] + cell[4] + cell[5];
  }

  static T computeRho(SpecializedCellBase const& cell)
  {
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
            + cell[5] + cell[6] + cell[7] + cell[8] + (T)1;
    return rho;
  }

  static void computeRhoU(SpecializedCellBase const& cell, T& rho, T u[2])
  {
    T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
    partial_rho(cell, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

    rho = lineX_P1 + lineX_0 + lineX_M1 + (T)1;
    T invRho= 1./rho;
    u[0]  = (lineX_P1 - lineX_M1)*invRho;
    u[1]  = (lineY_P1 - lineY_M1)*invRho;
  }

  static void computeRhoJ(SpecializedCellBase const& cell, T& rho, T j[2])
  {
    T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
    partial_rho(cell, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

    rho = lineX_P1 + lineX_0 + lineX_M1 + (T)1;
    j[0]  = (lineX_P1 - lineX_M1);
    j[1]  = (lineY_P1 - lineY_M1);
  }

  static void computeJ(SpecializedCellBase const& cell, T j[2] )
  {
    T lineX_P1, lineX_M1, lineY_P1, lineY_M1;

    lineX_P1  = cell[5] + cell[6] + cell[7];
    lineX_M1  = cell[1] + cell[2] + cell[3];
    lineY_P1  = cell[7] + cell[8] + cell[1];
    lineY_M1  = cell[3] + cell[4] + cell[5];

    j[0]  = (lineX_P1 - lineX_M1);
    j[1]  = (lineY_P1 - lineY_M1);
  }

  static void computeStress(SpecializedCellBase const& cell, T rho, const T u[2], T pi[3])
  {
    typedef descriptors::D2Q9<> L;
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices2D" is not sufficient
    using util::tensorIndices2D::xx;
    using util::tensorIndices2D::yy;
    using util::tensorIndices2D::xy;

    T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
    partial_rho(cell, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

    pi[xx] = lineX_P1+lineX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = lineY_P1+lineY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[xy] = -cell[1] + cell[3] - cell[5] + cell[7]   - rho*u[0]*u[1];
  }

  static void computeAllMomenta(SpecializedCellBase const& cell, T& rho, T u[2], T pi[3] )
  {
    typedef descriptors::D2Q9<> L;
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices2D" is not sufficient
    using util::tensorIndices2D::xx;
    using util::tensorIndices2D::yy;
    using util::tensorIndices2D::xy;

    T lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1;
    partial_rho(cell, lineX_P1, lineX_0, lineX_M1, lineY_P1, lineY_M1);

    rho = lineX_P1 + lineX_0 + lineX_M1 + (T)1;

    T rhoU0 = (lineX_P1 - lineX_M1);
    T rhoU1 = (lineY_P1 - lineY_M1);
    u[0]  = rhoU0/rho;
    u[1]  = rhoU1/rho;

    pi[xx] = lineX_P1 + lineX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rhoU0*u[0];
    pi[yy] = lineY_P1 + lineY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rhoU1*u[1];
    pi[xy] = -cell[1] + cell[3] - cell[5] + cell[7]        - rhoU0*u[1];
  }

  static void modifyVelocity(SpecializedCellBase const& cell, const T newU[2])
  {
    T rho, oldU[2];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,2>(oldU);
    const T newUSqr = util::normSqr<T,2>(newU);
    for (int iPop=0; iPop<9; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  //struct lbHelpers<D2Q9>

// Efficient specialization for D2Q9 lattice with force
template<typename T>
struct lbExternalHelpers<T, descriptors::D2Q9<descriptors::FORCE>> {

  static void addExternalForce(
    Cell<T,descriptors::D2Q9<descriptors::FORCE>>& cell,
    const T u[descriptors::D2Q9<descriptors::FORCE>::d], T omega, T amplitude)
  {
    const T* force = cell.template getFieldPointer<descriptors::FORCE>();
    const T  mu = amplitude*((T)1-omega/(T)2);

    cell[0] += mu *(T)4/(T)3  *( force[0] * (-  u[0]             ) +
                                 force[1] * (        -   u[1]    )   );
    cell[1] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] - 3*u[1] - 1) +
                                 force[1] * (-3*u[0] + 2*u[1] + 1)   );
    cell[2] += mu *(T)1/(T)3  *( force[0] * ( 2*u[0]          - 1) +
                                 force[1] * (        -   u[1]    )   );
    cell[3] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] + 3*u[1] - 1) +
                                 force[1] * ( 3*u[0] + 2*u[1] - 1)   );
    cell[4] += mu *(T)1/(T)3  *( force[0] * (-  u[0]             ) +
                                 force[1] * (        + 2*u[1] - 1)   );
    cell[5] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] - 3*u[1] + 1) +
                                 force[1] * (-3*u[0] + 2*u[1] - 1)   );
    cell[6] += mu *(T)1/(T)3  *( force[0] * ( 2*u[0]          + 1) +
                                 force[1] * (        -   u[1]    )   );
    cell[7] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] + 3*u[1] + 1) +
                                 force[1] * ( 3*u[0] + 2*u[1] + 1)   );
    cell[8] += mu *(T)1/(T)3  *( force[0] * (-  u[0]             ) +
                                 force[1] * (        + 2*u[1] + 1)   );
  }
};

// Efficient specialization for D2Q9 lattice and for forced D2Q9 lattice
//   (operations applying to the whole lattice)

template<typename T>
struct lbLatticeHelpers<T, descriptors::D2Q9<>> {

  static void swapAndStreamCell (
    Cell<T,descriptors::D2Q9<>> **grid,
    int iX, int iY, int nX, int nY, int iPop, T& fTmp )
  {
    fTmp                 = grid[iX][iY][iPop];
    grid[iX][iY][iPop]   = grid[iX][iY][iPop+4];
    grid[iX][iY][iPop+4] = grid[nX][nY][iPop];
    grid[nX][nY][iPop]   = fTmp;
  }

  static void swapAndStream2D (
    Cell<T,descriptors::D2Q9<>> **grid, int iX, int iY )
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iX-1, iY+1, 1, fTmp);
    swapAndStreamCell(grid, iX, iY, iX-1, iY,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iX-1, iY-1, 3, fTmp);
    swapAndStreamCell(grid, iX, iY, iX,   iY-1, 4, fTmp);
  }
};

template<typename T>
struct lbLatticeHelpers<T, descriptors::D2Q9<descriptors::FORCE>> {

  static void swapAndStreamCell (
    Cell<T,descriptors::D2Q9<descriptors::FORCE>> **grid,
    int iX, int iY, int nX, int nY, int iPop, T& fTmp )
  {
    fTmp                 = grid[iX][iY][iPop];
    grid[iX][iY][iPop]   = grid[iX][iY][iPop+4];
    grid[iX][iY][iPop+4] = grid[nX][nY][iPop];
    grid[nX][nY][iPop]   = fTmp;
  }

  static void swapAndStream2D (
    Cell<T,descriptors::D2Q9<descriptors::FORCE>> **grid, int iX, int iY )
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iX-1, iY+1, 1, fTmp);
    swapAndStreamCell(grid, iX, iY, iX-1, iY,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iX-1, iY-1, 3, fTmp);
    swapAndStreamCell(grid, iX, iY, iX,   iY-1, 4, fTmp);
  }

};

}  // namespace olb

#endif
