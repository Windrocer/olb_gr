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

#ifndef REFINEMENT_COUPLER_2D_H
#define REFINEMENT_COUPLER_2D_H

#include "grid2D.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
class Coupler2D {
protected:
  Grid2D<T,DESCRIPTOR>& _coarse;
  Grid2D<T,DESCRIPTOR>& _fine;

  const int  _coarseSize;
  const int  _fineSize;
  const bool _vertical;

  const Vector<T,2> _physOrigin;

  const Vector<int,3>& getFineLatticeR(int y) const;
  const Vector<int,3>& getCoarseLatticeR(int y) const;

  T getScalingFactor() const;
  T getInvScalingFactor() const;

private:
  std::vector<Vector<int,3>> _coarseLatticeR;
  std::vector<Vector<int,3>> _fineLatticeR;

public:
  Coupler2D(Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
            Vector<T,2> origin, Vector<T,2> extend);

};

template <typename T, typename DESCRIPTOR>
class FineCoupler2D : public Coupler2D<T,DESCRIPTOR> {
private:
  std::vector<Vector<T,1>>             _c2f_rho;
  std::vector<Vector<T,DESCRIPTOR::d>> _c2f_u;
  std::vector<Vector<T,DESCRIPTOR::q>> _c2f_fneq;

public:
  FineCoupler2D(Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
                Vector<T,2> origin, Vector<T,2> extend);

  void store();
  void interpolate();
  void couple();

};

template <typename T, typename DESCRIPTOR>
class CoarseCoupler2D : public Coupler2D<T,DESCRIPTOR> {
public:
  CoarseCoupler2D(Grid2D<T,DESCRIPTOR>& coarse, Grid2D<T,DESCRIPTOR>& fine,
                  Vector<T,2> origin, Vector<T,2> extend);

  void couple();

};


}

#endif
