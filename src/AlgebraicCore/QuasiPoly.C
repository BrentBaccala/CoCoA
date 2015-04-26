//   Copyright (c)  2014  John Abbott

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/QuasiPoly.H"
#include "CoCoA/error.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/VectorOperations.H"
#include "CoCoA/utils.H"


namespace CoCoA
{

  QuasiPoly::QuasiPoly(const std::vector<RingElem>& v):
      myConstituents(v)
  {
    // Simply check that the input v is valid...
    const char* const FnName = "QuasiPoly ctor";
    if (v.empty()) CoCoA_ERROR(ERR::Empty, FnName);
    const ring& R = owner(v[0]);
    const long n = len(v);
    for (long i=1; i < n; ++i)
      if (owner(v[i]) != R) CoCoA_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(R)) CoCoA_ERROR(ERR::NotPolyRing, FnName);
  }


  const std::vector<RingElem>& constituents(const QuasiPoly& p)
  {
    return p.myConstituents;
  }


  RingElem QuasiPoly::operator()(const MachineInt& n) const
  {
    return operator()(BigInt(n));
  }

  RingElem QuasiPoly::operator()(const BigInt& N) const
  {
    const long i = ConvertTo<long>(LeastNNegRemainder(N, len(myConstituents)));
    const PolyRing P = owner(myConstituents[0]);
    const RingHom EvalAtN = EvalHom(P, N);
    return EvalAtN(myConstituents[i]);
  }

  std::ostream& operator<<(std::ostream& out, const QuasiPoly& p)
  {
    out << "QuasiPoly(" << constituents(p) << ")";
    return out;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QuasiPoly.C,v 1.3 2014/07/31 14:45:18 abbott Exp $
// $Log: QuasiPoly.C,v $
// Revision 1.3  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.2  2014/07/14 13:16:46  abbott
// Summary: Minor cleaning
// Author: JAA
//
// Revision 1.1  2014/07/14 10:04:29  abbott
// Summary: Added new impl of quasipolynomials
// Author: JAA
//
