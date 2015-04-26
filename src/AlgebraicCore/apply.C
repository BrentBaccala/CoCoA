//   Copyright (c)  2008  John Abbott

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


#include "CoCoA/apply.H"
#include "CoCoA/DenseMatrix.H"


namespace CoCoA
{

  //PROTOTYPE IMPL!!!!!!!!
  //BUG BUG Result ought to be dense/sparse/diag according as M is....
  matrix apply(RingHom phi, ConstMatrixView M)
  {
    CoCoA_ASSERT(domain(phi) == RingOf(M));
    matrix NewM(NewDenseMat(codomain(phi), NumRows(M), NumCols(M)));
    for (long i=0; i < NumRows(M); ++i)
      for (long j=0; j < NumCols(M); ++j)
        SetEntry(NewM, i, j, phi(M(i,j)));
    return NewM;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/apply.C,v 1.3 2014/07/30 14:12:12 abbott Exp $
// $Log: apply.C,v $
// Revision 1.3  2014/07/30 14:12:12  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.2  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.1  2008/11/23 18:33:15  abbott
// First attempt at implementing an apply function for applying maps (RingHoms)
// to pure structures.
//
//
