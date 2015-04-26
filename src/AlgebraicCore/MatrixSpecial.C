//   Copyright (c)  2010  Anna Bigatti

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


#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixSpecial.H"
#include "CoCoA/PolyRing.H" // for jacobian
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/utils.H" // for len 
#include "CoCoA/VectorOperations.H" // for HasUniqueOwner

#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{


  //---------  jacobian matrix

  matrix jacobian_aux(const PolyRing& P, const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    const long LenPolys = len(polys);
    const long LenInds = len(inds);
    matrix ans = NewDenseMat(P, LenPolys, LenInds);
    for (long i=0; i < LenPolys; ++i)
      for (long j=0; j < LenInds; ++j)
        SetEntry(ans, i,j, deriv(polys[i], inds[j]));
    return ans;
  }


  matrix jacobian(const std::vector<RingElem>& polys)
  {
    if (polys.empty()) CoCoA_ERROR(ERR::BadArg, "jacobian(polys)");
    if (!HasUniqueOwner(polys)) CoCoA_ERROR(ERR::MixedRings, "jacobian(polys)");
    const PolyRing P = owner(polys[0]);
    return jacobian_aux(P, polys, indets(P));
  }
  

  matrix jacobian(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    if (len(inds)==0 && len(polys)==0) CoCoA_ERROR(ERR::BadArg, "jacobian");
    if (!HasUniqueOwner(inds)) CoCoA_ERROR(ERR::MixedRings, "jacobian");
    if (!HasUniqueOwner(polys)) CoCoA_ERROR(ERR::MixedRings, "jacobian");
    if (len(inds)==0)
      return NewDenseMat(owner(polys[0]), len(polys), 0);
    const PolyRing P(owner(inds[0]));
    if (polys.empty()) return NewDenseMat(P, 0, len(inds));
    if (owner(polys[0]) != P) CoCoA_ERROR(ERR::MixedRings, "jacobian");
    return jacobian_aux(P, polys, inds);
  }
  

  //---------  tensor matrix

  matrix TensorMat(ConstMatrixView A, ConstMatrixView B)
  {
    if (RingOf(A) != RingOf(B))
      CoCoA_ERROR(ERR::MixedRings, "TensorMat(A,B)");
    if (NumCols(A) > numeric_limits<long>::max() /NumCols(B)) //avoid overflow  
      CoCoA_ERROR(ERR::BadColIndex, "TensorMat(A,B)");
    if (NumRows(A) > numeric_limits<long>::max() /NumRows(B)) //avoid overflow  
      CoCoA_ERROR(ERR::BadRowIndex, "TensorMat(A,B)");
    matrix ans = NewDenseMat(RingOf(A), NumRows(A)*NumRows(B), NumCols(A)*NumCols(B));
    for (long iA=0; iA < NumRows(A); ++iA)
      for (long jA=0; jA < NumCols(A); ++jA)
        for (long iB=0; iB < NumRows(B); ++iB)
          for (long jB=0; jB < NumCols(B); ++jB)
            SetEntry(ans, iA*NumRows(B)+iB, jA*NumCols(B)+jB, A(iA, jA)*B(iB, jB));
    return ans;
  }
  
  //---------  Sylvester matrix

//   matrix sylvester(ConstRefRingElem f, ConstRefRingElem g, ConstRefRingElem x)
//   {
//     const PolyRing P = owner(x);
//     if (owner(f) != P || owner(g) != P)
//       CoCoA_ERROR(ERR::MixedRings, "sylvester(f, g, x)");
//     long index;
//     if (!IsIndet(index, x)) CoCoA_ERROR(ERR::NotIndet, "sylvester(f, g, x)");
//     const long degf = MaxExponent(f, index);
//     const long degg = MaxExponent(g, index);
    
//     matrix ans = NewDenseMat(P, degf+degg, degf+degg);
//     vector<RingElem> cf(degf+1, zero(P));
//     vector<RingElem> cg(degg+1, zero(P));

//     return ans;
//   }
  

// VandermondeMatrix
// HessianMatrix
// HilbertMatrix
// HilbertInverseMatrix
// ToeplitzMatrix
// WronskianMatrix


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixSpecial.C,v 1.10 2014/07/31 14:45:17 abbott Exp $
// $Log: MatrixSpecial.C,v $
// Revision 1.10  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.9  2014/07/30 14:06:39  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.8  2014/07/14 15:06:52  abbott
// Summary: Added include of UtilsTemplate.H
// Author: JAA
//
// Revision 1.7  2014/07/07 12:28:26  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.6  2014/04/30 16:08:58  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.5  2011/03/23 17:30:12  bigatti
// -- started Sylvester matrix
//
// Revision 1.4  2011/03/21 07:58:29  bigatti
// -- added TensorMat
// -- changed size into len
//
// Revision 1.3  2011/02/22 17:04:41  bigatti
// -- fixed jacobian
//
// Revision 1.2  2011/02/10 15:46:21  bigatti
// -- removed debug printing
//
// Revision 1.1  2011/02/10 15:30:14  bigatti
// -- first import: jacobian
//
