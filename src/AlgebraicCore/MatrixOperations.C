//   Copyright (c)  2005,2008  John Abbott

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


#include "CoCoA/MatrixOperations.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/IntOperations.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/apply.H"
#include "CoCoA/bool3.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/utils.H" // for len

#include <algorithm>
using std::min;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{

  // Naive dense matrix multiplication
  // Currently just creates a DenseMat to contain the answer.
  // BUG: must make this behave "intelligently" when multiplying two sparse matrices
  // (what does "sparse" mean?  e.g. consider  transpose(SparseMat)*SparseMat
  //  or even DiagMat(...)*AnyMat, or even ZeroMat*AnyMat,... lots of cases!!!??? BUG
  matrix operator*(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    if (NumCols(Mleft) != NumRows(Mright))
      CoCoA_ERROR(ERR::BadMatrixSize, "Mat1*Mat2");
    const ring R = RingOf(Mleft);
    if (RingOf(Mright) != R)
      CoCoA_ERROR(ERR::MixedRings, "Mat1*Mat2");
    const long N = NumCols(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mright);
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
      {
        RingElem tmp(R);
        for(long k=0; k < N; ++k)
          tmp += Mleft(i,k)*Mright(k,j);
        SetEntry(ans, i, j, tmp);
      }
    return ans;
  }


  matrix operator+(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const ring R = RingOf(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1+Mat2");
    if (NumCols(Mright) != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1+Mat2");
    if (RingOf(Mright) != R) CoCoA_ERROR(ERR::MixedRings, "Mat1+Mat2");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)+Mright(i,j));
    return ans;
  }


  matrix operator-(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const ring R = RingOf(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1-Mat2");
    if (NumCols(Mright) != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "Mat1-Mat2");
    if (RingOf(Mright) != R) CoCoA_ERROR(ERR::MixedRings, "Mat1-Mat2");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)-Mright(i,j));
    return ans;
  }


  // Innermost loop is ugly but rather faster than a "clean" implementation.
  void mul(matrix& lhs, ConstMatrixView M1, ConstMatrixView M2)
  {
    const ring& R = RingOf(lhs);
    if (NumCols(M1) != NumRows(M2))
      CoCoA_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
    if ( NumRows(lhs) != NumRows(M1) || NumCols(lhs) != NumCols(M2) )
      CoCoA_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
    if (RingOf(M1) != R || RingOf(M2) != R)
      CoCoA_ERROR(ERR::MixedRings, "mul(Mat1, Mat2)");
    // Use of the temporary ans avoids aliasing problems and makes the code exception safe.
    matrix ans(lhs->myZeroClone(RingOf(lhs), NumRows(M1), NumCols(M2)));
    RingElem sum(R), prod(R);
    for (long i=0; i < NumRows(M1); ++i)
      for (long j=0; j < NumCols(M2); ++j)
      {
        sum = 0;
        for (long k=0; k < NumCols(M1); ++k)
        {
          // Next 2 lines just do:  sum += M1(i,k) * M2(k,j);
          R->myMul(raw(prod), raw(M1(i,k)), raw(M2(k,j)));
          R->myAdd(raw(sum), raw(sum), raw(prod));
        }
        SetEntry(ans, i, j, sum);
      }

    // The answer is in ans; now swap it with the entries of lhs.
    swap(lhs, ans);
  }

  matrix operator*(ConstRefRingElem x, ConstMatrixView M)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "x*M");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, x*M(i,j));
    return ans;
  }
  
  matrix operator*(const BigRat& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const BigInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const MachineInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }


  // Separate impl in case ring is not commutative
  matrix operator*(ConstMatrixView M, ConstRefRingElem x)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "x*M");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)*x);
    return ans;
  }

  matrix operator*(ConstMatrixView M, const BigRat& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const BigInt& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const MachineInt& x)
  { return M * RingElem(RingOf(M),x); }


  matrix operator-(const ConstMatrixView& M)
  {
    return RingElem(RingOf(M),-1)*M;
  }


  matrix operator/(ConstMatrixView M, ConstRefRingElem x)
  {
    if (IsZeroDivisor(x)) CoCoA_ERROR(ERR::DivByZero, "Mat/x");
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (owner(x) != R) CoCoA_ERROR(ERR::MixedRings, "Mat/x");
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)/x);
    return ans;
  }

  matrix operator/(ConstMatrixView M, const BigRat& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const BigInt& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const MachineInt& x)
  { return M / RingElem(RingOf(M),x); }


  matrix power(ConstMatrixView M, long n)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (Nrows != Ncols) CoCoA_ERROR(ERR::BadMatrixSize, "power(M,n)");
    if (n == numeric_limits<long>::min()) CoCoA_ERROR(ERR::ExpTooBig, "power(M,n)");
    if (n < 0) return power(inverse(M), -n); // cannot overflow because we have excluded n == MinLong
    if (n == 0) return NewDenseMat(IdentityMat(R,Nrows));

    // An iterative implementation of binary powering.
    long bit = 1; while (bit <= n/2) bit <<= 1;
    matrix ans = NewDenseMat(M);
    while (bit > 1)
    {
      mul(ans,ans,ans);//ans *= ans;
      bit >>= 1;
      if (n&bit) mul(ans,ans,M);//ans *= M;
    }
    return ans;
  }


  matrix power(ConstMatrixView M, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_ERROR(ERR::ExpTooBig, "power(M,N)");
    return power(M, n);
  }


  // The square of the Frobenius norm of a matrix.
  RingElem FrobeniusNorm2(ConstMatrixView A)
  {
    RingElem FrNorm2 = zero(RingOf(A));
    for (long i=0; i < NumRows(A); ++i)
      for (long j=0; j < NumCols(A); ++j)
        FrNorm2 += A(i,j)*A(i,j);
    return FrNorm2;
  }


  // Compute the induced infty-norm of a matrix
  RingElem OperatorNormInfinity(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsOrderedDomain(R))
      CoCoA_ERROR(ERR::NotOrdDom, "OperatorNormInfinity(Mat)");
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    RingElem MaxRowSize = zero(R);

    for (long i=0; i < Nrows; ++i)
    {
      RingElem RowSize(zero(R));
      for (long j=0; j < Ncols; ++j) {	RowSize += abs(M(i,j)); }
      if (RowSize > MaxRowSize)
        MaxRowSize = RowSize;
    }
    return MaxRowSize;
  }

  RingElem OperatorNorm1(ConstMatrixView M)
  {
    return OperatorNormInfinity(transpose(M));
  }


  RingElem det(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det(Mat)");
    RingElem d(RingOf(M));
    M->myDet(d);
    return d;
  }


  long rank(const ConstMatrixView& M)
  {
    if (!IsIntegralDomain(RingOf(M)))
      CoCoA_ERROR(ERR::NotIntegralDomain, "rank(Mat)");
    return M->myRank();
  }


  matrix inverse(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "inverse(Mat)");
    return InverseByGauss(M);
  }


  matrix adjoint(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "adjoint(Mat)");
    if (!IsIntegralDomain(RingOf(M)))
      return AdjointByDetOfMinors(M);
    else
      return AdjointByInverse(M);
  }


  // Should restriction to full rank be in the name?
  matrix PseudoInverse(ConstMatrixView M)
  {
    // BUG??? Would it make sense to generalize to non fields???
    const ring R = RingOf(M);
    if (!IsField(R))
      CoCoA_ERROR(ERR::NotField, "PseudoInverse(Mat)");

    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    const long rk = rank(M);
    if (rk < Nrows && rk < Ncols) 
      CoCoA_ERROR(ERR::NYI, "PseudoInverse of non full rank matrix");

    // Easy case: a square full rank matrix
    if (Nrows == Ncols)
      return inverse(M);

    if (Nrows < Ncols)
      return transpose(M)*inverse(M*transpose(M));
    else
      return inverse(transpose(M)*M)*transpose(M);
  }


  matrix LinSolve(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R || NumRows(M) != NumRows(rhs))
      CoCoA_ERROR(ERR::BadArg, "LinSolve");
    if (IsField(R)) return LinSolveByGauss(M, rhs);
    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_ERROR(ERR::NYI, "LinSolve over non-field, non-gcddomain, non-polynomial-ring");
    return LinSolve(M,rhs); // never reached -- just to keep compiler quiet
  }


  matrix LinKer(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (IsField(R)) return LinKerByGauss(M);
    //    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    //    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_ERROR(ERR::NYI, "LinKer over non-field");
    return LinKer(M); // never reached -- just to keep compiler quiet
  }


  // WARNING!! Pivot selection strategy is simple rather than clever!
  long RankAndGauss(matrix& M, const int ToDoCols)
  {
    const ring R = RingOf(M);
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "gauss");  
    if (ToDoCols > NumCols(M)) CoCoA_ERROR(ERR::BadColIndex, "gauss");  
    const long Mrows = NumRows(M);

    long rank = 0;
    for (long j=0; j < ToDoCols; ++j)
    {
      // Look for a pivot in col j.
      long PivotRow=rank;
      while (PivotRow < Mrows && M(PivotRow, j) == 0)
        ++PivotRow;
      if (PivotRow == Mrows) continue; // col was zero, try next col.
      if (PivotRow != rank)  SwapRows(M, rank, PivotRow);
      M->myRowMul(rank, 1/M(rank, j));  // make pivot entry = 1
      for (long i=0; i < Mrows; ++i)
      {
        if (i == rank) continue;
        if (M(i, j) == 0) continue;
        M->myAddRowMul(i, rank, -M(i,j));
      }
      ++rank;
    }
    return rank;
  }


  matrix LinSolveByGauss(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R || NumRows(M) != NumRows(rhs))
      CoCoA_ERROR(ERR::BadArg, "LinSolveByGauss");
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "LinSolveByGauss");
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);
    const long RHScols = NumCols(rhs);
    matrix tmp = NewDenseMat(ConcatHor(M, rhs));

    // Do row redn by Gauss
    long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, RHScols); // initially full of zeroes
    long col=0;
    for (long i=0; i < rank; ++i)
    {
      while (tmp(i,col) == 0) { ++col; }
      for (long j=0; j < RHScols; ++j)
        SetEntry(ans, col, j, tmp(i, j+Mcols));
    }
    for (long i=rank; i < Mrows; ++i)
    {
      for (long j=0; j < RHScols; ++j)
        {
          if (tmp(i, j+Mcols) != 0)
            return NewDenseMat(R,0,0); // to indicate that no soln exists
        }
    }
    return ans;
  }


  matrix LinKerByGauss(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsField(R)) CoCoA_ERROR(ERR::NotField, "LinKerByGauss");
    matrix tmp = NewDenseMat(M);
  
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);

    // Do row redn by Gauss
    long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, Mcols-rank); // initially full of zeroes
    long row=0;
    long anscol=0;
    vector<long> PivotCols; // the i-th pivot is in col PivotCols[i]
    ConstMatrixView Z = ZeroMat(R, std::max(Mcols-Mrows, (long)0), Mcols);
    ConstMatrixView SqTmp = ConcatVer(tmp, Z); // make it square
    for (long j=0; j < Mcols; ++j) // we consider only Mcols x Mcols
      if (SqTmp(row,j) != 0) // j-th col with pivot
      {
        PivotCols.push_back(j);
        ++row;
      }
      else // j-th col without pivot
      {
        for (long i=0; i < len(PivotCols); ++i)  // "copy" j-th column
          SetEntry(ans, PivotCols[i], anscol, SqTmp(i, j));
        SetEntry(ans, j, anscol, -1);
        ++anscol;
      }
    return ans;
  }


  matrix LinSolveByHNF(ConstMatrixView M, ConstMatrixView rhs)
  {
    // HNF works only for PIDs: i.e. ZZ or k[x]
    // NB Could work in k[x,y,z] if M is univariate!
    CoCoA_ERROR(ERR::NYI, "LinSolveByHNF");
    return LinSolveByHNF(M,rhs); // never reached -- just to keep compiler quiet
  }

  matrix LinSolveByModuleRepr(ConstMatrixView M, ConstMatrixView rhs)
  {
    // Works in k[x,y,z] where k is a field.  Probably slow.
    CoCoA_ERROR(ERR::NYI, "LinSolveByModuleRepr");
    return LinSolveByModuleRepr(M,rhs); // never reached -- just to keep compiler quiet
  }

  /***************************************************************************/


  void det2x2(RingElem& d, ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det2x2(d,Mat)");
    if (NumRows(M) != 2)
      CoCoA_ERROR(ERR::BadRowIndex, "det2x2(d,Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    d = M(0,0)*M(1,1) - M(0,1)*M(1,0);
  }
  

  void det3x3(RingElem& d, ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "det3x3(d,Mat)");
    if (NumRows(M) != 3)
      CoCoA_ERROR(ERR::BadRowIndex, "det3x3(d,Mat)");
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    d = M(2,2) *(M(0,0)*M(1,1) - M(0,1)*M(1,0)) 
      - M(2,1) *(M(0,0)*M(1,2) - M(0,2)*M(1,0))
      + M(2,0) *(M(0,1)*M(1,2) - M(0,2)*M(1,1));
  }
  

  // Known defect: this algm is valid only over IntegralDomains
  void DetByGauss(RingElem& d, ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "DetByGauss(d,Mat)");
    const long N = NumRows(M);
    const ring R(RingOf(M));
    // ??? this function works only within integral domains!
    if (!IsIntegralDomain(R))
      CoCoA_ERROR(ERR::NotIntegralDomain, "DetByGauss(d,Mat)");
// Is it worth adding code for this special case?  General code works fine.
//     // Handle 0x0 matrix specially: its det is 1.
//     if (N == 0) { d = one(R); return; }
    const ring K((IsField(R) ? R : NewFractionField(R)));
    matrix Gss(apply(CanonicalHom(R,K), M));
    RingElem c(K);
    RingElem determinant(one(K));
    for (long col=0; col < N; ++col)
    {
      // Pick a good pivot row
      long PivotRow = -1;
      for (long r=col ; r < N; ++r)
      {
        if (!IsZero(Gss(r,col))) { PivotRow = r; break; }
      }

      if (PivotRow == -1)
      {
        d = zero(R);
        return;
      }
      if (PivotRow != col)
      {
        Gss->mySwapRows(PivotRow,col);
        determinant *= -1;
      }
      c = Gss(col,col);
      determinant *= c;
      for (long i=col+1; i < N; ++i)
        Gss->myAddRowMul(i, col, -Gss(i,col)/c);
    }
    if (R==K) swap(d, determinant);
    else d = num(determinant)/den(determinant);
  }


  long RankByGauss(std::vector<long>& IndepRows, ConstMatrixView M)
  {
    // ??? this function works only within integral domains!
    const ring R(RingOf(M));
    if (!IsIntegralDomain(R)) CoCoA_ERROR(ERR::NotIntegralDomain, "RankByGauss(v,Mat)  over non-integral domain");
    const ring K((IsField(R) ? R : NewFractionField(R)));
    const RingHom R2K(R==K ? IdentityHom(K) : EmbeddingHom(K));
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix Gss(NewDenseMat(K, Nrows, Ncols));

    // below copy M via RingHom into Gss
    {
      // This will eventually become a separate function
      for (long row=0; row < Nrows; ++row)
        for (long col=0; col < Ncols; ++col)
          SetEntry(Gss, row, col, R2K(M(row,col)));
    }
    IndepRows.clear();
    RingElem pivot(K);
    long row=0;
    for (long col=0; col < Ncols; ++col)
    {
      if (IsZero(Gss(row,col)))
      {
        long i=row+1;
        for ( ; i < Nrows; ++i)
          if (!IsZero(Gss(i,col)))
          {
            Gss->mySwapRows(i,row);
            break;
          }
        if (i==Nrows) continue;
      }
      IndepRows.push_back(row);
      pivot = Gss(row,col);
      for (long i=row+1; i < Nrows; ++i)
        Gss->myAddRowMul(i, row, -Gss(i,col)/pivot);
      ++row;
      if (row == Nrows) return row;
    }
    return row;
  }


  matrix InverseByGauss(ConstMatrixView M)
  {
    // this code works only if the base ring is a field
    CoCoA_ASSERT(IsField(RingOf(M)));
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "InverseByGauss(Mat)");

    const long N = NumRows(M);
    matrix Gss(NewDenseMat(M));
    matrix inv = NewDenseMat(IdentityMat(RingOf(M), N));
    RingElem c(RingOf(inv));
    for (long j=0; j < N; ++j)
    {
      if (IsZero(Gss(j,j)))
      {
        long i=j+1;
        for ( ; i < N; ++i)
          if (!IsZero(Gss(i,j))) break;
        if (i == N) CoCoA_ERROR(ERR::NotInvMatrix, "InverseByGauss(Mat)");
        Gss->mySwapRows(i,j);
        inv->mySwapRows(i,j);
      }
      c = 1/Gss(j,j);
      Gss->myRowMul(j, c);
      inv->myRowMul(j, c);
      for (long i=0; i < N; ++i)
        if (i != j)
        {
          c = -Gss(i,j);
          Gss->myAddRowMul(i, j, c); // AddRowMul(Gss, i, j, c);
          inv->myAddRowMul(i, j, c);
        }
    }
    return inv;
  }


//////////////////////////////////////////////////////////////////
// Bareiss

  namespace
  {
    RingElem DetByBareiss_ZZ(const ConstMatrixView& M)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)));
      const int n = NumRows(M);
      CoCoA_ASSERT(NumCols(M) == n);
      const ring R = RingOf(M);
      vector< vector< BigInt > > M2(n, vector<BigInt>(n));
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          IsInteger(M2[i][j], M(i,j)); // ignore return value (must be true)
      BigInt d(1);
      int sign = 1;
      for (int k=0; k < n-1; ++k)
      {
        // Find a non-zero pivot in k-th column
        int row = -1;
        for (int i=k; i < n; ++i)
          if (!IsZero(M2[i][k])) { row = i; break; }
        if (row == -1) return zero(R);
        if (row != k) { swap(M2[k], M2[row]); sign = -sign; }
        BigInt tmp; // temporary workspace used in inner loop below
        for (int i=k+1; i < n; ++i)
          for (int j=k+1; j < n; ++j)
          {
//  This block effectively does the following, but is usefully faster (about 4x)
//            M2[i][j] = (M2[i][j]*M2[k][k]-M2[i][k]*M2[k][j])/d;
            mpz_mul(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(M2[k][k]));
            mpz_mul(mpzref(tmp), mpzref(M2[i][k]), mpzref(M2[k][j]));
            mpz_sub(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(tmp));
            mpz_divexact(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(d));
          }
        d = M2[k][k];
      }
      return RingElem(R,sign*M2[n-1][n-1]);
    }
  } // end of anonymous namespace


  RingElem DetByBareiss(const ConstMatrixView& M)
  {
    const ring R = RingOf(M);
    if (IsZZ(R)) return DetByBareiss_ZZ(M); // a bit faster
    CoCoA_ASSERT(IsIntegralDomain(R));
    const int n = NumRows(M);
    CoCoA_ASSERT(NumCols(M) == n);
    matrix M2 = NewDenseMat(M);
    RingElem d = one(R);
    int sign = 1;
    for (int k=0; k < n-1; ++k)
    {
      // Find a non-zero pivot in k-th column
      int row = -1;
      for (int i=k; i < n; ++i)
        if (!IsZero(M2(i,k))) { row = i; break; }
      if (row == -1) return zero(R);
      if (row != k) { M2->mySwapRows(k, row); sign = -sign; }
      // Now use pivot row to reduce all lower rows
      for (int i=k+1; i < n; ++i)
        for (int j=k+1; j < n; ++j)
        {
          SetEntry(M2,i,j,(M2(i,j)*M2(k,k)-M2(i,k)*M2(k,j))/d );
        }
      d = M2(k,k);
    }
    if (sign == 1)
      return M2(n-1,n-1);
    return -M2(n-1,n-1);
  }


//////////////////////////////////////////////////////////////////

  matrix AdjointByInverse(ConstMatrixView M)
  {
    // This code works only if the matrix is invertible!
    // It is morally equivalent to swap(lhs, inverse(M)*det(M));
    CoCoA_ASSERT(IsIntegralDomain(RingOf(M))); // needed for this method, not mathematically necessary BUG BUG!!!
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    if (IsField(RingOf(M)))
    {
      RingElem d = det(M);
      if (!IsZero(d)) return d*inverse(M);
      else return AdjointByDetOfMinors(M);
    }
    FractionField K(NewFractionField(RingOf(M)));
    RingHom R2K(EmbeddingHom(K));
    matrix ans_K = adjoint(apply(R2K, M));
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix ans = NewDenseMat(RingOf(M), Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i,j, num(ans_K(i,j)));
    return ans;
  }


  // Simple (but probably not very fast).
  matrix AdjointByDetOfMinors(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    matrix adj = NewDenseMat(RingOf(M), n,n);
    vector<long> rows(n-1);
    for (long i=0; i < n-1; ++i) { rows[i] = i+1; }

    vector<long> cols(n-1);
    for (long i=0; i < n; ++i)
    {
      if (i > 0) rows[i-1] = i-1;
      for (long j=0; j < n-1; ++j) { cols[j] = j+1; }
      for (long j=0; j < n; ++j)
      {
        if (j > 0) cols[j-1] = j-1;
        SetEntry(adj, j,i, power(-1, i+j)*det(submat(M,rows,cols)));
      }
    }
    return adj;
  }


  bool IsZero(const ConstMatrixView& M)
  { return M == ZeroMat(RingOf(M), NumRows(M), NumCols(M)); }


  bool IsZeroRow(const ConstMatrixView& M, long i)
  {
    M->myCheckRowIndex(i, "IsZeroRow(M)");
    return M->myIsZeroRow(i);
  }


  bool IsZeroCol(const ConstMatrixView& M, long j)
  {
    M->myCheckColIndex(j, "IsZeroCol(M)");
    return M->myIsZeroCol(j);
  }


  bool IsSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsSymmetric");
    return M->IamSymmetric();
  }


  bool IsAntiSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsAntiSymmetric");
    return M->IamAntiSymmetric();
  }


  bool IsDiagonal(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))
      CoCoA_ERROR(ERR::NotSquareMatrix, "IsDiagonal");
    return M->IamDiagonal();
  }

  
  bool IsMat0x0(const ConstMatrixView& M)
  {
    return NumRows(M) == 0 && NumCols(M) == 0;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOperations.C,v 1.6 2014/08/26 12:55:58 abbott Exp $
// $Log: MatrixOperations.C,v $
// Revision 1.6  2014/08/26 12:55:58  abbott
// Summary: Cleaned up DetByGauss; added DetByBareiss
// Author: JAA
//
// Revision 1.5  2014/07/30 14:06:24  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.4  2014/07/08 08:35:16  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.3  2014/07/07 12:23:22  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.2  2014/04/17 13:38:47  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.1  2014/04/11 15:42:37  abbott
// Summary: Renamed from MatrixArith
// Author: JAA
//
// Revision 1.46  2014/04/08 15:40:50  abbott
// Summary: Replaced test for IsZero by IsZeroDivisor
// Author: JAA
//
// Revision 1.45  2014/01/16 16:10:56  abbott
// Removed a blank line.
//
// Revision 1.44  2012/11/23 17:28:26  abbott
// Modified LinSolve so that it returns a 0x0 matrix when no soln exists
// (previously it threw an exception).
//
// Revision 1.43  2012/10/16 09:45:44  abbott
// Replaced RefRingElem by RingElem&.
//
// Revision 1.42  2012/10/03 15:25:05  abbott
// Replaced swap by assignment in DetByGauss; new impl of swap did not work
// in that instance (since one value was a temporary).
//
// Revision 1.41  2012/07/10 12:59:34  bigatti
// -- added two lines to keep compiler quiet
//
// Revision 1.40  2012/07/10 09:48:41  bigatti
// -- fixes of some naive errors
//
// Revision 1.39  2012/07/10 09:23:30  bigatti
// -- separated gauss code from LinSolveByGauss
// -- added LinKerByGauss
//
// Revision 1.38  2012/06/19 15:43:27  abbott
// Added division of a matrix by a scalar.
//
// Revision 1.37  2012/06/11 08:20:33  abbott
// Added multiplication on the right by a scalar.
//
// Revision 1.36  2012/06/10 22:57:31  abbott
// Added negation for matrices -- same as doing (-1)*M.
//
// Revision 1.35  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.34  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.33  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.32  2012/05/04 15:39:29  abbott
// Corrected a comment.
//
// Revision 1.31  2012/04/27 14:49:33  abbott
// Added LinSolve family (incl. LinSolveByGauss, LinSolveByHNF, LinSolveByModuleRepr).
//
// Revision 1.30  2012/04/16 09:21:19  abbott
// Added missing include directive.
//
// Revision 1.29  2012/04/13 16:24:35  abbott
// Added solve and SolveByGauss.
//
// Revision 1.28  2012/04/11 14:03:24  abbott
// Very minor change: slight improvement to readability.
//
// Revision 1.27  2012/03/16 14:42:46  bigatti
// -- fixed AdjointByDetOfMinors
// -- fixed adjoint (field + det(M)=0)
//
// Revision 1.26  2012/02/10 10:26:39  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.25  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.24  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.23  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.22  2011/05/13 16:47:20  abbott
// Added power fn for matrices: partial impl, cannot yet handle negative powers.
//
// Revision 1.21  2011/05/03 13:48:02  abbott
// Now using CanonicalHom inside DetByGauss.
// Cleaner and avoids a mysterious compiler warning.
//
// Revision 1.20  2011/03/22 16:44:19  bigatti
// -- fixed check in det
//
// Revision 1.19  2011/03/16 15:41:06  bigatti
// -- minor cleaning
//
// Revision 1.18  2011/03/10 11:25:46  bigatti
// -- now using long instead of size_t, and len(v) instead of v.size()
//
// Revision 1.17  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.16  2011/03/01 14:13:24  bigatti
// -- added f*M
//
// Revision 1.15  2011/02/28 14:08:49  bigatti
// -- added det3x3
// -- using apply mapping matrix (in DetByGauss)
//
// Revision 1.14  2011/02/10 15:27:06  bigatti
// -- commented #include vector  (included in MatrixArith.H)
//
// Revision 1.13  2011/02/09 16:48:27  bigatti
// -- added + and - for matrices
//
// Revision 1.12  2009/12/23 18:55:16  abbott
// Removed some useless comments.
//
// Revision 1.11  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.10  2009/06/25 16:59:42  abbott
// Minor improvement to some error messages (better coherence & comprehensibility).
//
// Revision 1.9  2008/07/09 16:09:11  abbott
// Removed pointless bogus function declarations.
//
// Revision 1.8  2008/04/22 14:42:03  abbott
// Two very minor changes.
//
// Revision 1.7  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.6  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.5  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.4  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.3  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/30 15:54:15  bigatti
// -- fixed index too big in RankByGauss(ConstMatrix M)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.8  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.5  2006/11/27 16:18:32  cocoa
// -- moved classes declarations from .H to .C (DenseMatrix, DiagMatrix,
//    FieldIdeal, SpecialMatrix)
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/08/17 09:39:07  cocoa
// -- added: elimination ordering matrix for non-homogeneous input
//
// Revision 1.2  2006/07/17 16:58:05  cocoa
// -- added: NewMatrixElim(size_t NumIndets, std::vector<size_t> IndetsToElim)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/02 14:39:20  cocoa
// -- Changed "not" into "!" becuase of M$Windoze (by M.Abshoff)
//
// Revision 1.7  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.6  2006/04/10 13:20:43  cocoa
// -- fixed buglets for Elimination orderings
//
// Revision 1.5  2006/04/05 16:45:29  cocoa
// -- added comment and CoCoA_ASSERT in NewPositiveMatrix
// -- added IsPositiveGrading
//
// Revision 1.4  2006/04/05 14:49:20  cocoa
// -- fixed: NewPositiveMatrix (tested and used in OrdvArith.C)
//
// Revision 1.3  2006/01/19 15:48:49  cocoa
// -- fixed RankByGauss by Stefan Kaspar
//
// Revision 1.2  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.1  2005/03/11 18:38:32  cocoa
// -- first import
//
