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


#include "CoCoA/utils_gmp.H"

namespace CoCoA
{

  int mpq_cmpabs(const mpq_t q1, const mpq_t q2)
  {
    if (mpq_sgn(q1) == 0) { if (mpq_sgn(q2) == 0) return 0; else return -1; }
    if (mpq_sgn(q2) == 0) return 1;

    if (mpz_cmp_ui(mpq_denref(q1),1) == 0 && mpz_cmp_ui(mpq_denref(q2),1) == 0)
      return mpz_cmpabs(mpq_numref(q1), mpq_numref(q2));

    const unsigned long logn1 = mpz_sizeinbase(mpq_numref(q1),2); // logn1-1 <= log(n1) < logn1
    const unsigned long logd1 = mpz_sizeinbase(mpq_denref(q1),2); //
    const unsigned long logn2 = mpz_sizeinbase(mpq_numref(q2),2); //
    const unsigned long logd2 = mpz_sizeinbase(mpq_denref(q2),2); //

    if (logn1+logd2-2 >= logn2+logd1) return 1;
    if (logn1+logd2 <= logn2+logd1-2) return -1;
    // IDEA: could compare more accurate logs of n1*d2 and n2*d1.
    mpz_t n1d2; mpz_init(n1d2);
    mpz_t n2d1; mpz_init(n2d1);
    mpz_mul(n1d2,mpq_numref(q1),mpq_denref(q2));
    mpz_mul(n2d1,mpq_numref(q2),mpq_denref(q1));
    const int ans = mpz_cmpabs(n1d2,n2d1);
    mpz_clear(n2d1);
    mpz_clear(n1d2);
    return ans;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/utils_gmp.C,v 1.1 2014/06/13 12:05:58 abbott Exp $
// $Log: utils_gmp.C,v $
// Revision 1.1  2014/06/13 12:05:58  abbott
// Summary: new GMP fn for CmpAbs of rationals
// Author: JAA
//
//
