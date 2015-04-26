//   Copyright (c)  2011  John Abbott

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/IntOperations.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

using namespace CoCoA;

// Handy macro for making assertions.
#define TEST_ASSERT(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)


template <typename T1, typename T2>
void test(const vector<T1>& v1, const vector<T2>& v2)
{
  const int n1 = len(v1);
  const int n2 = len(v2);

  for (int i=0; i < n1; ++i)
  {
    for (int j=0; j < n2; ++j)
    {
      const BigInt V1i(v1[i]);
      const BigInt V2j(v2[j]);
      TEST_ASSERT(sign(v1[i]) == sign(V1i));
      TEST_ASSERT(sign(v2[j]) == sign(V2j));
      TEST_ASSERT(cmp(v1[i],v2[j]) == cmp(V1i,V2j));
    }
  }
}


template <typename T>
vector<T>  MakeVec()
{
  vector<T> v;
  if (numeric_limits<T>::is_signed)
  {
    v.push_back(numeric_limits<T>::min());
    v.push_back(-1);
  }
  v.push_back(0);
  v.push_back(1);
  v.push_back(numeric_limits<T>::max()/2);
  v.push_back(numeric_limits<T>::max());
  return v;
}

void program()
{
  // This test checks the comparison fn "cmp" on built-in integral types.
  GlobalManager CoCoAFoundations;

  vector<signed char> SignedChar = MakeVec<signed char>();
  vector<unsigned char> UnsignedChar = MakeVec<unsigned char>();

  vector<signed short> SignedShort = MakeVec<signed short>();
  vector<unsigned short> UnsignedShort = MakeVec<unsigned short>();

  vector<signed int> SignedInt = MakeVec<signed int>();
  vector<unsigned int> UnsignedInt = MakeVec<unsigned int>();

  vector<signed long> SignedLong = MakeVec<signed long>();
  vector<unsigned long> UnsignedLong = MakeVec<unsigned long>();

  test(SignedChar, SignedChar);
  test(SignedChar, UnsignedChar);
  test(SignedChar, SignedShort);
  test(SignedChar, UnsignedShort);
  test(SignedChar, SignedInt);
  test(SignedChar, UnsignedInt);
  test(SignedChar, SignedLong);
  test(SignedChar, UnsignedLong);

  test(UnsignedChar, SignedChar);
  test(UnsignedChar, UnsignedChar);
  test(UnsignedChar, SignedShort);
  test(UnsignedChar, UnsignedShort);
  test(UnsignedChar, SignedInt);
  test(UnsignedChar, UnsignedInt);
  test(UnsignedChar, SignedLong);
  test(UnsignedChar, UnsignedLong);

  test(SignedShort, SignedChar);
  test(SignedShort, UnsignedChar);
  test(SignedShort, SignedShort);
  test(SignedShort, UnsignedShort);
  test(SignedShort, SignedInt);
  test(SignedShort, UnsignedInt);
  test(SignedShort, SignedLong);
  test(SignedShort, UnsignedLong);

  test(UnsignedShort, SignedChar);
  test(UnsignedShort, UnsignedChar);
  test(UnsignedShort, SignedShort);
  test(UnsignedShort, UnsignedShort);
  test(UnsignedShort, SignedInt);
  test(UnsignedShort, UnsignedInt);
  test(UnsignedShort, SignedLong);
  test(UnsignedShort, UnsignedLong);

  test(SignedInt, SignedChar);
  test(SignedInt, UnsignedChar);
  test(SignedInt, SignedShort);
  test(SignedInt, UnsignedShort);
  test(SignedInt, SignedInt);
  test(SignedInt, UnsignedInt);
  test(SignedInt, SignedLong);
  test(SignedInt, UnsignedLong);

  test(UnsignedInt, SignedChar);
  test(UnsignedInt, UnsignedChar);
  test(UnsignedInt, SignedShort);
  test(UnsignedInt, UnsignedShort);
  test(UnsignedInt, SignedInt);
  test(UnsignedInt, UnsignedInt);
  test(UnsignedInt, SignedLong);
  test(UnsignedInt, UnsignedLong);

  test(SignedLong, SignedChar);
  test(SignedLong, UnsignedChar);
  test(SignedLong, SignedShort);
  test(SignedLong, UnsignedShort);
  test(SignedLong, SignedInt);
  test(SignedLong, UnsignedInt);
  test(SignedLong, SignedLong);
  test(SignedLong, UnsignedLong);

  test(UnsignedLong, SignedChar);
  test(UnsignedLong, UnsignedChar);
  test(UnsignedLong, SignedShort);
  test(UnsignedLong, UnsignedShort);
  test(UnsignedLong, SignedInt);
  test(UnsignedLong, UnsignedInt);
  test(UnsignedLong, SignedLong);
  test(UnsignedLong, UnsignedLong);

}


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  BuildInfo::PrintAll(cerr);
  return 1;
}
