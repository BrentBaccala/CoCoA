//   Copyright (c)  2007,2010  John Abbott

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpF5.H"
#include "CoCoA/symbol.H"
#include "CoCoA/VectorOperations.H"

using namespace CoCoA;
using namespace std;

int FindReducerIndex(const PPMonoidElem& pp, const vector<RingElem>& v)
{
  const long nelems = len(v);
  for (long i=0; i < nelems; ++i)
    if (IsDivisible(pp, LPP(v[i])))
      return i;
  return -1;
}


RingElem NRLM(ConstRefRingElem f, const vector<RingElem>& v)
{
  if (IsZero(f)) return f;
  const SparsePolyRing P = owner(f);
  RingElem LMfOverLMvi(P), NRf(f);

  int i = FindReducerIndex(LPP(NRf), v);
  while (i != -1)
  {
    P->myDivLM(raw(LMfOverLMvi), raw(NRf), raw(v[i]));
    NRf -= LMfOverLMvi * v[i];
    if(IsZero(NRf)) return NRf;
    i = FindReducerIndex(LPP(NRf), v);
  }
  return NRf;
}


RingElem NormalRemainder(ConstRefRingElem f, const vector<RingElem>& v)
{
  if (IsZero(f)) return f;
  const SparsePolyRing P = owner(f);
  RingElem LMfOverLMvi(P), tmpNR(f), ansNR(P), LM(P);

  tmpNR = NRLM(f, v);
  while (!IsZero(tmpNR))
  {
    P->myMoveLM(raw(LM), raw(tmpNR));
    P->myAppendClear(raw(ansNR), raw(LM)); // now LM is 0
    tmpNR = NRLM(tmpNR, v);
  }
  return ansNR;
}

void program()
{
  GlobalManager CoCoAFoundations;
  SparsePolyRing Qx = NewPolyRing(RingQQ(), SymbolRange("x",0,5)); //Q[x0..x5]
  const vector<RingElem>& x = indets(Qx);

  vector<RingElem> g;
  g.push_back(power(x[1],5) - power(x[3],5));
  g.push_back(power(x[1],3) - x[1]*x[2]*x[3] + x[2]*x[2]*x[3]);
  ideal I = ideal(Qx, g);
  g.clear();

  cout << "gens(I)   =   " << gens(I) << endl;
  vector<RingElem> GB = TidyGens(I);  // it is the Groebner Basis of I
  cout << "TidyGens(I) = " << GB << endl;
  cout << "When I is an ideal of polynomials TidyGens returns its Groebner Basis." << endl << endl;

  F5(g, I);
  cout << "F5(I) = " << g << endl;

  RingElem f(power(x[1],12) + power(x[2],6) + power(x[3],6));
  cout << "-- f = " << f << endl;
  cout << "NormalRemainder(f, gens(I)) = " <<  NormalRemainder(f, gens(I)) << endl;
  cout << "NormalRemainder(f, GB)      = " <<  NormalRemainder(f, GB) << endl;
  cout << "NormalRemainder(f, g)       = " <<  NormalRemainder(f, GB) << endl;
  cout << "NF(f, I)                    = " <<  NF(f, I) << endl;
  cout << endl;
  
  RingElem h(power(x[1],3)*g[0] + power(x[2],4)*g[1]);
  cout << "-- h = " << h << endl;
  cout << "NormalRemainder(h, gens(I)) = " <<  NormalRemainder(h, gens(I)) << endl;
  cout << "NormalRemainder(h, GB)      = " <<  NormalRemainder(h, GB) << endl;
  cout << "NormalRemainder(h, g)       = " <<  NormalRemainder(h, GB) << endl;
  cout << "NF(h, I)                    = " <<  NF(h, I) << endl;
}


int main()
{
  try
  {
    program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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

