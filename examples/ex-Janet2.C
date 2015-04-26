// Copyright (c) 2013  John Abbott -- Author: Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
#include <bitset>
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to using some functions which are using the Janet Basis.  \n";

const string LongDescription =
  " PrintNonMultVar, PrintMultVar, IsPommaretBasis, StandardRepresentation(f) \n";
//----------------------------------------------------------------------

void output(vector<RingElem> vec)
{
  for(vector<RingElem>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      cout << *iter << endl;
    }
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  ring Q = RingQQ();

  cout <<"/////////////////////////////////////Cyclic 4//////////////////////////////////////////" << endl;
  SparsePolyRing polyRingCyc5 = NewPolyRing(Q,4);
  vector<RingElem> x = indets(polyRingCyc5);
  vector<RingElem>cyclic5;
  cyclic5.push_back(x[0]+x[1]+x[2]+x[3]);
  cyclic5.push_back(x[0]*x[1] + x[2]*x[1] + x[2]*x[3] + x[0]*x[3]);
  cyclic5.push_back(x[1]*x[2]*x[0] + x[2]*x[3]*x[1] +  x[0]*x[3]*x[1] + x[0]*x[2]*x[3]);
  cyclic5.push_back(x[0]*x[1]*x[2]*x[3] - 1);
  // cyclic5.push_back(power(x[0], 7) - x[0] - 1);
  // cyclic5.push_back(power(x[0], 3) * x[1] - 4 * x[0] + 1);



  cout << "---------computing the Janet Base with a list of polynomials:------------------------" << endl;
  JBMill result(ExtendedJanetBasis(cyclic5));
  //----------ATTENTION we using ExtendedJanetBasis!!! (this returns the JBMill)
  //  result = ExtendedJanetBasis(cyclic5); 

  cout << "-----nonmultiplicative variables of the Janet-Base of cylic 4 (easy way)" << endl;
  result.myPrintNonMultVar();

  cout << endl << endl ;

  cout << "-----nonmultiplicative variables of the Janet-base of cyclic 4 (complicate way)" << endl;
  map<PPMonoidElem, vector<bool> > multVars = JBNonMultVar(result);
  for (map<PPMonoidElem, vector<bool> >::iterator i = multVars.begin(); i != multVars.end(); ++i)
  {
    cout << i->first << endl;
    long counter = 0;
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }

  cout << "-----multiplicative variables of the Janet-base of cyclic 4 (easy way)" << endl;
  result.myPrintMultVar();

  cout << endl << endl ;

  cout << "-----multiplicative variables of the Janet-base of cyclic 4 (complicate way)" << endl;
  map<PPMonoidElem,vector<bool> > nonMultVars = JBMultVar(result);
  for (map<PPMonoidElem,vector<bool> >::iterator i = nonMultVars.begin(); i != nonMultVars.end(); ++i)
  {
    cout << i->first << endl;
    long counter = 0;
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }

  cout << endl << endl;
  cout << "-----checking if the Janet-base of cyclic 4 is also a Pommaret basis" << endl;

  cout << "Janet-base of cyclic 4 is also a Pommaret basis: " << JBIsPommaretBasis(result) << endl;

  cout << endl << endl;
  RingElem f = 324*x[1]*x[1]*x[0]*x[0] + 234*power(x[0],4) - 123;

  cout << "involutive standard representation of the element f (easy way)" << endl;
  cout << "f = " << f << endl;
  JBOutputStandardRepresentation(result, f);

  cout << "involutive standard representation of the element f (complicate way)" << endl;
  pair<map<PPMonoidElem, RingElem>, RingElem> representation = JBStandardRepresentation(result, f);
  map<PPMonoidElem, RingElem> sRep = representation.first;
  RingElem res = representation.second;
  cout << "Involutive Standard Representation" << endl;
  cout << res << " =" << endl;
  cout << endl;
  for(map<PPMonoidElem, RingElem>::iterator iter(sRep.begin()); iter != sRep.end(); ++iter)
    {
      cout << iter->second <<  endl;
    }
  cout << "rest = " << res << endl;
  cout << "----------------------------------" << endl;


  cout << "the hilbert polynomial of cyclic 4 (with variable x[0])" << endl;
  
  ring Z = RingQQ();
  PolyRing hilbPolyRing = NewPolyRing(Z, symbols("t"));

  cout << "hilbert Polynomial = " << JBHilbertPol(result, indet(hilbPolyRing, 0)) << endl;
  // cout << "hilbert Polynomial = " << JBHilbertPol(result t) << endl;
  cout << "the hilbert function of cyclic 4 (with s = 10)" << endl;
  std::cout << "hilbfunc(" << 10 << ")=" << JBHilbertFunc(result, BigInt(10)) << endl;; //s must be of type ZZ 

  cout << "hilbert function version 2 with functional expression" << endl;
  JBHilbertFunc(result);

  cout << "the rational function of the hilbert series of cyclic 4" << endl;
  FractionField QR = NewFractionField(polyRingCyc5);
  RingHom embeddedHom(EmbeddingHom(QR));
  RingElem qx = embeddedHom(x[0]);
  
  cout << JBHilbertSeries(result, qx) << endl;

  cout << "generators of the first syzygy" << endl;
  FGModule firstSyz = JBSyzygy(result);

  vector<ModuleElem> gensFirstSyz = gens(firstSyz);
  for(vector<ModuleElem>::iterator iter = gensFirstSyz.begin(); iter != gensFirstSyz.end(); ++iter)
    {
      cout << *iter << endl;
    }

}

//----------------------------------------------------------------------
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
  return 1;
}
