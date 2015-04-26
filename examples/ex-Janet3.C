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
  cout <<"///////////////A homogenous ideal with a Pommaret-Base//////////////////////////////////////////" << endl;
  SparsePolyRing polyRing6 = NewPolyRing(Q,6);
  vector<RingElem> x = indets(polyRing6);
  vector<RingElem> polys;

  polys.push_back(  power(x[0],2) + power(x[1],2) + power(x[2],2) );
  polys.push_back(  x[0]*x[3] + x[1]*x[4] + x[2]*x[5]);
  polys.push_back(  50*power(x[0],2) - 2*x[0]*x[4] + 14*x[0]*x[5] + power(x[1],2) - 14*x[1]*x[2] + 2*x[1]*x[3] + 49*power(x[2],2) - 14*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  29*power(x[0],2) - 10*x[0]*x[1] - 4*x[0]*x[2] - 4*x[0]*x[4] + 10*x[0]*x[5] + 5*power(x[1],2) - 20*x[1]*x[2] + 4*x[1]*x[3] - 2*x[1]*x[5] + 26*power(x[2],2) - 10*x[2]*x[3] + 2*x[2]*x[4] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  9*power(x[0],2) - 6*x[0]*x[5] + 9*power(x[2],2) + 6*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  9*power(x[0],2) + 6*x[0]*x[5] + 9*power(x[2],2) - 6*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );

  JBMill res7(ExtendedJanetBasis(polys));
  if(JBIsPommaretBasis(res7))
  {
    cout << "pom" << endl;
  }
  else
  {
    cout << "not pom" << endl;
  }
  if(JBIsHomogenous(res7))
  {
    std::cout << "homog" << std::endl;
  }
  else
  {
    std::cout << "not homog" << std::endl;
  }
  if(JBIsMonomialIdeal(res7))
  {
    std::cout << "monomial" << std::endl;
  }
  else
  {
    std::cout << "not monomial" << std::endl;
  }


  cout << "---------computing the Janet-Base/Pommaret-Base:------------------------" << endl;
  JBMill result(ExtendedJanetBasis(polys));
  //----------ATTENTION we using ExtendedJanetBasis!!! (this returns the JBMill)
  //  result = ExtendedJanetBasis(polys); 
  // cout << "---------printMultVar-------" << endl;
  // result.PrintMultVar();
  cout << "-----computing the dimension of P/I" << endl;
  cout << "dimension = " << JBDim(result) << endl << endl;

  cout << "-----computing the depth of the ideal" << endl;
  cout << "depth = " << JBDepth(result) << endl << endl;

  cout << "-----computing the projective dimension" << endl;
  cout << "projDim = " << JBProjDim(result) << endl;
 cout << "-----test if I is Cohen Macaulay" << endl;
  if (JBIsCohenMacaulay(result))
  {
    cout << "ideal is Cohen Macaulay" << endl;
  }
  else
  {
    cout << "ideal isn't Cohen Macaulay" << endl;
  }

  cout << "-----regular Sequence" << endl;
  std::vector<RingElem> regularSeq = JBRegularSequence(result);
  cout << "ok... not my fault..." << endl;
  for (std::vector<RingElem>::iterator i = regularSeq.begin(); i != regularSeq.end(); ++i)
  {
    if(i != regularSeq.begin())
    {
      cout << ", " << *i;
    }
    else
    {
      cout << *i;
    }
  }


  cout << "-----maxStronglyIndependentSet" << endl;
  std::vector<RingElem> stronglySet = JBMaxStronglyIndependentSet(result);
  for (std::vector<RingElem>::iterator i = stronglySet.begin(); i != stronglySet.end(); ++i)
  {
    if(i != stronglySet.begin())
    {
      cout << ", " << *i;
    }
    else
    {
      cout << *i;
    }
  }
  cout << endl;

  cout << "---------regularity" << endl;
  cout << "regularity = " << JBRegularity(result) << endl;

  cout << "----------satiety" << endl;
  cout << "satiety = " << JBSatiety(result) << endl;

  cout << "----------castelnuovo-mumford-regularity" << endl;
  cout << "cmr = " << JBCastelnuovoMumfordRegularity(result) << endl;

  cout << "----------saturation" << endl;
  std::vector<RingElem> satur = JBSaturation(result);
  long counter = 0;
  for (std::vector<RingElem>::iterator i = satur.begin(); i != satur.end(); ++i)
  {
    ++counter;
    cout << *i << endl;
  }

  cout << "----------classes stuff" << endl;

  RingElem clsRingElem = x[0]*x[1];

  cout << "cls of " << LPP(clsRingElem) << " = " << JBCls(result, LPP(clsRingElem)) <<  endl;

  cout << "minimal class of result = " << JBMinCls(result) << endl;

  cout << "all elements with cls 3 in the current basis" << endl;
  std::vector<RingElem> clsElems = JBElementsWithClass(result, BigInt(3));
  for (std::vector<RingElem>::iterator i = clsElems.begin(); i != clsElems.end(); ++i)
  {
    cout << *i << endl;
  }

  cout << endl << endl;


  cout << "degree of all pommaret classes" << endl;
  std::map<long, long> classes = result.myDegPommaretClasses();

  for(std::map<long,long>::iterator i = classes.begin(); i != classes.end(); ++i)
  {
    cout << "class = " << i->first << ", degree = " << i->second <<  endl;
  }
  cout << endl;


  cout << "extremalBettiNumbers" << endl;

  map<pair<long, long>, long> bettis = JBExtremalBettiNumbers(result);
  for (map<pair<long, long>, long>::iterator i = bettis.begin(); i != bettis.end(); ++i)
  {
    cout << "(" << i->first.first << ", " << i->first.second << ")" << "=" << i->second << endl;
  }

  cout << "degree of pommaret class at index 2 = " << JBDegPommaretClass(result, BigInt(1)) << endl;

  std::vector<RingElem> testCompDecomp1;
  testCompDecomp1.push_back(x[0] * x[0]); 
  testCompDecomp1.push_back(x[0] * x[1]); 
  testCompDecomp1.push_back(x[1] * x[1]); 
  testCompDecomp1.push_back(x[2] * x[3] * x[3]); 

  JBMill monResult(ExtendedJanetBasis(testCompDecomp1));



  cout << "complementary Decomposition" << endl;

  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > cD = monResult.myComplementaryDecomposition();
  for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = cD.begin(); i < cD.end(); ++i)
  {
    cout << "mon = " << i->first << endl;
    long counter = 0;
    cout << "vars = ";
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter == true)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }
  cout << endl;

  cout << "standard Pairs" << endl;

  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > sP = monResult.myStandardPairs();
  for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = sP.begin(); i < sP.end(); ++i)
  {
    cout << "mon = " << i->first << endl;
    long counter = 0;
    cout << "vars = ";
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter == true)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }
  cout << endl;


  cout << "Noether Normalization" << endl;
  SparsePolyRing NoetherRing = NewPolyRing(Q,4);
  x.clear();
  x = indets(NoetherRing);
  vector<RingElem> NoetherPolys;

  NoetherPolys.push_back(x[2]*x[2]*x[3] - x[0]*x[1]*x[2]*x[2]);
  NoetherPolys.push_back(x[1]*x[2]*x[3] - x[0]*x[3]*x[3]);
  NoetherPolys.push_back(x[2]*x[2]*x[3] - x[0]*x[1]*x[1]*x[2]*x[3]);

  JBMill NoetherMill(ExtendedJanetBasis(NoetherPolys));
  cout << "Janet Base" << endl;
  vector<RingElem> NoetherJanetBase = JBReturnJB(NoetherMill);
  for(std::vector<RingElem>::const_iterator i = NoetherJanetBase.begin(); i != NoetherJanetBase.end(); ++i)
  {
    cout << *i << endl;
  }
  cout << "Noether Normalization" << endl;
  pair<RingHom, std::vector<bool> >  NoetherResult = NoetherMill.myNoetherNormalization();
  cout << NoetherResult.first << endl;


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
