#include "CoCoA/TmpJBEnv.H"


using namespace CoCoA;

std::vector<RingElem> JBMill::myReturnJB() const
{
  std::vector<RingElem> output;
    for(std::list<JanetTriple>::const_iterator it = mySetT.begin(); it != mySetT.end(); ++it)
    {
      output.push_back(it->myGetPol());
    }
  return output;
}


std::vector<RingElem> JBMill::myReturnGB() const
{
  std::vector<RingElem> output;
  // an element is part of the reduced GB if anc = lpp
  for(std::list<JanetTriple>::const_iterator it = mySetT.begin(); it != mySetT.end(); ++it)
  {
    if(LPP(it->myGetPol()) == it->myGetAnc())
    {
      output.push_back(it->myGetPol());
    }
  }
  return output;
}


void JBMill::myPrintMultVar() const
{
  std::cout << "computation of multiplicative variables:" << std::endl;
  std::map<PPMonoidElem,std::vector<bool> > MultVars = myComputeMultVars();
  myOutputVar(MultVars, true);
}

void JBMill::myPrintNonMultVar() const
{
  std::cout << "computation of nonmultiplicative variables:" << std::endl;
  std::map<PPMonoidElem,std::vector<bool> > MultVars = myComputeMultVars();
  myOutputVar(MultVars, false);
}


void JBMill::myOutputVar(std::map<PPMonoidElem,std::vector<bool> > MultVars, bool OutputMultVar) const
{
  for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter(MultVars.begin()); iter != MultVars.end(); ++iter)
  {
    std::cout << "LT(pol) = " << (*iter).first << std::endl;
    if(OutputMultVar)
    {
      std::cout << "multiplicative variables: ";
    }
    else
    {
      std::cout << "nonmultiplicative variables: ";
    }
    bool begin(true);
    const long VecSize = len((*iter).second);
    // printing index of every iter, which is equal to OutputMultVar
    for(long VecIter(0); VecIter != VecSize; ++VecIter)
    {
      if((*iter).second[VecIter] == OutputMultVar)
      {
        if(!begin)
        {
          std::cout << ", ";
        }
        begin = false;
        std::cout << VecIter;
      }
    }
    std::cout << std::endl;
    std::cout << "---------------------------------------" << std::endl;
  }
}


std::map<PPMonoidElem,std::vector<bool> > JBMill::myComputeMultVars() const
{
  //initialization
  JanetIterator iter(myJTree);
  std::vector<int> CurrentNonMultVars;
  long CountIndets = NumIndets(myOptions.myRingOptions.myPolyRing);
  std::map<PPMonoidElem,std::vector<bool> > MultVars;
  for(std::list<JanetTriple>::const_iterator VecIter(mySetT.begin()); VecIter != mySetT.end(); ++VecIter)
  {
    MultVars.insert(std::pair<PPMonoidElem, std::vector<bool> >(LPP(VecIter->myGetPol()), std::vector<bool>(CountIndets, true)));
  }
  //computation
  myRekComputeMultVar(MultVars, iter, CurrentNonMultVars, 0);
  //output
  return MultVars;
}

std::map<PPMonoidElem, std::vector<bool> > JBMill::myComputeNonMultVars() const
{
  std::map<PPMonoidElem, std::vector<bool> > MultVars = myComputeMultVars();
  //reversing multiplicative variables -> nonMultVars
  for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVars.begin(); i != MultVars.end(); ++i)
  {
    i->second = myReverseBoolVec(i->second);  
  }
  return MultVars;
}

std::vector<bool> JBMill::myNonMultVarsOf(RingElem elem) const
{
  std::map<PPMonoidElem, std::vector<bool> > NonMultVars = myComputeNonMultVars();
  std::map<PPMonoidElem, std::vector<bool> >::iterator i(NonMultVars.find(LPP(elem)));
  if(i != NonMultVars.end())
  {
    return i->second;
  }
  return std::vector<bool>(NumIndets(myOptions.myRingOptions.myPolyRing), true);
}


std::vector< std::pair<RingElem, std::vector<bool> > > JBMill::myComputeNonMultVarsWithRingElem() const
{
  std::map<PPMonoidElem, std::vector<bool> > MultVarsWithPPM = myComputeNonMultVars();
  std::vector< std::pair<RingElem, std::vector<bool> > > result;
  if(IsOne(LPP(mySetT.begin()->myGetPol())))
  {
    result.push_back(std::pair<RingElem, std::vector<bool> >(mySetT.begin()->myGetPol(), (MultVarsWithPPM.begin())->second));
    return result;
  }
  for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVarsWithPPM.begin(); i != MultVarsWithPPM.end(); ++i)
  {
    JanetTriple* gPtr(myJTree.myJDivisor(i->first));
    result.push_back(std::pair<RingElem, std::vector<bool> >(gPtr->myGetPol(), i->second));
  }
  return result;
}

std::vector< std::pair<RingElem, std::vector<bool> > > JBMill::myComputeMultVarsWithRingElem() const
{
  std::map<PPMonoidElem, std::vector<bool> > MultVarsWithPPM = myComputeMultVars();
  std::vector< std::pair<RingElem, std::vector<bool> > > result;
  if(IsOne(LPP(mySetT.begin()->myGetPol())))
  {
    result.push_back(std::pair<RingElem, std::vector<bool> >(mySetT.begin()->myGetPol(), (MultVarsWithPPM.begin())->second));
    return result;
  }
  for(std::map<PPMonoidElem, std::vector<bool> >::iterator i = MultVarsWithPPM.begin(); i != MultVarsWithPPM.end(); ++i)
  {
    JanetTriple* gPtr(myJTree.myJDivisor(i->first));
    result.push_back(std::pair<RingElem, std::vector<bool> >(gPtr->myGetPol(), i->second));
  }
  return result;
}



void JBMill::myRekComputeMultVar(std::map<PPMonoidElem, std::vector<bool> >& MultVars, JanetIterator iter, std::vector<int> CurrentNonMultVars, int CurVar) const
{
  do //until highest node in degree direction
  {
    // copy the current non mult vars
    std::vector<int> CopyVars = CurrentNonMultVars;
    // if this isn't the highest degree in current var add var to CopyVars
    if(iter.myDisNextDeg() != 0)
    {
      CopyVars.push_back(CurVar);
    }
    // if there is a node in variable-direction call myRekComputeMultVar
    if(iter.myDisNextVar())
    {
      int distance(iter.myDisNextVar());
      JanetIterator TmpIter(iter);
      TmpIter.myNextVar();
      myRekComputeMultVar(MultVars, TmpIter, CopyVars, CurVar + distance);
    }
  }
  while(iter.myNextDeg());

  // highest degree node for current variable -> this variable is multiplicative
  if(!(iter.myDisNextVar()))
  {
    //everything is multiplicative except the CurrentNonMultVars
    for(std::vector<int>::iterator VecIter(CurrentNonMultVars.begin()); VecIter != CurrentNonMultVars.end(); ++VecIter)
    {
      MultVars[LPP(iter.myGetPol())][*VecIter]= false;
    }
  }
}


std::vector<bool> JBMill::myPommaretMultVar(PPMonoidElem pp, long CountIndets) const
{
  std::vector<bool> MultVars(CountIndets, false);
  for(long iter(CountIndets - 1); iter != -1; --iter)
  {
    MultVars[iter] = true;
    if(exponent(pp, iter) != 0)
    {
      break;
    }
  }
  return MultVars;
}

bool JBMill::IamPommaretBasisDecide()
{
  if(IsUncertain3(IamPommaret))
  {
    IamPommaret = true;
    std::map<PPMonoidElem,std::vector<bool> > MultVars = myComputeMultVars();
    for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter(MultVars.begin()); iter != MultVars.end(); ++iter)
    {
      std::vector<bool> IPBPommaretMultVar(myPommaretMultVar((*iter).first, NumIndets(myOptions.myRingOptions.myPolyRing)));
      const long VecSize = len((*iter).second);
      for(long VecIter(0); VecIter != VecSize; ++VecIter)
      {
        if((*iter).second[VecIter] != IPBPommaretMultVar[VecIter])
        {
          IamPommaret = false;
          break;
        }
      }
    } 
  }
  return IsTrue3(IamPommaret);
}

bool JBMill::IamPommaretBasis() const
{
  if(IsUncertain3(IamPommaret))
  {
    CoCoA_ERROR("it is uncertain if the ideal is a pommaret basis", "IamPommaretBasis");
  }
  return IsTrue3(IamPommaret);
}

bool JBMill::IamHomogenousDecide()
{
  if(IsUncertain3(IamHomog))
  {
    IamHomog = true;
    std::vector<RingElem> gb(myReturnGB());
    for(std::vector<RingElem>::iterator iter = gb.begin(); iter != gb.end(); ++iter)
    {
      if(!IsHomog(*iter))
      {
        IamHomog = false;
        break;
      }
    }
  }
  return IsTrue3(IamHomog);
}

bool JBMill::IamHomogenous() const
{
  if(IsUncertain3(IamHomog))
  {
    CoCoA_ERROR("it is uncertain if the ideal is homogenous", "IamHomogenous");
  }
  return IsTrue3(IamHomog);
}



bool JBMill::IamMonomialIdealDecide()
{
  if(IsUncertain3(IamMonomial))
  {
    IamMonomial = true;
    std::vector<RingElem> jb(myReturnJB());
    for(std::vector<RingElem>::iterator iter = jb.begin(); iter != jb.end(); ++iter)
    {
      if(!IsMonomial(*iter))
      {
        IamMonomial = false;
        break;
      }
    }
  }
  return IsTrue3(IamMonomial);
}

bool JBMill::IamMonomialIdeal() const
{
  if(IsUncertain3(IamMonomial))
  {
    CoCoA_ERROR("it is uncertain if the ideal is a monomial ideal", "IamMonomialIdeal");
  }
  return IsTrue3(IamMonomial);
}


//TODO: Throw error if r has the wrong type!!!
std::pair<std::map<PPMonoidElem, RingElem>, RingElem> JBMill::myStandardRepresentation(RingElem r) const
{
  std::map<PPMonoidElem, RingElem> StandardRep;
  //initializing reduction history -> nothing is reduced every factor must be zero
  for(std::list<JanetTriple>::const_iterator iter(mySetT.begin()); iter != mySetT.end(); ++iter)
  {
    StandardRep.insert(std::pair<PPMonoidElem, RingElem>(LPP(iter->myGetPol()), zero(myOptions.myRingOptions.myPolyRing)));
  }
  RingElem h(r);
  RingElem res(myOptions.myRingOptions.myPolyRing);
  //testing if basis is one
///  if(IsOne(LPP(mySetT.begin()->myGetPol())))
  if(myOptions.myRingOptions.myPPMValue->myIsOne(raw(myOptions.myRingOptions.myPolyRing->myLPP(raw(mySetT.begin()->myGetPol())))))
  {
///    StandardRep[LPP(mySetT.begin()->myGetPol())] = h;
    StandardRep[myOptions.myRingOptions.myPolyRing->myLPP(raw(mySetT.begin()->myGetPol()))] = h;
    h = 0;
  }
  //perform involutive full reduction
///  while(!IsZero(h))
  while(!myOptions.myRingOptions.myPolyRing->myIsZero(raw(h)))
  {
///    JanetTriple* gPtr(myJTree.myJDivisor(LPP(h)));
    JanetTriple* gPtr(myJTree.myJDivisor(myOptions.myRingOptions.myPolyRing->myLPP(raw(h))));
    if(gPtr == 0)
    {
///      res += monomial(myOptions.myRingOptions.myPolyRing,LC(h),LPP(h));
///      h -= monomial(myOptions.myRingOptions.myPolyRing,LC(h),LPP(h));
      myOptions.myRingOptions.myPolyRing->myAddMul(raw(res), raw(h), raw(one(myOptions.myRingOptions.myPolyRing)));
      myOptions.myRingOptions.myPolyRing->myDeleteLM(raw(h));
    }
    else
    {
      //if we can reduce something add factor to the according element in the map
      RingElem factor(myOptions.myRingOptions.myPolyRing);
      myOptions.myRingOptions.myPolyRing->myDivLM(raw(factor), raw(h), raw(gPtr->myGetPol()));
      myOptions.myRingOptions.myPolyRing->myAddMul(raw(h), raw(-factor), raw(gPtr->myGetPol()));
      myOptions.myRingOptions.myPolyRing->myAdd(raw(factor), raw(factor), raw(StandardRep.find(LPP(gPtr->myGetPol()))->second));
      StandardRep.erase(LPP(gPtr->myGetPol()));
      StandardRep.insert(std::pair<PPMonoidElem,RingElem>(LPP(gPtr->myGetPol()), factor));
    }
  }
  // store in std::pair, first part is the reduction history, second part is the rest
  std::pair<std::map<PPMonoidElem, RingElem>, RingElem> pair(StandardRep, res);
  return pair;
}

std::vector<std::pair<RingElem, RingElem> > JBMill::myStandardRepresentationWithoutRest(RingElem r) const
{
  std::vector<std::pair<RingElem, RingElem> > result; 
  std::map<PPMonoidElem, RingElem> StandardRep(myStandardRepresentation(r).first);
  for(std::map<PPMonoidElem, RingElem>::iterator i = StandardRep.begin(); i != StandardRep.end(); ++i)
  {
    JanetTriple* gPtr(myJTree.myJDivisor(i->first));
    result.push_back(std::pair<RingElem, RingElem>(gPtr->myGetPol(), i->second));
  }
  return result;
}

std::vector<RingElem> JBMill::myStandardRepresentationWithoutRestShort(RingElem r) const
{
  std::vector<RingElem> result; 
  std::map<PPMonoidElem, RingElem> StandardRep(myStandardRepresentation(r).first);
  for(std::map<PPMonoidElem, RingElem>::iterator i = StandardRep.begin(); i != StandardRep.end(); ++i)
  {
    result.push_back(i->second);
  }
  return result;
}


void JBMill::myOutputStandardRepresentation(RingElem r) const
{
  std::pair<std::map<PPMonoidElem, RingElem>, RingElem> representation = myStandardRepresentation(r);
  std::map<PPMonoidElem, RingElem> StandardRep = representation.first;
  RingElem res = representation.second;
  std::cout << "Involutive Standard Representation" << std::endl;
  std::cout << r << " =" << std::endl;
  std::cout << std::endl;
  for(std::map<PPMonoidElem, RingElem>::iterator iter(StandardRep.begin()); iter != StandardRep.end(); ++iter)
  {
    std::cout << iter->second <<  std::endl;
  }
  std::cout << "rest = "<< res << std::endl;
  std::cout << "----------------------------------" << std::endl;
}


RingElem JBMill::myJNormalForm(const RingElem& elem) const
{
  RingElem h(elem);
  RingElem res(zero(myOptions.myRingOptions.myPolyRing));
  if(IsOne(LPP(mySetT.begin()->myGetPol())))
  {
    h = 0;
  }
  while(!IsZero(h))
  {
    JanetTriple* gPtr(myJTree.myJDivisor(myOptions.myRingOptions.myPolyRing->myLPP(raw(h))));
    if(gPtr == 0)
    {
      res += monomial(myOptions.myRingOptions.myPolyRing, LC(h), LPP(h));
      h -= monomial(myOptions.myRingOptions.myPolyRing, LC(h), LPP(h));
    }
    else
    {
      myOptions.myRingOptions.myPolyRing->myReductionStep(raw(h), raw(gPtr->myGetPol()));
    }
  }
  if(IsFractionField(CoeffRing(myOptions.myRingOptions.myPolyRing)) && !IsZero(h))
  {
    myOptions.myRingOptions.myPolyRing->myDivByCoeff(raw(res),raw(LC(res)));
    res = ClearDenom(res);
  }
  return res;
}

void JBMill::myDegreeTQ()
{
  JBSets sets(myOptions.myRingOptions.myPolyRing, myOptions.myCriteria);
  //initialization of janet-triples
  for(std::vector<RingElem>::iterator iter = myInput.begin(); iter != myInput.end(); ++iter)
  {
    JanetTriple triple(*iter, myOptions.myRingOptions.myPolyRing->myLPP(raw(*iter)));
    sets.myInsertSetQ(triple);
  }
  std::list<JanetTriple>::iterator iter(sets.myMinimizeAndInsertSetT(myJTree));
  //exit if there is a constant input
  if(IsOne(myOptions.myRingOptions.myPolyRing->myLPP(raw(iter->myGetPol()))))
  {
    myJTree.myDelete();
    mySetT.push_back(JanetTriple(one(myOptions.myRingOptions.myPolyRing)));
    return ;
  }
  // adding first element to janet tree
  JanetTree tree(sets.myOneLeafTree(iter, 0, one(myOptions.myRingOptions.myPPMValue)));
  myJTree.myAddAtBegin(tree);
  //main loop
  while(sets.myHeadReduceSetQ(myJTree))
  {
    long SizeT(sets.mySizeSetT());
    //insert element to setT and minimize setT
    std::list<JanetTriple>::iterator iter(sets.myMinimizeAndInsertSetT(myJTree));
    // exit if we get a constant element
    if(IsOne(myOptions.myRingOptions.myPolyRing->myLPP(raw(iter->myGetPol()))))
    {
      myJTree.myDelete();
      mySetT.push_back(JanetTriple(one(myOptions.myRingOptions.myPolyRing)));
      return ;
    }
    //if we 'lost' an element in setT rebuild janet tree
    //otherwise add new element to setT and tailreduce setT
    if(sets.mySizeSetT() <= SizeT)
    {
      myJTree.myDelete();
      for(JBSets::MultisetIterator IterT(sets.myBeginSetT()); IterT != sets.myEndSetT(); ++IterT)
      {
        std::list<JanetTriple>::iterator ListIter(*IterT);
        sets.myJTailNormalForm(myJTree, ListIter);
        sets.myJInsert(ListIter, myJTree);
      }
    } 
    else
    {
      sets.myJInsert(iter, myJTree);  
      sets.myTailReduceSetT(myJTree, iter);
    }  
  }
  //copy everything to mySetT in JBMill
  myJTree.myDelete();
  for(JBSets::MultisetIterator IterT(sets.myBeginSetT()); IterT != sets.myEndSetT(); ++IterT)
  {
    mySetT.push_back(*(*IterT));
  }
  //rebuilding JanetTree (Iterators have to point to mySetT and not to sets.mySetT)
  for(std::list<JanetTriple>::iterator iter(mySetT.begin()); iter != mySetT.end(); ++iter)
  {
    sets.myJInsertWithoutProlong(myJTree, iter);
  }
}

const std::vector<RingElem> CoCoA::JanetBasis(const std::vector<RingElem>& PolyList, const std::bitset<3> crit, ResOutputFlag res, StrategyFlag algorithm)
{
  if(PolyList.empty())
  {
    CoCoA_ERROR("Empty input", "JBMill");
  }

  JBEnv enviroment(owner(*(PolyList.begin())));
  JBFlag flag(enviroment);
  if(crit.test(0))
  {
    flag.myCriteria.set(0);
  }
  if(crit.test(1))
  {
    flag.myCriteria.set(1);
  }
  if(crit.test(2))
  {
    flag.myCriteria.set(2);
  }
  JBMill mill(PolyList, flag);
  switch(algorithm)
  {
    case  TQDegree: mill.myDegreeTQ();
          break;
    case  TQBlockHigh: mill.myBlockTQ(true);
          break;
    case  TQBlockLow: mill.myBlockTQ(false);
          break;
  }
  if(res == JB)
  {
    return mill.myReturnJB();
  }
  return mill.myReturnGB();
}

JBMill CoCoA::ExtendedJanetBasis(const std::vector<RingElem>& PolyList, const std::bitset<3> crit, StrategyFlag algorithm)
{
  if(PolyList.empty())
  {
    CoCoA_ERROR("Empty input", "JBMill");
  }

  JBEnv enviroment(owner(*(PolyList.begin())));
  JBFlag flag(enviroment);
  if(crit.test(0))
  {
    flag.myCriteria.set(0);
  }
  if(crit.test(1))
  {
    flag.myCriteria.set(1);
  }
  if(crit.test(2))
  {
    flag.myCriteria.set(2);
  }
  JBMill mill(PolyList, flag);
  switch(algorithm)
  {
    case  TQDegree: mill.myDegreeTQ();
          break;
    case  TQBlockHigh: mill.myBlockTQ(true);
          break;
    case  TQBlockLow: mill.myBlockTQ(false);
          break;
  }
  mill.IamPommaretBasisDecide();
  mill.IamHomogenousDecide();
  mill.IamMonomialIdealDecide();
  return mill;
}


void JBMill::myBlockTQ(bool UpdateHigh)
{
  //initialization
  JBSets sets(myOptions.myRingOptions.myPolyRing, myOptions.myCriteria, UpdateHigh);
  for(std::vector<RingElem>::iterator iter = myInput.begin(); iter != myInput.end(); ++iter)
  {
    if(IsZero(*iter))
    {
      continue;
    }
    if(IsFractionField(CoeffRing(myOptions.myRingOptions.myPolyRing)))
    {
      myOptions.myRingOptions.myPolyRing->myDivByCoeff(raw(*iter),raw(LC(*iter)));
      *iter = ClearDenom(*iter);
    }
    JanetTriple triple(*iter, myOptions.myRingOptions.myPolyRing->myLPP(raw(*iter)));
    sets.myInsertSetQ(triple);
  }
  std::list<JanetTriple>::iterator iter(sets.myMinimizeAndInsertSetT(myJTree));
  //return if constant input
  if(IsOne(myOptions.myRingOptions.myPolyRing->myLPP(raw(iter->myGetPol()))))
  {
    myJTree.myDelete();
    mySetT.push_back(JanetTriple(one(myOptions.myRingOptions.myPolyRing)));
    return ;
  }
  JanetTree tree(sets.myOneLeafTree(iter, 0, one(myOptions.myRingOptions.myPPMValue)));
  myJTree.myAddAtBegin(tree);
  //main loop
  while(!sets.IamEmptySetQ())
  {
    //construct mySetP
    sets.myMoveFromQtoP(myJTree);
    // update mySetP
    bool ModifySetT = sets.IamJUpdate();
    // if mySetP empty continue
    if(sets.myBeginSetP() == sets.myEndSetP())
    {
      continue;
    }
    if(ModifySetT)
    {
      //recompute janettree
      long SizeT(sets.mySizeSetT());
      sets.myMinimizeSetT();
      if(sets.mySizeSetT() < SizeT)
      {
        myJTree.myDelete();
        for(JBSets::MultisetIterator IterT(sets.myBeginSetT()); IterT != sets.myEndSetT(); ++IterT)
        {
          std::list<JanetTriple>::iterator ListIter(*IterT);
          sets.myJInsert(ListIter, myJTree);
        }
      }
      sets.myTailReduceSetTAll(myJTree);
    }
    JBSets::MultisetIterator IterP = sets.myEndSetP();
    //constant input -> retun
    if(IsOne((*(sets.myBeginSetP()))->myGetPol()))
    {
      myJTree.myDelete();
      mySetT.clear();
      mySetT.push_back(JanetTriple(one(myOptions.myRingOptions.myPolyRing)));
      return ;
    }
    // insert elements from mySetP into janet tree 
    do
    {
      --IterP;
      std::list<JanetTriple>::iterator ListIter(*IterP);
      sets.myJTailNormalForm(myJTree, ListIter);
      sets.myJInsert(ListIter, myJTree);
    }while(IterP != sets.myBeginSetP());
    //tail reduce set T if necessary
    if(!ModifySetT)
    {
      sets.myTailReduceSetTAll(myJTree);
    }
    sets.myInsertSetPInSetT(); 
  }
  //output
  myJTree.myDelete();
  mySetT.clear(); 
  //copy elements from sets.mySetT to mySetT
  for(JBSets::MultisetIterator IterT(sets.myBeginSetT()); IterT != sets.myEndSetT(); ++IterT)
  {
    mySetT.push_back(*(*IterT));
  }
  //recompute janet tree
  for(std::list<JanetTriple>::iterator iter(mySetT.begin()); iter != mySetT.end(); ++iter)
  {
    sets.myJInsertWithoutProlong(myJTree, iter);
  }
}

RingElem JBMill::myBinLike(PolyRing ring, RingElem PolAbove, long IntBelow) const
{
  RingElem result = one(ring);
  for (long j = 1; j <= IntBelow; ++j)
  {
    result *= (PolAbove + 1 - j)/j;
  }
  return result;
}

long JBMill::myCountTrues(std::vector<bool> vec) const
{
  long result(0);
  for(std::vector<bool>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
  {
    if(*iter)
    {  
      ++result;
    }
  }
  return result;
}

RingElem JBMill::myHilbertPol(RingElem s) const
{
  PolyRing P = owner(s);
  RingElem result(P);

  result = myBinLike(P, s + NumIndets(myOptions.myRingOptions.myPolyRing) - 1, NumIndets(myOptions.myRingOptions.myPolyRing) - 1);

  std::map<PPMonoidElem, std::vector<bool> > MultVar(myComputeMultVars());

  for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
  {
    long deg = StdDeg(iter->first);
    long NumMultVars = myCountTrues(iter->second);
    result -= myBinLike(P, s - deg + NumMultVars - 1, NumMultVars - 1);
  }
  return result;
}


void JBMill::myHilbertFunc() const
{
  ring Q = RingQQ();
  PolyRing HilbPolyRing = NewPolyRing(Q, symbols("t"));
  RingElem FuncExpression(myHilbertPol(indet(HilbPolyRing, 0)));
  long maxDeg = StdDeg(FuncExpression);
  for (int i = 0; i < maxDeg; ++i)
  {
    std::cout << "H(" << i << ") = " << myHilbertFunc(BigInt(i)) << std::endl;
  }
  std::cout << "H(t) = " << FuncExpression << "    for t >= " << maxDeg << std::endl;
}

BigInt JBMill::myHilbertFunc(BigInt m) const
{
  BigInt k(NumIndets(myOptions.myRingOptions.myPolyRing));
  BigInt r = m + k - 1;
  BigInt s = k - 1;
  BigInt res;
  if((r < s) || (r < 0) || (s < 0))
  {
    if((r == -1) && (s == -1))
    {
      res = 1;
    }
  }
  else
  {
    res = binomial(r, s);
  }


  std::map<PPMonoidElem, std::vector<bool> > MultVar(myComputeMultVars());
  for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
  {
    BigInt deg = BigInt(StdDeg(iter->first));
    BigInt NumMultVars = BigInt(myCountTrues(iter->second));
    BigInt r = m - deg + NumMultVars - 1;
    BigInt s = NumMultVars - 1;
    if((r < s) || (r < 0) || (s < 0))
    {
      if((r == -1) && (s == -1))
      {
        res -= 1;
      }
    }
    else
    {
      res -= binomial(r, s);
    }
  }
  return res;
}

RingElem JBMill::myHilbertSeries(RingElem s) const
{
  RingElem res(owner(s));
  if(IsFractionField(owner(s)))
  {
    std::map<PPMonoidElem, std::vector<bool> > MultVar(myComputeMultVars());
    for(std::map<PPMonoidElem,std::vector<bool> >::iterator iter = MultVar.begin(); iter != MultVar.end(); ++iter)
    {
      long deg = StdDeg(iter->first);
      long NumMultVars = myCountTrues(iter->second);
      RingElem numerator(power(s,deg));
      RingElem denominator(power(1 - s, NumMultVars));
      res += numerator/denominator;
    }
  }
  else
  {
    CoCoA_ERROR(ERR::NotElemFrF, "s must be in a Fraction Field!!");
  }
  return res;
}

FreeModule JBMill::myMakeNewFreeModuleForSyz() const
{
  std::vector<RingElem> jb = myReturnJB();
  return NewFreeModuleForSyz(jb);
}

FGModule JBMill::mySyzygy() const
{
  std::vector<ModuleElem> GenSyz;
  FreeModule FModule = myMakeNewFreeModuleForSyz(); //generates the module
  std::map<PPMonoidElem, std::vector<bool> > listMultVars = myComputeMultVars(); //computes the multiplicative variables for the janet basis
  std::vector<ModuleElem> GenFModule = gens(FModule); //generators of the module
  std::vector<ModuleElem>::iterator GenFModuleIter = GenFModule.begin(); 
  for(std::list<JanetTriple>::const_iterator ListIter = mySetT.begin(); ListIter != mySetT.end(); ++ListIter) //iterates over all elements in the janet basis
  {
    std::vector<bool> MultVars = (listMultVars.find(LPP(ListIter->myGetPol())))->second; //multiplicative variables of the current element
    long VarPos(0); //variable position, used in the for-loop below
    for(std::vector<bool>::iterator MultVar = MultVars.begin(); MultVar != MultVars.end(); ++MultVar) //iterates over all variables
    {
      if(!(*MultVar)) //tests if the variables is nonmultiplicative
      {
        std::map<PPMonoidElem, RingElem> StandardRep = myStandardRepresentation(indet(myOptions.myRingOptions.myPolyRing, VarPos) * ListIter->myGetPol()).first; //computes the standardrepresentation of nm-var * element
        ModuleElem s = indet(myOptions.myRingOptions.myPolyRing, VarPos) * (*GenFModuleIter);
        std::vector<ModuleElem>::const_iterator StandardRepModuleIter = GenFModule.begin();
        for(std::map<PPMonoidElem, RingElem>::const_iterator StandardRepIter = StandardRep.begin(); StandardRepIter != StandardRep.end(); ++StandardRepIter)
        {
          s = s - (StandardRepIter->second) * (*StandardRepModuleIter);
          ++StandardRepModuleIter;
        }
        GenSyz.push_back(s);
      }
      ++VarPos;
    }
    ++GenFModuleIter;
  }
  return submodule(FModule, GenSyz);
}

long JBMill::myDim() const
{
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = myComplementaryDecompositionLeadingIdeal();
  long deg(0);
  for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = CompDecomp.begin(); i != CompDecomp.end(); ++i)
  {
    if(deg < myCountTrues(i->second))
    {
      deg = myCountTrues(i->second);
    }
  }
  return deg;
}

long JBMill::myDepth() const
{
  long res(0);
  if(IamPommaretBasis() && IamHomogenous())
  {
    std::map<PPMonoidElem, std::vector<bool> > MultVars = myComputeMultVars();
    long depth(NumIndets(myOptions.myRingOptions.myPolyRing));
    for(std::map<PPMonoidElem, std::vector<bool> >::iterator iter = MultVars.begin(); iter != MultVars.end(); ++iter)
    {
      long CurrentMultVars(myCountTrues(iter->second));
      if(depth > CurrentMultVars)
      {
        depth = CurrentMultVars;
      }
    }
    res = depth;
    // --res;
  }
  else
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis or the ideal isn't homogenous", "myDepth");
  }
  return res;
}

long JBMill::myProjDim() const
{
  return NumIndets(myOptions.myRingOptions.myPolyRing) - myDepth(); //depth throws the error if not pommaret or homogenous
}

std::vector<RingElem> JBMill::myRegSeq() const//mod I!!!
{
  std::vector<RingElem> x = indets(myOptions.myRingOptions.myPolyRing);
  std::vector<RingElem> res;
  for(long iter = 0; iter !=  myDepth(); ++iter)//throws an error if not a pommaret basis or homogenous
  {
    res.push_back(x[iter]);
  }
  return res;
}

std::vector<RingElem> JBMill::myMaxStronglyIndependentSet() const
{
  std::vector<RingElem> res;
  if(IsOne(LPP(mySetT.begin()->myGetPol())))
  {
    return res;
  }
  ring quotientRing = NewQuotientRing(myOptions.myRingOptions.myPolyRing, ideal(myOptions.myRingOptions.myPolyRing, myReturnGB()) );
  const RingHom phi = QuotientingHom(quotientRing);
  if(IamPommaretBasis())
  {
    long CountIndets = NumIndets(myOptions.myRingOptions.myPolyRing);
    std::vector<RingElem> x = indets(myOptions.myRingOptions.myPolyRing);
    for(long iter = CountIndets - myDim(); iter != CountIndets; ++iter)
    {
      res.push_back(phi(x[iter]));
    }
  }
  else
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis or the ideal isn't homogenous", "myMaxStronglyIndependentSet");
  }
  return res;
}

bool JBMill::IamCohenMacaulay() const
{
  if(myDepth() == myDim()) //throws an error if it is not a pommaret basis or not homogenous
  {
    return true;
  }

  return false;
}

long JBMill::myRegularity() const
{
  if(!(IsStdDegRevLex(ordering(myOptions.myRingOptions.myPPMValue))))
  {
    CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
  }
  if(!IamHomogenous() || !IamPommaretBasis())
  {
    CoCoA_ERROR("The ideal isn't homogenous or pommaret", "myRegularity");
  }

  long MaxDeg(0);
  for(std::list<JanetTriple>::const_iterator iter = mySetT.begin(); iter != mySetT.end(); ++iter)
  {
    if(MaxDeg < StdDeg(iter->myGetPol()))
    {
      MaxDeg = StdDeg(iter->myGetPol());
    }
  }
  return MaxDeg;
}

long JBMill::myCastelnuovoMumfordRegularity() const
{
  return myRegularity(); //throws an error if isn't a pommaret basis
}

std::vector<RingElem> JBMill::mySaturation() const
{
  if(!(IamPommaretBasis() && IamHomogenous()))
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis or the ideal isn't homogenous", "mySaturation");
  }
  if(!(IsStdDegRevLex(ordering(myOptions.myRingOptions.myPPMValue))))
  {
    CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
  }
  // initialization
  std::multimap<long, RingElem> classes = myPommaretClasses();
  long CountIndets = NumIndets(myOptions.myRingOptions.myPolyRing);
  std::vector<RingElem> res;
  //handling H_1
  for (std::multimap<long, RingElem>::iterator iter = classes.lower_bound(1); iter != classes.upper_bound(1); ++iter)
  {
    res.push_back(iter->second / IndetPower(myOptions.myRingOptions.myPolyRing, CountIndets - 1, exponent(LPP(iter->second), CountIndets - 1)));
  }
  //handling rest of H_n
  for (int i = 2; i <= CountIndets; ++i)
  {
    for (std::multimap<long, RingElem>::iterator iter = classes.lower_bound(i); iter != classes.upper_bound(i); ++iter)
    {
      res.push_back(iter->second);
    }
  }
  return res;
}

std::multimap<long, RingElem> JBMill::myPommaretClasses() const
{
  if(!IamPommaretBasis())
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis", "myPommaretClasses");
  }
  std::map<PPMonoidElem, std::vector<bool> > MultVars = myComputeMultVars();
  std::multimap<long, RingElem> res;
  for(std::map<PPMonoidElem, std::vector<bool> >::iterator iter = MultVars.begin(); iter != MultVars.end(); ++iter)
  {
    long cls = myCountTrues(iter->second);
    RingElem RingElement(myOptions.myRingOptions.myPolyRing);
    for(std::list<JanetTriple>::const_iterator i = mySetT.begin(); i != mySetT.end(); ++i)
    {
      if(LPP(i->myGetPol()) == iter->first)
      {
        RingElement = i->myGetPol();
        break;
      }
    }
    res.insert(std::pair<long, RingElem>(cls, RingElement));
  }
  return res;
}

long JBMill::myDegPommaretClass(long index) const
{
  if(!IamPommaretBasis())
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis", "myDegPommaretClass");
  }
  std::multimap<long, RingElem> PomClasses(myPommaretClasses());
  long res(-1);
  for(std::multimap<long, RingElem>::iterator iter = PomClasses.lower_bound(index); iter != PomClasses.upper_bound(index); ++iter)
  {
    if(res < StdDeg(iter->second))
    {
      res = StdDeg(iter->second);
    }
  }
  return res;
}

std::map<long, long> JBMill::myDegPommaretClasses() const
{
  if(!IamPommaretBasis())
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis", "myDegPommaretClasses");
  }
  long CountIndets(NumIndets(myOptions.myRingOptions.myPolyRing));
  std::map<long, long> res;
  for(long iter = 0; iter != CountIndets; ++iter)
  {
    res.insert(std::pair<long, long>(iter + 1, myDegPommaretClass(iter + 1)));
  }
  return res;
}
  

long JBMill::mySatiety() const
{
  if(!(IamPommaretBasis() && IamHomogenous()))
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis or the ideal isn't homogenous", "mySatiety");
  }
  if(!(IsStdDegRevLex(ordering(myOptions.myRingOptions.myPPMValue))))
  {
    CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
  }
  return myDegPommaretClass(1);
}
 
void JBMill::myComplementaryDecompositionRecPart(std::vector< std::pair<PPMonoidElem, std::vector<bool> > >& output, JanetIterator JIter) const
{
  JanetIterator DegCounter(JIter);
  do
  {
    JanetIterator TmpIter(DegCounter);
    if( 0 != TmpIter.myNextVar())
    {
      myComplementaryDecompositionRecPart(output, TmpIter);
    }
  } while(DegCounter.myNextDeg() != 0);
  JanetIterator HighestNode = JIter.myGotoHighestNode(); //HighestNode == my
  std::vector<bool> n = myComputeMultVar(HighestNode); //mult or nonmult?
  n[JIter.myCurrentVar()] = false;
  for(long iter = JIter.myCurrentVar() + 1; iter != NumIndets(myOptions.myRingOptions.myPolyRing); ++iter)
  {
    n[iter] = true;
  }
  PPMonoidElem highestNodePPelem(HighestNode.myGetMonomial());
  if(JIter.myDisNextVar() == 0)
  {
    output.push_back(std::pair<PPMonoidElem, std::vector<bool> >(JIter.myGetMonomial(), n));
  }
  while(JIter.myCurrentDeg() < exponent(highestNodePPelem, JIter.myCurrentVar()))
  {
    long OldDeg = JIter.myCurrentDeg();
    ++OldDeg;
    PPMonoidElem pp = JIter.myGetMonomial();
    JIter.myNextDeg();
    long NewDeg = JIter.myCurrentDeg();
    for (int i = 1; i <= NewDeg - OldDeg; ++i)
    {
      output.push_back(std::pair<PPMonoidElem, std::vector<bool> >(pp * IndetPower(owner(pp), JIter.myCurrentVar(), i), n));
    }
  }
}

std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecomposition() const
{
  if(!IamMonomialIdeal())
  {
    CoCoA_ERROR("ideal isn't monomial", "myComplementaryDecomposition");
  }
  return myComplementaryDecompositionLeadingIdeal();
}

std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecompositionPolynomial() const
{
  return myComplementaryDecompositionLeadingIdeal();  
}

std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myComplementaryDecompositionLeadingIdeal() const
{
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > output;
  if(!IsOne(LPP(mySetT.begin()->myGetPol())))
  {
    JanetIterator iter(myJTree);
    myComplementaryDecompositionRecPart(output, iter);
  }
  return output;
}


std::vector< std::pair<PPMonoidElem, std::vector<bool> > > JBMill::myStandardPairs() const
{
  if(!IamMonomialIdeal())
  {
    CoCoA_ERROR("ideal isn't monomial", "myStandardPairs");
  }
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = myComplementaryDecomposition();
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > res;
  for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin(); iter != CompDecomp.end(); ++iter)
  {
    bool NBarEmpty = true;
    for(long i = 0; i != NumIndets(myOptions.myRingOptions.myPolyRing); ++i)
    {
      NBarEmpty = NBarEmpty && !(iter->second[i] && (0 != exponent(iter->first, i)));     
    }
    if(NBarEmpty)
    {
      res.push_back(*iter);
    }
    else
    {
      PPMonoid pp = PPM(myOptions.myRingOptions.myPolyRing);
      PPMonoidElem ResultPPMElem(one(pp));
      for(long i = 0; i != NumIndets(myOptions.myRingOptions.myPolyRing); ++i)
      {
        if(!(iter->second)[i])
        {
          ResultPPMElem *= IndetPower(pp, i, exponent(iter->first, i));
        }
      }
      res.push_back(std::pair<PPMonoidElem, std::vector<bool> >(ResultPPMElem, iter->second));
    }
  }
  return res;
}

std::pair<RingHom, std::vector<bool> > JBMill::myNoetherNormalization() const
{
  //construnction of decomposition
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = myComplementaryDecompositionPolynomial();
  //construction of nu
  std::vector<bool> nu;
  if (len(CompDecomp) > 1)
  {
    std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin();
    std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator TmpIter(iter);
    ++iter;
    nu = myUnionBoolVectors(TmpIter->second, iter->second);
    ++iter;
    for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = iter; i != CompDecomp.end(); ++i)
    {
      nu = myUnionBoolVectors(nu, i->second);  
    }
  } 
  else
  {
    std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator iter = CompDecomp.begin();
    nu = iter->second;     
  }
  //construction of d
  long d = len(nu);
  for (std::vector<bool>::reverse_iterator i = nu.rbegin(); i != nu.rend(); ++i)
  {
    if(*i == false)
    {
      break;
    }
    else
    {
      --d;
    }
  }
  if(myDim() == d)
  {
    return std::pair<RingHom, std::vector<bool> >(IdentityHom(myOptions.myRingOptions.myPolyRing), nu);
  }
  RingElem z = myGreatestMultVar(nu);
  std::vector<RingElem> x = indets(myOptions.myRingOptions.myPolyRing);
  //choose p
  RingElem p(myOptions.myRingOptions.myPolyRing);

  long AddZ(1);
  long ChoosePoly(0);
  while(true)
  {
    long TmpChoosePoly(ChoosePoly);
    for(std::list<JanetTriple>::const_iterator i = mySetT.begin(); i != mySetT.end(); ++i)
    {
      if(IamMultipleOfMultVars(LPP(i->myGetPol()), nu) == true)
      {
        if(TmpChoosePoly == 0)
        { 
          p = i->myGetPol();
          --TmpChoosePoly;
          break;        
        }
        else
        {
          --TmpChoosePoly;
        }
      }
    }
    if((TmpChoosePoly != -1) && (ChoosePoly != 0))
    {
      CoCoA_ERROR("Not able to compute noether normalization", "myNoetherNormalization");
    }
    //construction of homomorphism
    std::vector<RingElem> xImages;
    long n = NumIndets(myOptions.myRingOptions.myPolyRing);
    if(IsZero(p) ==  true)
    {
      for (long i = 0; i < n; ++i)
      {
        xImages.push_back(zero(myOptions.myRingOptions.myPolyRing));
      } 
      RingHom varphi = PolyAlgebraHom(myOptions.myRingOptions.myPolyRing, myOptions.myRingOptions.myPolyRing, xImages);
      return std::pair<RingHom, std::vector<bool> >(varphi, nu);
    }
    for (long i = 0; i < n; ++i)
    {
      if((indet(myOptions.myRingOptions.myPolyRing, i) == z) || exponent(LPP(p), i) == 0)
      {
        xImages.push_back(indet(myOptions.myRingOptions.myPolyRing, i)); 
      }
      else
      {
        RingElem poly(indet(myOptions.myRingOptions.myPolyRing, i) - z);
        for (long i = 1; i < AddZ; ++i)
        {
          poly = poly - z;
        }
        if(IsZero(poly) == true)
        {
          xImages.clear();
          for (long i = 0; i < n; ++i)
          {
            xImages.push_back(zero(myOptions.myRingOptions.myPolyRing));
          } 
          RingHom varphi = PolyAlgebraHom(myOptions.myRingOptions.myPolyRing, myOptions.myRingOptions.myPolyRing, xImages);
        }
        xImages.push_back(poly);
      }
    }
    RingHom varphi = PolyAlgebraHom(myOptions.myRingOptions.myPolyRing, myOptions.myRingOptions.myPolyRing, xImages);
    //recursion step
    std::vector<RingElem> VarphiImage;
    for(std::list<JanetTriple>::const_iterator i = mySetT.begin(); i != mySetT.end(); ++i)
    {
      VarphiImage.push_back(varphi(i->myGetPol()));
    }
    JBMill TmpMill = ExtendedJanetBasis(VarphiImage);
    std::pair<RingHom, std::vector<bool> > noether = TmpMill.myNoetherNormalization();
    //output
    if(IsZero(noether.first(indet(myOptions.myRingOptions.myPolyRing, 0))) == false)
    {
      return std::pair<RingHom, std::vector<bool> >(noether.first(varphi), noether.second);    
    }
    else
    {
      if(AddZ == 10)
      {
        AddZ = 0;
        ++ChoosePoly;
      }
      else
      {
        ++AddZ;
      }
    }
  }
}

bool JBMill::IamMultipleOfMultVars(PPMonoidElem ppm, std::vector<bool> v) const
{
  long n = NumIndets(myOptions.myRingOptions.myPolyRing);
  for (long i = 0; i < n; ++i)
  {
    if(exponent(ppm, i) > 0)
    {
      if(v[i] == false)
      {
        return false;
      }
    }
  }
  return true;
}


RingElem JBMill::myGreatestMultVar(std::vector<bool> v) const
{
  long n = NumIndets(myOptions.myRingOptions.myPolyRing);
  long MaxVar(-1);
  for(long i = 0; i < n; ++i)
  {
    if(v[i] == true)
    {
      if(MaxVar != -1)
      {
        if(LPP(indet(myOptions.myRingOptions.myPolyRing, i)) > LPP(indet(myOptions.myRingOptions.myPolyRing, MaxVar)))
        {
          MaxVar = i;
        }
      }
      else
      {
        MaxVar = i;
      }
    }
  }
  return indet(myOptions.myRingOptions.myPolyRing, MaxVar);
}

std::vector<bool> JBMill::myUnionBoolVectors(std::vector<bool> v1, std::vector<bool> v2) const
{
  std::vector<bool> result;
  std::vector<bool>::iterator iter2 = v2.begin();
  for (std::vector<bool>::iterator iter1 = v1.begin(); iter1 != v1.end(); ++iter1)
  {
    if(*iter1 == true)
    {
      result.push_back(true);
    }
    else
    {
      if(*iter2 == true)
      {
        result.push_back(true);
      }
      else
      {
        result.push_back(false);
      }
    }
    ++iter2;
  }
  return result;
}


std::vector<bool> JBMill::myComputeMultVar(JanetIterator iter) const
{
  std::map<PPMonoidElem, std::vector<bool> > MultVars = myComputeMultVars();
  std::map<PPMonoidElem, std::vector<bool> >::iterator MapIter = MultVars.find(iter.myGetMonomial());
  return MapIter->second;
}

std::vector<bool> JBMill::myComputeNonMultVar(JanetIterator iter) const
{
  return myReverseBoolVec(myComputeMultVar(iter));
}

std::vector<bool> JBMill::myReverseBoolVec(std::vector<bool> vec) const
{
  std::vector<bool> res;
  for(std::vector<bool>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
  {
    res.push_back(!(*iter));
  }
  return res;
}

long JBMill::myCls(PPMonoidElem elem) const
{
  long CountIndets = NumIndets(owner(elem));
  for (long i = CountIndets; i != 0; --i)
  {
    if(exponent(elem, i - 1) != 0)
    {
      return NumIndets(owner(elem)) - i + 1;
    }
  }
  return CountIndets;
}

long JBMill::myOldCls(PPMonoidElem elem) const
{
  return NumIndets(owner(elem)) - myCls(elem) - 1;
}


long JBMill::myMinCls() const
{
  long indets = NumIndets(myOptions.myRingOptions.myPolyRing);
  long min(indets);
  for(std::list<JanetTriple>::const_iterator i = mySetT.begin(); i != mySetT.end(); ++i)
  {
    if(min > myCls(LPP(i->myGetPol())))
    {
      min = myCls(LPP(i->myGetPol()));
    }
  }
  return min;
}

std::vector<RingElem> JBMill::myElementsWithClass(long InputCls) const
{
  std::vector<RingElem> res;
  for(std::list<JanetTriple>::const_iterator i = mySetT.begin(); i != mySetT.end(); ++i)
  {
    if(InputCls == myCls(LPP(i->myGetPol())))
    {
      res.push_back(i->myGetPol());
    }
  }
  return res;
}


std::vector<RingElem> JBMill::mySocle() const
{
  if(!(IamPommaretBasis()))
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis", "mySocle");
  }
  if(!(IamHomogenous()))
  {
    CoCoA_ERROR("The ideal isn't homogenous", "mySocle");
  }

  if(!(IamCohenMacaulay()))
  {
    CoCoA_ERROR("Ideal isn't CohenMacaulay", "mySocle");    
  }
  if(!(IsStdDegRevLex(ordering(myOptions.myRingOptions.myPPMValue))))
  {
    CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
  }

  long min = myMinCls();
  std::vector<RingElem> socle = myElementsWithClass(min);
  std::vector<RingElem> ResultWithoutResidueClass;
  for (std::vector<RingElem>::iterator i = socle.begin(); i != socle.end(); ++i)
  {
    // RingElem rest = myStandardRepresentation(*i).second;
    bool NoConstant(true);
    for (int iter = 0; iter < (min - 1); ++iter)
    {
      RingElem rest = myStandardRepresentation((*i) * monomial(myOptions.myRingOptions.myPolyRing,one(myOptions.myRingOptions.myCoefRing  ), indet(myOptions.myRingOptions.myPPMValue, iter))).second;
      while(!IsZero(rest))
      {
        if(IsConstant(rest))
        {
          NoConstant = false;
          break;
        }
        rest = rest - monomial(myOptions.myRingOptions.myPolyRing,LC(rest),LPP(rest));
      }
      if(!NoConstant)
      {
        break;
      }
    }
    if(NoConstant)
    {
      ResultWithoutResidueClass.push_back((*i) / monomial(myOptions.myRingOptions.myPolyRing,one(myOptions.myRingOptions.myCoefRing  ), indet(myOptions.myRingOptions.myPPMValue, min - 1)));
    }
  }
  long dimension = myDim();
  std::vector<RingElem> x = indets(myOptions.myRingOptions.myPolyRing);
  std::vector<RingElem> GenSet;
  for (long i = 0; i < dimension; ++i)
  {
    GenSet.push_back(x[NumIndets(myOptions.myRingOptions.myPolyRing) - 1 - i]);
  }
  ring quotientRing = NewQuotientRing(myOptions.myRingOptions.myPolyRing, ideal(myOptions.myRingOptions.myPolyRing, GenSet) );
  const RingHom phi = QuotientingHom(quotientRing);
  std::vector<RingElem> res;
  for (std::vector<RingElem>::iterator i = ResultWithoutResidueClass.begin(); i != ResultWithoutResidueClass.end(); ++i)
  {
    res.push_back(phi(*i));
  }
  return res;
}

std::map<std::pair<long, long>, long> JBMill::myExtremalBettiNumbers() const
{
  if(!(IamPommaretBasis()))
  {
    CoCoA_ERROR("Janet basis isn't a Pommaret basis", "myExtremalBettiNumbers");
  }
  if(!(IamHomogenous()))
  {
    CoCoA_ERROR("Janet basis isn't homogenous", "myExtremalBettiNumbers");
  }
  if(!(IsStdDegRevLex(ordering(myOptions.myRingOptions.myPPMValue))))
  {
    CoCoA_ERROR(ERR::PPOrder, "It must be the degrevlex ordering!!!");
  }

  long CurPDim = myProjDim();
  long n(NumIndets(myOptions.myRingOptions.myPolyRing));
  std::map<long, long> PommCls(myDegPommaretClasses());
  long deg(-1);
  long cls(-1);
  for(long i = 1; i != n + 1; ++i)
  {
    if(deg < PommCls[i])
    {
      deg = PommCls[i];
      cls = i;
    }
  }
  std::map<std::pair<long, long>, long> result;
  std::vector<RingElem>  ElementsWithClass(myElementsWithClass(cls));
  result.insert(std::pair<std::pair<long, long>, long>(std::pair<long,long>(n - cls, deg + n - cls), myCountElementsWithDegree(ElementsWithClass, deg)));
  while( n - cls < CurPDim)
  {
    deg = 0;
    long OldCls(cls);
    for(long i = 1; i != OldCls; ++i)
    {
      if(deg < PommCls[i])
      {
        deg = PommCls[i];
        cls = i;
      }
    }
    ElementsWithClass = myElementsWithClass(cls);
    result.insert(std::pair<std::pair<long, long>, long>(std::pair<long,long>(n - cls, deg + n - cls), myCountElementsWithDegree(ElementsWithClass, deg)));
  }
  return result;
}

long JBMill::myCountElementsWithDegree(const std::vector<RingElem>& vec, long degree) const
{
  long counter(0);
  for (std::vector<RingElem>::const_iterator i = vec.begin(); i != vec.end(); ++i)
  {
    if(StdDeg(*i) == degree)
    {
      ++counter;
    }
  }
  return counter;
}
