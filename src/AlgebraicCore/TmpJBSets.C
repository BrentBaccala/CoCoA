#include "CoCoA/TmpJBSets.H"
#include "CoCoA/utils.H"

using namespace CoCoA;

JBSets::MultisetIterator JBSets::myBeginSetT() const
{
  return mySetT.begin();
}

JBSets::MultisetIterator JBSets::myEndSetT() const
{
  return mySetT.end();
}

JBSets::MultisetIterator JBSets::myBeginSetQ() const
{
  return mySetQ.begin();
}

JBSets::MultisetIterator JBSets::myEndSetQ() const
{
  return mySetQ.end();
}

JBSets::MultisetIterator JBSets::myBeginSetP() const
{
  return mySetP.begin();
}

JBSets::MultisetIterator JBSets::myEndSetP() const
{
  return mySetP.end();
}


void JBSets::myInsertSetQ(JanetTriple& triple)
{
  //at the moment we don't look if we have already this element -> maybe we have to change this!!!
  myBasicSet.push_front(triple);
  mySetQ.insert(myBasicSet.begin());
}


//TODO naming of iterators
Iterator JBSets::myMinimizeAndInsertSetT(JanetTree& JTree)
{
  //choosing the smalles element from mySetQ
  std::pair<JBSets::MultisetIterator,JBSets::MultisetIterator> IterPair(mySetQ.equal_range(*(mySetQ.begin())));
  JBSets::MultisetIterator iter(mySetQ.begin());
  for(JBSets::MultisetIterator ForIter = IterPair.first; ForIter != IterPair.second; ++ForIter)
  {
    if(NumTerms((*ForIter)->myGetPol()) <= NumTerms((*iter)->myGetPol()))
    {
      if(NumTerms((*ForIter)->myGetPol()) < NumTerms((*iter)->myGetPol()))
      {
        iter = ForIter;
      }
      else
      {
        if(IsOrderedDomain(CoeffRing(myPolyRing)))
        {
          if(abs(LC((*ForIter)->myGetPol())) < abs(LC((*iter)->myGetPol())))
          {
            iter = ForIter;
          }
        }
      }
    }
  }
  //deletes the choosen element from mySetQ
  Iterator BasicSetIter(*iter);
  mySetQ.erase(iter);

  //perform a tail reduction at the choosen element
  myJTailNormalForm(JTree,BasicSetIter);
  PPMonoidElem ppelem(LPP(BasicSetIter->myGetPol()));

  //if choosen element isn't a prolongation check if there is a bigger element in mySetT which is divided by the choosen elemnt
  //if there is such an element every element which is bigger then the choosen element is moved to mySetQ 
  if(ppelem == BasicSetIter->myGetAnc())
  {
    MultisetIterator iter(mySetT.upper_bound(BasicSetIter));
    bool DeleteElems(false);
    while(iter != mySetT.end())
    {
      if(IsDivisible(myPolyRing->myLPP(raw(*((*iter)->myGetPolPtr()))),ppelem))
      {
        DeleteElems = true;
        break;
      }
      ++iter;
    }
    if(DeleteElems)
    {
      iter = mySetT.upper_bound(BasicSetIter);
      if(iter != mySetT.end())
      {
        mySetQ.insert(iter,mySetT.end());
        mySetT.erase(iter,mySetT.end());
      }
    }
  }

  //insert the choosen element to mySetT
  mySetT.insert(BasicSetIter);
  return BasicSetIter;
}


long JBSets::myHeadReduceSetQ(JanetTree& JTree)
{
  //copy mySetQ into a temporary multiset
  std::multiset<Iterator,CompareIterator> TempSet(mySetQ);
  //clearing mySetQ
  mySetQ.clear();
  for(MultisetIterator iter = TempSet.begin(); iter != TempSet.end(); ++iter)
  {
    bool ExitViaBreak(false);
    PPMonoidElem IterLPP(LPP(*((*iter)->myGetPolPtr()))); //we using IterLPP only two times...
    JanetTriple* GStart = JTree.myJDivisor(IterLPP);
    if(GStart != 0)
    {
      //testing if we can apply an involutive criteria
      if((IterLPP != (*iter)->myGetAnc()) && (IamCriteria(*(*iter),*GStart)))
      {
        // if true element reduces to zero -> remove element from myBasicSet
        myBasicSet.erase(*iter);
        continue;
      }
      else
      {
        //involutive headreduction
        myReductionCog->myAssignReset(*((*iter)->myGetPolPtr()));
        myReductionCog->myReduce(*(GStart->myGetPolPtr()));
        while(!IsActiveZero(myReductionCog))
        {
          JanetTriple* g = JTree.myJDivisor(ActiveLPP(myReductionCog));
          if(g != 0)
          {
            myReductionCog->myReduce(*(g->myGetPolPtr()));
          }
          else
          {
            ExitViaBreak = true;
            break;
          }
        }
        if(ExitViaBreak == true)
        {
          //coefficient stuff
          myReductionCog->myRelease(*((*iter)->myGetPolPtr())); 
          if(IsFractionField(CoeffRing(myPolyRing)))
          {
            myPolyRing->myDivByCoeff(raw(*((*iter)->myGetPolPtr())),raw(LC(*((*iter)->myGetPolPtr()))));
            (*iter)->mySetPol(ClearDenom(*((*iter)->myGetPolPtr())));
          }
          //adepting ancestor and nonmultiplicative variables
          (*iter)->mySetAnc(myPolyRing->myLPP(raw(*(*iter)->myGetPolPtr())));
          (*iter)->myClearNM();
          mySetQ.insert(*iter);
        }
        else
        {
          // element reduces to zero -> remove element from myBasicSet
          myBasicSet.erase(*iter);
        }
      }
    }
    else
    {
      // inserting head reduced element to mySetQ
      mySetQ.insert(*iter);
    }
  }
  return len(mySetQ);
}


void JBSets::myJTailNormalForm(JanetTree& JTree, Iterator iter)
{
  //skipping the leading-monomial (we want to perform a !tail!-reduction)
  myReductionCog->myAssignReset(*(iter->myGetPolPtr()));
  myReductionCog->myMoveToNextLM();
  while(!IsActiveZero(myReductionCog))
  {
    JanetTriple* g = JTree.myJDivisor(ActiveLPP(myReductionCog));
    if(g != 0)
    {
      myReductionCog->myReduce(*(g->myGetPolPtr()));
    }
    else
    {
      myReductionCog->myMoveToNextLM();
    }
  }
  myReductionCog->myRelease(*(iter->myGetPolPtr()));
  if(IsFractionField(CoeffRing(myPolyRing)))
  {
    myPolyRing->myDivByCoeff(raw(*(iter->myGetPolPtr())),raw(LC(*(iter->myGetPolPtr()))));
    iter->mySetPol(ClearDenom(*(iter->myGetPolPtr())));
  }
}


void JBSets::myTailNormalForm(Iterator ToRed, Iterator red)
{
  PPMonoidElem lpp(LPP(red->myGetPol()));
  myReductionCog->myAssignReset(*(ToRed->myGetPolPtr()));
  myReductionCog->myMoveToNextLM();
  while((!IsActiveZero(myReductionCog))&&(lpp <= ActiveLPP(myReductionCog)))
  {
    if(ActiveLPP(myReductionCog) == lpp)
    {
      myReductionCog->myReduce(*(red->myGetPolPtr()));
    }
    else
    {
      myReductionCog->myMoveToNextLM();
    }
  }
  myReductionCog->myRelease(*(ToRed->myGetPolPtr()));
}



void JBSets::myTailReduceFirstElemSetQ(JanetTree& JTree)
{
  myJTailNormalForm(JTree,*(mySetQ.begin()));
}

void JBSets::myTailReduceSetT(JanetTree& JTree, const std::list<JanetTriple>::iterator& TripleIter)
{
  //choosing the starting point in mySetT
  MultisetIterator beginner(mySetT.upper_bound(TripleIter));
  for(MultisetIterator iter = beginner; iter != mySetT.end(); ++iter)
  {
    myJTailNormalForm(JTree, *iter);
  }
}

void JBSets::myTailReduceSetTAll(JanetTree& JTree)
{
  for(MultisetIterator iter = mySetT.begin(); iter != mySetT.end(); ++iter)
  {
    myJTailNormalForm(JTree, *iter);
  }
}

long JBSets::mySizeSetT()
{
  return len(mySetT);
}

bool JBSets::IamCriteria(const JanetTriple& p, const JanetTriple& g) const 
{
  if(myCriteria.test(0))
  {
    //first involutive criteria
    if( (p.myGetAnc() * g.myGetAnc()) == myPolyRing->myLPP(raw(*p.myGetPolPtr())))
    {
      return true;
    }
  }
  if(myCriteria.test(1))
  {
    //second involutive criteria
    PPMonoidElem ppm(lcm(p.myGetAnc(), g.myGetAnc()));
    if(IsDivisible(ppm, myPolyRing->myLPP(raw(*p.myGetPolPtr()))))
    {
      if(ppm != myPolyRing->myLPP(raw(*p.myGetPolPtr())))
      {
        return true;
      }
    }
  }
  if(myCriteria.test(2))
  {
    //third involutive criteria
    if(IamCriteria3(p.myGetAnc(),g.myGetAnc()))
    {
      return true;
    }
  }
  return false;
}

bool JBSets::IamCriteria3(const PPMonoidElem& AncP, const PPMonoidElem& AncG) const
{  
  MultisetIterator iter(mySetT.begin());
  PPMonoidElem lcmpg(lcm(AncP, AncG));
  while((iter != mySetT.end()) && (myPolyRing->myLPP(raw(*((*iter)->myGetPolPtr()))) < lcmpg))
  {
    PPMonoidElem h(myPolyRing->myLPP(raw(*((*iter)->myGetPolPtr()))));
    PPMonoidElem lcmhp(lcm(h, AncP));
    PPMonoidElem lcmhg(lcm(h, AncG));
    if(IsDivisible(lcmpg, lcmhp) && IsDivisible(lcmpg, lcmhg))
    {
      if((lcmhp != lcmpg) && (lcmhg != lcmpg))
      {
        return true;
      }
    }
    ++iter;
  }
  return false;
}

void JBSets::myMoveFromQtoP(JanetTree& JTree)
{
  MultisetIterator iter(mySetQ.begin());
  long deg(StdDeg(LPP((*iter)->myGetPol())));
  //iterating over all elements in mySetQ with degree deg
  do //  until (iter != mySetQ.end()) && (StdDeg(LPP((*iter)->myGetPol())) == deg is false
  {
    bool ExitViaBreak(false);
    PPMonoidElem IterLPP(LPP(*((*iter)->myGetPolPtr()))); //we using IterLPP only two times...
    // performing a full j-reduction
    JanetTriple* GStart = JTree.myJDivisor(IterLPP);
    if(GStart != 0)
    {
      if((IterLPP != (*iter)->myGetAnc()) && (IamCriteria(*(*iter),*GStart)))
      {
        myBasicSet.erase(*iter);
        ++iter;
      }
      else
      {
        myReductionCog->myAssignReset(*((*iter)->myGetPolPtr()));
        myReductionCog->myReduce(*(GStart->myGetPolPtr()));
        while(!IsActiveZero(myReductionCog))
        {
          JanetTriple* g = JTree.myJDivisor(ActiveLPP(myReductionCog));
          if(g != 0)
          {
            myReductionCog->myReduce(*(g->myGetPolPtr()));
          }
          else
          {
            ExitViaBreak = true;
            break;
          }
        }
        if(ExitViaBreak == true)
        {
          myReductionCog->myRelease(*((*iter)->myGetPolPtr()));  
          (*iter)->mySetAnc(myPolyRing->myLPP(raw(*(*iter)->myGetPolPtr())));
          (*iter)->myClearNM();
          myJTailNormalForm(JTree, *iter);
          mySetP.insert(*iter);
          ++iter;
        }
        else
        {
          myBasicSet.erase(*iter);
          ++iter;
        }
      }
    }
    else
    {
      myJTailNormalForm(JTree, *iter);
      mySetP.insert(*iter);
      ++iter;
    }

  }
  while((iter != mySetQ.end()) && (StdDeg(LPP((*iter)->myGetPol())) == deg));
  mySetQ.erase(mySetQ.begin(),iter);
}


void JBSets::myJFullNormalForm(JanetTree& JTree, Iterator iter)
{
  myReductionCog->myAssignReset(*(iter->myGetPolPtr()));
  myReductionCog->myMoveToNextLM();
  while(!IsActiveZero(myReductionCog))
  {
    JanetTriple* g = JTree.myJDivisor(ActiveLPP(myReductionCog));
    if(g != 0)
    {
      myReductionCog->myReduce(*(g->myGetPolPtr()));
    }
    else
    {
      myReductionCog->myMoveToNextLM();
    }
  }
  myReductionCog->myRelease(*(iter->myGetPolPtr()));
  if(IsFractionField(CoeffRing(myPolyRing)))
  {
    myPolyRing->myDivByCoeff(raw(*(iter->myGetPolPtr())),raw(LC(*(iter->myGetPolPtr()))));
    iter->mySetPol(ClearDenom(*(iter->myGetPolPtr())));
  }
}

void JBSets::myJUpdateHighPart(std::list<Iterator>& stack)
{
  do
  {
    //choose element from mySetP
    JBSets::MultisetIterator ChooseIter;
    ChooseIter = mySetP.end();
    --ChooseIter;
    myChoosingJUpdate(ChooseIter);
    stack.push_front(*ChooseIter);
    mySetP.erase(ChooseIter);
    MultisetIterator PIter(mySetP.begin());
    //reduce every element in mySetP w.r.t. stack.begin()
    while(PIter != mySetP.end())
    {
      if(LPP((*PIter)->myGetPol()) != LPP((*(stack.begin()))->myGetPol()))
      {
        ++PIter;
      }
      else
      {
        myReductionStepJUpdate(PIter, stack.begin());
      }
    }
    myRecomputeOrderingSetP();
  } while(!mySetP.empty());
}

void JBSets::myJUpdateLowPart(std::list<Iterator>& stack)
{
  do
  {
    //choosing element from mySetP
    JBSets::MultisetIterator ChooseIter;
    ChooseIter = mySetP.begin();
    myChoosingJUpdate(ChooseIter);
    stack.push_front(*ChooseIter);
    mySetP.erase(ChooseIter);
    std::list<Iterator>::iterator StackIter(stack.begin());
    ++StackIter;
    //everything which is bigger then the begin of stack move to mySetP
    while((StackIter != stack.end()) && LPP((*StackIter)->myGetPol()) > LPP((*(stack.begin()))->myGetPol()))
    {
      std::list<Iterator>::iterator TempStackIter(StackIter);
      ++StackIter;
      mySetP.insert(*TempStackIter);
      stack.erase(TempStackIter);
    }
    //reduce every element in mySetP w.r.t. stack.begin()
    MultisetIterator PIter(mySetP.begin());
    while(PIter != mySetP.end())
    {
      if(LPP((*PIter)->myGetPol()) > LPP((*(stack.begin()))->myGetPol()))
      {
        PPMonoidElem ppm(LPP((*PIter)->myGetPol()));
        myTailNormalForm(*PIter, *(stack.begin()));
        if(IsZero((*PIter)->myGetPol()))
        {
          myBasicSet.erase(*PIter);
//???useless          MultisetIterator TempIter(PIter);
          ++PIter;
        }
        else
        {
          if(ppm != LPP((*PIter)->myGetPol()))
          { 
            (*PIter)->mySetAnc(LPP((*PIter)->myGetPol()));
            (*PIter)->myClearNM();
          }
          ++PIter;
        }
      }
      else
      {
        myReductionStepJUpdate(PIter, stack.begin());
      }
    }
    myRecomputeOrderingSetP();
  } while(!mySetP.empty());
}

void JBSets::myChoosingJUpdate(JBSets::MultisetIterator& ChooseIter)
{
  //choosing strategy:
  // Top priority: number of terms
  // second priority: abs(leading coefficient) (if possible)
  std::pair<JBSets::MultisetIterator,JBSets::MultisetIterator> IterPair(mySetP.equal_range(*(ChooseIter)));
  for(JBSets::MultisetIterator ForIter = IterPair.first; ForIter != IterPair.second; ++ForIter)
  {
    if(NumTerms((*ForIter)->myGetPol()) <= NumTerms((*ChooseIter)->myGetPol()))
    {
      if(NumTerms((*ForIter)->myGetPol()) < NumTerms((*ChooseIter)->myGetPol()))
      {
        ChooseIter = ForIter;
      }
      else
      {
        if(IsOrderedDomain(CoeffRing(myPolyRing)))
        {
          if(abs(LC((*ForIter)->myGetPol())) < abs(LC((*ChooseIter)->myGetPol())))
          {
            ChooseIter = ForIter;
          }
        }
      }
    }
  }
}

void JBSets::myReductionStepJUpdate(MultisetIterator& PIter, std::list<Iterator>::iterator ReductionIter)
{
  //test if we can apply an involutive criteria
  if(IamCriteria(**PIter,**(ReductionIter)))
  {
    myBasicSet.erase(*PIter);
    MultisetIterator TempIter(PIter);
    ++PIter;
    mySetP.erase(TempIter);
  }
  else
  {
    //head reduction (only one time this is possible!)
    //ATTENTION THIS STEP DESTROY THE ORDER IN SetP
    myPolyRing->myReductionStep(raw(*((*PIter)->myGetPolPtr())),raw(*((*(ReductionIter))->myGetPolPtr())));
    if(!IsZero((*PIter)->myGetPol()))
    {
      //adepting ancestor etc.
      if(IsFractionField(CoeffRing(myPolyRing)))
      {
        myPolyRing->myDivByCoeff(raw(*((*PIter)->myGetPolPtr())),raw(LC(*((*PIter)->myGetPolPtr()))));
        myPolyRing->myClearDenom(raw(*((*PIter)->myGetPolPtr())),raw(*((*PIter)->myGetPolPtr())));
      }
      (*PIter)->mySetAnc(LPP(*((*PIter)->myGetPolPtr())));
      (*PIter)->myClearNM();
      ++PIter;
    }
    else
    {
      //reduced to zero -> erase element completly
      myBasicSet.erase(*PIter);
      MultisetIterator TempIter(PIter);
      ++PIter;
      mySetP.erase(TempIter);
    }
  }
}


void JBSets::myRecomputeOrderingSetP()
{
  //recompute the ordering of mySetP
  std::list<Iterator> TempList;
  TempList.insert(TempList.begin(),mySetP.begin(),mySetP.end());
  mySetP.clear();
  mySetP.insert(TempList.begin(),TempList.end());
}

bool JBSets::IamJUpdate()
{
  if(mySetP.empty())
  {
    return false;
  }
  //update process
  std::list<Iterator> stack;
  if(myHigh)
  {
    myJUpdateHighPart(stack);
  }
  else
  {
    myJUpdateLowPart(stack);
  }
  mySetP.insert(stack.begin(), stack.end());
  long deg(StdDeg((*(mySetP.begin()))->myGetPol()));
  MultisetIterator BeginDelIter(mySetP.begin());
  while(BeginDelIter != mySetP.end())
  {
    if(StdDeg((*BeginDelIter)->myGetPol()) > deg)
    {
      mySetQ.insert(BeginDelIter,mySetP.end());
      mySetP.erase(BeginDelIter,mySetP.end());
      break;
    }
    else
    {
      ++BeginDelIter;
    }
  }
  //check if we have to modify mySetT   
  for(MultisetIterator iter = mySetP.begin(); iter != mySetP.end(); ++iter)
  {
    PPMonoidElem lpp(LPP((*iter)->myGetPol()));
    if(lpp == (*iter)->myGetAnc())
    {
      for(std::multiset<Iterator,CompareIterator>::reverse_iterator SetTIter = mySetT.rbegin(); SetTIter != mySetT.rend(); ++SetTIter)
      {
        if(IsDivisible(LPP((*SetTIter)->myGetPol()),lpp))
        {
          return true;
        }
        else
        {
          if(lpp > LPP((*SetTIter)->myGetPol()))
          {
            break;
          }
        }
      }
    }
  }
  return false;
}



JanetTree JBSets::myOneLeafTree(const std::list<JanetTriple>::iterator& ps, long i, PPMonoidElem w)
{
  // if we say ps we mean myPolyRing->myLPP(raw(ps->myGetPol()))
  long n = NumIndets(myPolyRing);
  //navigates to the first index where deg_i(ps) != 0
  while((i < n - 1) && (exponent(myPolyRing->myLPP(raw(ps->myGetPol())),i) == 0)) //we don't reach n because the algorithms which use myOneLeafTree check this
  {
    ++i;
  } 
  //janet tree starts at deg_i(w) at variable_index i
  JanetTree LeafTree(myPolyRing, exponent(w, i) , i);
  JanetIterator iter(LeafTree);
  JanetIterator TempIter(iter);
  PPMonoidElem u(myPolyRing->myLPP(raw(ps->myGetPol())));
  while((u != w) && (i < n)) 
  {
    // add node in deg-direction
    if(exponent(w,i) < exponent(u,i)) 
    {
      iter.mySetNextDeg(exponent(u,i) - exponent(w,i));
      w *= IndetPower(PPM(myPolyRing),i,exponent(u,i) - exponent(w,i));
      iter.myNextDeg();
      TempIter = iter;
    }
    // add node in var-direction
    if(u != w && (iter.myCurrentVar() < (n - 1)))
    {
      if(i == n - 1)   
      {
        //should raise an error???????
        return JanetTree(myPolyRing, exponent(w, i), i);
      }
      // let i = iter.myCurrentVar()
      // deg_i(ps) = 0 -> skipping this variable in the janet tree
      if((TempIter.myCurrentVar()) != (iter.myCurrentVar())) 
      {
        TempIter.mySetNextVar(iter.myCurrentVar() - TempIter.myCurrentVar() + 1);
        iter = TempIter;
      }
      else
      {
        iter.mySetNextVar(1);
      }
      ++i;
      iter.myNextVar();
    }
  }
  // the node was a internal node, but it must be a leaf node
  iter.myChangeToLeafNode(ps);
  return LeafTree;
}

void JBSets::myJInsert(std::list<JanetTriple>::iterator& ps, JanetTree& JTree) //call-by-reference
{
  //initialization
  long i(0);
  JanetIterator iter(JTree);
  long n(NumIndets(myPolyRing));
  // main loop
  while(i < n)
  {
    //moving in deg direction until deg(iter) >= deg(ps) or there is no next node in deg direction
    while((iter.myCurrentDeg() < exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i)) && (iter.myDisNextDeg()))
    {
      iter.myNextDeg();
    }
    //if deg(iter) < deg(ps) (that means that there isn't a next node in deg direction)
    //'insert' the new element here and exit loop (i = n)
    if(iter.myCurrentDeg() < exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i))
    {
      myJProlong(iter, i);
      for(long j = i; j != n; ++j) 
      {
        ps->myDelNM(j);
      }
      JanetTree tree(myOneLeafTree(ps, i, iter.myGetMonomial() * indet(PPM(myPolyRing), i)));
      iter.myConnectJanetTreeDeg(tree);
      i = n ;
    }
    else
    { 
      // deg(iter) >= deg(ps)
      // first we adept nonmult. vars
      if ((iter.myDisNextDeg()) || (exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i) < iter.myCurrentDeg()))
      {
        if(!ps->IamNM(i))
        {  
          JanetTriple triple(ps->myGetPol() * monomial(myPolyRing, 1, indet(PPM(myPolyRing), i)), ps->myGetAnc());
          myInsertSetQ(triple);
          ps->mySetNM(i);
        }
      }
      else
      {
        ps->myDelNM(i);
      }

      // deg(iter) == deg(ps) and there is a node in variable direction 
      // (second condition is unnecessary but I feel better with it ;-) )
      if((exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i) == iter.myCurrentDeg())  && (iter.myDisNextVar()))
      {
        bool NewVarNode(false);
        long DisNextVarNode(0);
        //check if we need a new var node between the cur var and the next var
        for(long j = i + 1; j < (iter.myDisNextVar() + iter.myCurrentVar()); ++j)
        {
          NewVarNode = (exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), j) != 0);
          if(NewVarNode)
          {
            DisNextVarNode = j - i;
            break;
          }
        }
        if(NewVarNode)
        {
          iter.mySetNextVar(DisNextVarNode);
        }
        // we have to set all variables between the current and the old node to multiplicative... 
        for(long j = (i + 1); j < (i + iter.myDisNextVar()); ++j) //THE SOLUTION OF THE MYSTHERIOUS BUG! (after 3 month of searching!)
        {
          ps->myDelNM(j);
        }
        //goto next var
        i += iter.myNextVar();
      }
      else //exponent(LPP(ps.myGetPol()), i) > iter.myCurrentDeg()
      {
        // we have to create an new node with degree equal to deg(ps)
        // this node is an internal node!
        // because we create a new node we can construct the whole way to ps
        for(long j = i + 1; j != n; ++j)
        {
          ps->myDelNM(j);
        }
        if(exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i) < iter.myCurrentDeg())
        {
        iter.myPrevDeg();
          iter.mySetNextDeg(exponent(myPolyRing->myLPP(raw(*(ps->myGetPolPtr()))), i) - iter.myCurrentDeg());
          iter.myNextDeg();
        }

        JanetTree tree(myOneLeafTree(ps, i + 1, iter.myGetMonomial()));
        iter.mySetNextVar(tree.myGetBeginVar() - iter.myCurrentVar());
        iter.myNextVar();
        iter.myConnectJanetTreeDeg(tree);
        //leaving loop
        i = n;
      }
    }
  }
}

void JBSets::myJProlong(JanetIterator& iter, const long& index)
{
  //navigating in to each node
  while(iter.myDisNextDeg())
  {
    if(iter.myDisNextVar())
    {
      JanetIterator TempIter(iter);
      TempIter.myNextVar();
      myJProlong(TempIter, index);
    }
    iter.myNextDeg();
  }
  if(iter.myDisNextVar())
  {
    JanetIterator TempIter(iter);
    TempIter.myNextVar();
    myJProlong(TempIter, index);
  }
  else
  {
    //if we reach a leaf node prolong
    if(iter.myGetTriplePtr() != 0)
    {
      if(!iter.IamNM(index)) 
      {
        JanetTriple triple(iter.myGetPol() * monomial(myPolyRing, 1, indet(PPM(myPolyRing), index)), iter.myGetAnc());
        myInsertSetQ(triple);
        iter.IamSetNM(index);
      }
    }
  }
}


void JBSets::myJInsertWithoutProlong(JanetTree& JTree ,const std::list<JanetTriple>::iterator& ps)
{
  //like myJInsert skipping prolongations
  long i(0);
  JanetIterator iter(JTree);
  long n(NumIndets(myPolyRing));
  while( i < n)
  {
    while((iter.myCurrentDeg() < exponent(LPP(*(ps->myGetPolPtr())),i)) && (iter.myDisNextDeg()))
    {
      iter.myNextDeg();
    }
    if(iter.myCurrentDeg() < exponent(LPP(*(ps->myGetPolPtr())),i))
    {
      JanetTree tree(myOneLeafTree(ps, i, iter.myGetMonomial() * indet(PPM(myPolyRing),i)));
      iter.myConnectJanetTreeDeg(tree);
      i = n ;
    }
    else
    {
      if(exponent(LPP(*(ps->myGetPolPtr())),i) == iter.myCurrentDeg()  && (iter.myDisNextVar()))
      {
        bool NewVarNode(false);
        long DisNextVarNode(0);
        for(long j = i + 1; j < (iter.myDisNextVar() + iter.myCurrentVar()); ++j)
        {
          NewVarNode = (exponent(LPP(*(ps->myGetPolPtr())),j) != 0);
          if(NewVarNode)
          {
            DisNextVarNode = j - i;
            break;
          }
        }
        if(NewVarNode)
        {
          iter.mySetNextVar(DisNextVarNode);
        }
        i += iter.myNextVar();
      }
      else //exponent(LPP(ps.myGetPol()),i) > iter.myCurrentDeg()
      {
        if(exponent(LPP(*(ps->myGetPolPtr())),i) < iter.myCurrentDeg())
        {
          iter.myPrevDeg();
          iter.mySetNextDeg(exponent(LPP(*(ps->myGetPolPtr())),i) - iter.myCurrentDeg());
          iter.myNextDeg();
        }
        JanetTree tree(myOneLeafTree(ps, i + 1, iter.myGetMonomial()));
        iter.mySetNextVar(tree.myGetBeginVar() - iter.myCurrentVar());
        iter.myNextVar();
        iter.myConnectJanetTreeDeg(tree);
        i = n;
      }
    }
  }
}


void JBSets::myMinimizeSetT()
{
  long deg(StdDeg((*(mySetP.begin()))->myGetPol()));
  MultisetIterator StartIter(mySetT.begin());
  //searching the first element t in mySetT where deg >= deg(t)
  while((StartIter != mySetT.end()) && (deg >= StdDeg((*StartIter)->myGetPol())))
  {
    ++StartIter;
  }
  mySetQ.insert(StartIter, mySetT.end());
  mySetT.erase(StartIter, mySetT.end());
}


void JBSets::myTailReduceSetP(JanetTree& JTree)
{
  for(MultisetIterator iter = mySetP.begin(); iter != mySetP.end(); ++iter)
  {
    myJTailNormalForm(JTree,*iter);
  }
}

void JBSets::myInsertSetPInSetT()
{
  mySetT.insert(mySetP.begin(),mySetP.end());
  mySetP.clear();
}
