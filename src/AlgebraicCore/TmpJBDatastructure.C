#include "CoCoA/TmpJBDatastructure.H"
#include "CoCoA/utils.H"


using namespace CoCoA;

void JanetTriple::mySetNM(ConstRefPPMonoidElem var)
{
  long i;
  //check if var is only an indet
  if(IsIndet(i, var))
  {
    myNonMultiplicatives[i] = true;
  }
}

void JanetTriple::myDelNM(ConstRefPPMonoidElem var)
  {
    long  i;
    //check if var is only an indet
    if(IsIndet(i, var))
    {
      myNonMultiplicatives[i] = false;
    }
  }

bool JanetTriple::IamNM(ConstRefPPMonoidElem var)
{
  long i;
  //check if var is only an indet
  if(IsIndet(i, var))
  {
    return myNonMultiplicatives[i];
  }
  // if not return false
  return false;
}

std::vector<PPMonoidElem> JanetTriple::myGetNM()
{
  std::vector<PPMonoidElem> vec;
  long length = len(myNonMultiplicatives);
  for(long i = 0 ; i != length; ++i) //iterating over all variables
  {
    if(myNonMultiplicatives[i])
    {
      //vec.push_back(x_i)
      vec.push_back(indet(PPM(owner(myPolynom)),i)); 
    }
  }
  return vec;
}

void JanetHandle::myDecrUse()
{
  if(--*myUse == 0)
  {
    delete myPtr;
    delete myUse;
  }
}

void JanetTree::myAddAtBegin(JanetTree& tree)
{
  if (myBeginVar == tree.myBeginVar)
    {
      //if the current tree has the same beginning var as the new tree, then we cann
      //connect these trees easily
      JanetIterator iter(*this);
      iter.myConnectJanetTreeDeg(tree);
    }
  else
  {
    if(tree.myBeginDeg == 0)
    {
      //dangerours if tree.myBeginVar < myBeginVar we have a negative var_dis in the tree
      (*myArm.begin())->mySetDisNextVar(tree.myBeginVar - myBeginVar);
      (*myArm.begin())->mySetNextVarArm(*tree.myGetRootArm());
    }
    else
    {
      JanetIterator iter(*this);
      //create new var node such that this var matches the beginVar of tree
      iter.mySetNextVar(tree.myBeginVar - myBeginVar);
      iter.myNextVar();
      //then connect
      iter.myConnectJanetTreeDeg(tree);
    }
  }
}

 void JanetIterator::myChangeToLeafNode(const nodeData triple)
{
  JanetLeafNodeImpl leaf(triple);
  (*myCurIter) = JanetHandle(leaf);
  //if there is a arm which hanging around we must delete this arm:
  myCurArm->erase(++myCurIter,myCurArm->end());
  myCurIter = myCurArm->end();
  --myCurIter;
}

void JanetIterator::myChangeToInternalNode()
{
  //overwrite the current node
  (*myCurIter) = JanetHandle(JanetInternalNodeImpl());
}

long JanetIterator::myNextDeg()
{
  // check distance to next degree
  long dis = (*myCurIter)->myGetDisNextDeg();
  // if distance is not 0 then we add the distance to the current distance
  // and moving the vector iterator
  if(dis)
  {
    myMonomial[myCurVar] += dis;
    ++myCurIter;
    return dis;
  }
  // if not there is no moving -> moving distance is zero
  return 0;
}

long JanetIterator::myNextVar()
{
  // check distance to next var
  long dis((*myCurIter)->myGetDisNextVar());
  // if distance is not 0 then we 'initilize the degree at the current var'
  // and jumping to the next arm
  if(dis)
  {
    myCurVar += dis;
    CoCoA_ASSERT(myCurVar < NumIndets(myCurTree->myGetPolyRing()));
    myMonomial[myCurVar] = 0;
    myCurArm =(*myCurIter)->myNextVarArm();
    myCurIter = myCurArm->begin();
    return dis;
  }
  return 0;
}

void JanetIterator::mySetNextDeg(long dis)
{
  //nothing to do if dis == 0
  if(dis > 0)
  {
    //dis to next node
    long deg = (*myCurIter)->myGetDisNextDeg();
    if(deg)
    {
      // it there is a next deg we have to check if the new deg would be between 
      // the current deg and the following one
      if(deg > dis)
      {
        // if this is the case we add the new node between both nodes
        myCurArm->insert(++myCurIter,JanetHandle(JanetInternalNodeImpl()));
        (*(--myCurIter))->mySetDisNextDeg(deg - dis);
        (*(--myCurIter))->mySetDisNextDeg(dis);
      }
      else
      {
        // if not we erase the deg-direction starting with the node after current node
        // (in deg direction)
        myCurArm->erase((++myCurIter),myCurArm->end());
        myCurIter = myCurArm->end();
        (*(--myCurIter))->mySetDisNextDeg(dis);
        myCurArm->push_back(JanetHandle(JanetInternalNodeImpl()));
      }
    }
    else
    {
      //if there is no next deg we can append the new node easily
      if((*myCurIter)->IamLeafNode())
      {
        *myCurIter = JanetHandle(JanetInternalNodeImpl());
      }
      (*myCurIter)->mySetDisNextDeg(dis);
      myCurArm->push_back(JanetHandle(JanetInternalNodeImpl()));
    }
  }
}

void JanetIterator::mySetNextVar(long dis)
{
  //same as above but with var.
  if(dis)
  {
    long var((*myCurIter)->myGetDisNextVar());
    if(var)
    {
      if(var > dis)
      {
        JanetInternalNodeImpl node;
        std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
        node.mySetNextVarArm(*varArm);
        node.mySetDisNextVar(var - dis);
        varArm->erase(varArm->begin(),varArm->end());
        varArm->push_back(node);
        (*myCurIter)->mySetDisNextVar(dis);
      }
      else 
      {
        (*myCurIter)->mySetDisNextVar(dis);
        std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
        varArm->erase(varArm->begin(),varArm->end());
        varArm->push_back(JanetHandle(JanetInternalNodeImpl()));
      }
    }
    else
    {
      if((*myCurIter)->IamLeafNode())
      {
        *myCurIter = JanetHandle(JanetInternalNodeImpl());
      }
      (*myCurIter)->mySetDisNextVar(dis);
      std::list<JanetHandle>* varArm((*myCurIter)->myNextVarArm());
      varArm->push_back(JanetHandle(JanetInternalNodeImpl()));
    }
  }
}

void JanetIterator::mySetNextDeg(const nodeData& triple,long dis)
{
  // same as above, but no it is necessary that the current node is a leaf node
  CoCoA_ASSERT(!((*myCurIter)->IamLeafNode()));
  if(dis)
  {
    long deg = (*myCurIter)->myGetDisNextDeg();
    if(deg)
    {
      myCurArm->erase((++myCurIter),myCurArm->end());
      myCurIter = myCurArm->end();
      --myCurIter;
      myCurArm->push_back(JanetHandle(JanetLeafNodeImpl(triple)));
      (*myCurIter)->mySetDisNextDeg(dis);
    }
    else
    {
      if((*myCurIter)->IamLeafNode())
      {
        *myCurIter = JanetHandle(JanetLeafNodeImpl(triple));
      }
      (*myCurIter)->mySetDisNextDeg(dis);
      myCurArm->push_back(JanetHandle(JanetLeafNodeImpl(triple)));
    }
  }
}

void JanetIterator::myConnectJanetTreeDeg(JanetTree &tree)
{
  // check if current position and the beginning of the tree matches
  if((tree.myGetBeginDeg() == myMonomial[myCurVar]) && (tree.myGetBeginVar() == myCurVar))
  {
    // if yes everything is nice and we can combine the trees (we combine the last node
    // of the current tree and the first node of tree)
    // maybe here is a mistake because we don't check wether the current node is useless
    if((*myCurIter)->IamLeafNode())
    {  
      (*myCurIter) = JanetHandle(JanetInternalNodeImpl());
    }
    std::list<JanetHandle>* rootArm(tree.myGetRootArm());
    std::list<JanetHandle>::iterator tempIter(rootArm->begin());
    (*myCurIter)->mySetNextVarArm(*((*tempIter)->myNextVarArm()));
    (*myCurIter)->mySetDisNextDeg((*tempIter)->myGetDisNextDeg());
    std::list<JanetHandle>::iterator beg(rootArm->begin()), end(rootArm->end());
    for(std::list<JanetHandle>::iterator it = ++beg; it != end; ++it)
    {
      myCurArm->push_back(*it);
    }
  }
  else
  {
    // also simpling combining. But no we have a distande between the curNode and the beginning of
    // the new tree
    long extraDis(0);
    std::list<JanetHandle>* arm(tree.myGetRootArm());
    std::list<JanetHandle>::iterator beg(arm->begin()), end(arm->end());
    // checks if there is at the beginning a node in var direction and if the beginning node is a leaf node
    // if this isn't the case then we ignore the beginning node (to provide the minimality of the tree)
    if(!((*beg)->myGetDisNextVar()) && !(*beg)->IamLeafNode())
    {
      extraDis = (*beg)->myGetDisNextDeg();
      ++beg;
    }
    (*myCurIter)->mySetDisNextDeg(tree.myGetBeginDeg() + extraDis - myMonomial[myCurVar]);
    myCurArm->erase(++myCurIter,myCurArm->end());
    myCurIter = myCurArm->end();
    --myCurIter;
    for(std::list<JanetHandle>::iterator it = beg; it != end; ++it)
    {
      myCurArm->push_back(*it);
    }
  }
}

void JanetIterator::myDelNextVar()
{
  // is this the last node in var direction?
  if((*myCurIter)->myGetDisNextVar())
  {
    // if not we delete the complete deg direction
    std::list<JanetHandle>::iterator it = ((*myCurIter)->myNextVarArm())->begin();
    long dis((*it)->myGetDisNextVar());
    if(dis)
    {
      // if this isn't the last node in var direction we replace the current node
      // with the next var-node
      (*myCurIter)->mySetDisNextVar((*myCurIter)->myGetDisNextVar() + dis);
      (*myCurIter)->mySetNextVarArm(*((*it)->myNextVarArm()));
    }
    else
    {
      (*myCurIter)->mySetDisNextVar(0);
      ((*myCurIter)->myNextVarArm())->clear();
    }
  }
}

const JanetTriple& JanetIterator::myGetTriple() const
{
  CoCoA_ASSERT((*myCurIter)->myGetTriplePtr() != 0); 
  return *((*myCurIter)->myGetTriplePtr());
} 

bool JanetIterator::IamSetTriple(const JanetTriple& newTriple)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    *triple = newTriple;
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamSetPol(ConstRefRingElem newPol)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    triple->mySetPol(newPol);
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamSetAnc(ConstRefPPMonoidElem newAnc)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    triple->mySetAnc(newAnc);
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamSetNM(ConstRefPPMonoidElem var)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    triple->mySetNM(var);
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamSetNM(const long &index)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    triple->mySetNM(index);
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamDelNM(ConstRefPPMonoidElem var)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr();
    triple->myDelNM(var);
    return true;
  }
  else
  {
    return false;
  }
}

bool JanetIterator::IamDelNM(const long &index)
{
  // if this isn't a leaf node return false
  if((*myCurIter)->IamLeafNode())
  {
    JanetTriple* triple = (*myCurIter)->myGetTriplePtr(); 
    triple->myDelNM(index);
    return true;
  }
  else
  {
    return false;
  }
}

void JanetIterator::myReturnToBegin()
{
  // everything to the beginning
  myCurArm = myCurTree->myGetRootArm();
  myCurIter =myCurArm->begin();
  myCurVar = myCurTree->myGetBeginVar();
  myMonomial = std::vector<long>(NumIndets(myCurTree->myGetPolyRing()));
  myMonomial[myCurVar] = myCurTree->myGetBeginDeg();
}

long JanetIterator::myPrevDeg()
{
  // returns to the previous deg. If we are already at the beginning do nothing
  if( myCurIter == myCurArm->begin())
  {
    return 0;
  }
  else
  {
    --myCurIter;
    myMonomial[myCurVar] -= (*myCurIter)->myGetDisNextDeg();
    return (*myCurIter)->myGetDisNextDeg();
  }
}

void JanetTriple::myClearNM()
{
  // delete all nonmultiplicative variables -> all entries false
  for(std::vector<bool>::iterator iter =myNonMultiplicatives.begin(); iter != myNonMultiplicatives.end(); ++iter)
  {
    *iter = false;
  }
}

JanetTriple* JanetTree::myJDivisor(ConstRefPPMonoidElem w) const
{
  // straight forward implementation
  std::vector<long> expv(NumIndets(PPM(myPolyRing)));
  PPM(myPolyRing)->myExponents(expv, raw(w));
  std::list<JanetHandle>::const_iterator iter(myArm.begin());
  long curVar(myBeginVar);
  long curDeg(myBeginDeg);
  long numIndets(NumIndets(owner(w)));
  while(numIndets >= curVar)
  {
    while(((*iter)->myGetDisNextDeg()) && (curDeg < expv[curVar]))//exponent(w,curVar)))
    {
      if(curDeg + (*iter)->myGetDisNextDeg() > expv[curVar])
      {
        return 0;
      }
      curDeg += (*iter)->myGetDisNextDeg();
      ++iter;
    }
    if((*iter)->myGetDisNextVar())
    {
      curVar += (*iter)->myGetDisNextVar();
      iter = ((*iter)->myNextVarArm())->begin();
      curDeg = 0;
    }
    else
    {
      if((*iter)->myGetDisNextDeg())
      {
        return 0;
      }
      else
      {
        return (*iter)->myGetTriplePtr();
      }
    }
  }
  return 0; //only of safety. A return command going to be used in the while loop
}

JanetIterator JanetIterator::myGotoHighestNode()
{
  // goto the highest node: If we are in the highest node go to the next var if possible
  // if not return the current node
  JanetIterator iter(*this);
  while(true)
  {
    while(iter.myNextDeg() != 0)
    {
    }
    if(iter.myDisNextVar() == 0)
    {
      return iter;
    }
    else
    {
      iter.myNextVar();
    }
  }
  return iter;    
}
