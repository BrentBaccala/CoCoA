//   Copyright (c) 2010-2013 Giovanni Lagorio, John Abbott, Anna M. Bigatti
//
//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "BuiltInFunctions.H"

using namespace std;
using namespace boost;
using namespace boost::iostreams;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;

namespace CoCoA {
namespace InterpreterNS {

//----------------------------------------------------------------------
// builtIns is essentially a global variable but with better controlled initialization.
std::vector<NameFunPair>& builtIns()
{
  static std::vector<NameFunPair> GlobalVar;
  return GlobalVar;
}

namespace{ // anonymous  =============== IsVectorBigInt =================
  bool IsVectorBigInt(std::vector<BigInt>& BigIntVec, const intrusive_ptr<LIST> l)
  {
    vector<BigInt> v;
    LIST::ContainerType::size_type size = l->size();
    for (unsigned long i=0; i<size; ++i)
      if (const boost::intrusive_ptr<INT> n = boost::dynamic_pointer_cast<INT>(l->getValue(i)))
        v.push_back(n->theBigInt);
      else
        return false;
    swap(v, BigIntVec);
    return true;
  }  
} // anonymous namespace




DECLARE_STD_BUILTIN_FUNCTION(NumEvalUPoly, 2)
{ // JAA
  const BigRat q = runtimeEnv->evalArgAs<RAT>(ARG(1))->theBigRat;
  intrusive_ptr<LIST> CoeffList = runtimeEnv->evalArgAs<LIST>(ARG(0));
  LIST::ContainerType::size_type size = CoeffList->size();
  vector<BigInt> C; C.reserve(size);
  if (!IsVectorBigInt(C, CoeffList)) CoCoA_ERROR("Bad CoeffList","NumEvalUPoly");
  if (IsZero(q)) return new RAT(BigRat(C.front(),1));
  BigInt numer;
  HornerRecursiveIterQQ2(numer, C, num(q), den(q));
  return new RAT(BigRat(numer, power(den(q),size-1)));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(not, 1) {
  intrusive_ptr<BOOL> arg = runtimeEnv->evalArgAs<BOOL>(ARG(0));
  return Value::from(!(arg->theBool));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetStackSize, 1) {
  intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(0));
  const BigInt NewSize = N->theBigInt;
  long n;
  if (NewSize < 2 || !IsConvertible(n, NewSize)) throw RuntimeException("Ridiculous stack size", ARG(0).exp);
  return new INT(runtimeEnv->ResizeStack(n, ARG(0).exp));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsInteger, 2) {
	intrusive_ptr<LeftValue> resultRef = intrusive_ptr_cast<LeftValue>(runtimeEnv->evalArg(ARG(0), RuntimeEnvironment::EVAL_BY_REF));
	intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
	BigInt N;
	if (!IsInteger(N, poly->theRingElem)) return Value::from(false);
  if (resultRef->assignmentNeedsOwnership()) resultRef->obtainOwnership();
  const intrusive_ptr<const Expression> resultRefExp = ARG(0).exp;
  resultRef->assign(new INT(N), resultRefExp->getBegin(), resultRefExp->getEnd(), runtimeEnv);
	return Value::from(true);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsRational, 2) {
	intrusive_ptr<LeftValue> resultRef = intrusive_ptr_cast<LeftValue>(runtimeEnv->evalArg(ARG(0), RuntimeEnvironment::EVAL_BY_REF));
	intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
	BigRat qq;
	bool isRational(IsRational(qq, poly->theRingElem));
	if (isRational) {
		if (resultRef->assignmentNeedsOwnership()) resultRef->obtainOwnership();
		const intrusive_ptr<const Expression> resultRefExp = ARG(0).exp;
		resultRef->assign(new RAT(qq), resultRefExp->getBegin(), resultRefExp->getEnd(), runtimeEnv);
	}
	return Value::from(isRational);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(aliases, 0) {
	runtimeEnv->topLevelAliases->dump(runtimeEnv->getOutputStream());
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(pause, 0) {
#ifdef C5IDE
	runtimeEnv->interpreter->singleStepExecution = true;
#else // #ifdef C5IDE
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
#endif // #ifdef C5IDE
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

#ifdef C5IDE
DECLARE_STD_BUILTIN_FUNCTION(sleep, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(arg);
	long l;
	if (!IsConvertible(l, N->theBigInt) || l<0 || l>1000)
		throw RuntimeException("The argument must be in the range [0..1000]", arg.exp);
	this_thread::sleep(boost::posix_time::seconds(static_cast<int>(l)));
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(debug, 1) {
	const Argument &arg = ARG(0);
	string filename = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	try {
		runtimeEnv->interpreter->singleStepExecution = true;
		runtimeEnv->interpreter->readAndExecute(filename, false, false);
	} catch (const RuntimeException &) {
		throw;
	} catch (const BaseException &be) {
		throw RuntimeException(be.reason, invocationExpression->getBegin(), invocationExpression->getEnd());
	}
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION
#endif // #ifdef C5IDE

#if 0  // 2 functions for accessing zip files, currently disabled
DECLARE_STD_BUILTIN_FUNCTION(ZipFileList, 1) {
	const Argument &arg = ARG(0);
	const string filename = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	int err;
	struct zip *zf = ::zip_open(filename.c_str(), ZIP_CHECKCONS, &err);
	if (!zf)
		throw RuntimeException("Cannot read Zip file \""+filename+"\"", arg.exp);
	intrusive_ptr<LIST> result(new LIST);
	const int numEntries = zip_get_num_files(zf);
	for(int a=0; a<numEntries; ++a)
		result->addValue(new STRING(zip_get_name(zf, a, ZIP_FL_UNCHANGED)));
	::zip_close(zf);
	return result;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ZipRead, 2) {
	const Argument &argFilename = ARG(0);
	const Argument &argEntryName = ARG(1);
	const string filename = runtimeEnv->evalArgAs<STRING>(argFilename)->theString;
	const string entryName = runtimeEnv->evalArgAs<STRING>(argEntryName)->theString;
	int err;
	struct zip *zf = ::zip_open(filename.c_str(), ZIP_CHECKCONS, &err);
	if (!zf)
		throw RuntimeException("Cannot read Zip file \""+filename+"\"", argFilename.exp);
	const int index = zip_name_locate(zf, entryName.c_str(), /* flags */ 0);
	if (index<0) {
		::zip_close(zf);
		throw RuntimeException("Cannot find an entry named \""+entryName+" inside the Zip file \""+filename+"\"", argEntryName.exp);
	}
	struct zip_file *entry = ::zip_fopen_index(zf, index, /* flags */ 0);
	if (!entry) {
		::zip_close(zf);
		throw RuntimeException("Error reading the entry named \""+entryName+" inside the Zip file \""+filename+"\"", argEntryName.exp);
	}
	string result;
	for(;;) {
		char buf[1024];
		int nRead = ::zip_fread(entry, buf, sizeof(buf));
		if (nRead<=0) {
			::zip_fclose(entry);
			::zip_close(zf);
			return new STRING(result);
		}
		result.append(buf, buf+nRead);
	}
}
END_STD_BUILTIN_FUNCTION
#endif

DECLARE_STD_BUILTIN_FUNCTION(len, 1) {
	intrusive_ptr<const RightValue> v=runtimeEnv->evalArgAs<RightValue>(ARG(0));
	if (dynamic_pointer_cast<const MAT>(v))
    throw RuntimeException("len(matrix) not allowed, use NumRows(matrix) instead", ARG(0).exp);
	if (dynamic_pointer_cast<const MatrixRowValue>(v))
    throw RuntimeException("len(MatrixRow) not allowed, use NumCols(matrix) instead", ARG(0).exp);
	int which;
	v = runtimeEnv->evalArgAsT1orT2orT3<STRING, LIST, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(BigInt(RefTo<string>(v).length()));
	case 2: return Value::from(BigInt(intrusive_ptr_cast<const LIST>(v)->size()));
  case 3: return Value::from(BigInt(NumTerms(RefTo<RingElem>(v))));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_ARITYCHECK_FUNCTION(sum) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(sum) {
	invocationExpression->checkNumberOfArgs(1,2);
	const Argument &a0 = ARG(0);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(a0);
	LIST::ContainerType::size_type size = list->size();
  intrusive_ptr<RightValue> result = INT::zero;
  if (invocationExpression->args.size()==1)
  {
    if (size==0) return result;
    result = list->getValue(--size);
  }
  else
    result = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(LIST::ContainerType::size_type a=size; a!=0; /**/)
    try	{
			result = runtimeEnv->binaryOperatorDispatch(list->getValue(--a), result, RuntimeEnvironment::opPlusMap, a0.exp->getBegin(), a0.exp->getEnd());
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		throw RuntimeException("Some elements have incompatible types for the sum", a0.exp);
	}
	return result;
}

DECLARE_ARITYCHECK_FUNCTION(product) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(product) {
	invocationExpression->checkNumberOfArgs(1,2);
	const Argument &a0 = ARG(0);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(a0);
	LIST::ContainerType::size_type size = list->size();
	intrusive_ptr<RightValue> result = INT::one;
  if (invocationExpression->args.size()==1)
  {
    if (size==0) return result;
    result = list->getValue(--size);
  }
  else
    result = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(LIST::ContainerType::size_type a=size; a!=0; /**/)
		try	{
			result = runtimeEnv->binaryOperatorDispatch(list->getValue(--a), result, RuntimeEnvironment::opStarMap, a0.exp->getBegin(), a0.exp->getEnd());
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		throw RuntimeException("Some elements have incompatible types for the product", a0.exp);
	}
	return result;
}

DECLARE_ARITYCHECK_FUNCTION(NewList) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(NewList) {
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(0));
  long n;
  if (!IsConvertible(n, N->theBigInt) || n<0)
     throw RuntimeException("invalid length", ARG(0).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
  intrusive_ptr<RightValue> val = INT::zero;
  if (invocationExpression->args.size()==2)
    val = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(long i=0; i<n; ++i)  returnValue->addValue(val);
	return returnValue;
}

DECLARE_ARITYCHECK_FUNCTION(first) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(first) {
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
  if (invocationExpression->args.size()==1)
  {
    if (list->size() == 0) throw RuntimeException("list is empty", ARG(0).exp);
    return list->getValue(0);
  }
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, N->theBigInt) || n<0)
    throw RuntimeException("invalid length", ARG(1).exp);
  if ((unsigned long)n > list->size())
    throw RuntimeException("value greater than list length", ARG(1).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
	for(long i=0; i<n; ++i)  returnValue->addValue(list->getValue(i));
	return returnValue;
}

DECLARE_ARITYCHECK_FUNCTION(last) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(last) {
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
  ///////// here
  LIST::ContainerType::size_type s = list->size();
  if (invocationExpression->args.size()==1) 
  if (invocationExpression->args.size()==1)
  {
    if (list->size() == 0) throw RuntimeException("list is empty", ARG(0).exp);
    return list->getValue(s-1);
  }
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, N->theBigInt) || n<0)
    throw RuntimeException("invalid length", ARG(1).exp);
  if ((unsigned long)n>s)
    throw RuntimeException("value greater than list length", ARG(1).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
	for(unsigned long i=s-n; i<s; ++i)  returnValue->addValue(list->getValue(i));
	return returnValue;
}

DECLARE_STD_BUILTIN_FUNCTION(ScalarProduct, 2) {
	const Argument &arg1 = ARG(0);
	intrusive_ptr<LIST> l1 = runtimeEnv->evalArgAs<LIST>(arg1);
	const Argument &arg2 = ARG(1);
	intrusive_ptr<LIST> l2 = runtimeEnv->evalArgAs<LIST>(arg2);
	const LIST::ContainerType::size_type size = l1->size();
	const CharPointer &begin = arg1.exp->getBegin();
	const CharPointer &end = arg2.exp->getEnd();
	if (l2->size()!=size)
		throw RuntimeException("The arguments must be lists of the same size", begin, end);
	intrusive_ptr<RightValue> result = INT::zero;
	try	{
		for(LIST::ContainerType::size_type a=0; a<size; ++a)
			result = runtimeEnv->binaryOperatorDispatch(
						result,
						runtimeEnv->binaryOperatorDispatch(l1->getValue(a), l2->getValue(a), RuntimeEnvironment::opStarMap, begin, end),
						RuntimeEnvironment::opPlusMap,
						begin,
						end);
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		throw RuntimeException("Some elements have incompatible types for the scalar product", begin, end);
	}
	return result;
}
END_STD_BUILTIN_FUNCTION

DECLARE_ARITYCHECK_FUNCTION(concat) { (void)nArg; return true; }

DECLARE_BUILTIN_FUNCTION(concat) {
	intrusive_ptr<LIST> returnValue = new LIST();
	BOOST_FOREACH(const Argument &arg, invocationExpression->args) {
		intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(arg);
		int size = list->size();
		for(int a=0; a<size; ++a)
			returnValue->addValue(list->getValue(a));
	}
	return returnValue;
}

DECLARE_STD_BUILTIN_FUNCTION(reversed, 1) {
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
	intrusive_ptr<LIST> returnValue = new LIST();
	for(LIST::ContainerType::size_type a=list->size(); a!=0; /**/)
		returnValue->addValue(list->getValue(--a));
	return returnValue;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(reverse, 1) {
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	const LIST::ContainerType::size_type size = list->size();
	if (size) {
		LIST::ContainerType::size_type left = 0, right = size-1;
		while(left < right) {
			const intrusive_ptr<RightValue> tmp(list->getValue(left));
			list->setValue(left++, list->getValue(right));
			list->setValue(right--, tmp);
		}
	}
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(deg) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(deg) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(deg(a->theRingElem));
  intrusive_ptr<RINGELEM> b = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  long i;
  if (!IsIndet(i,b->theRingElem))
    throw RuntimeException("not an indet", ARG(1).exp);
  return Value::from(deg(a->theRingElem, i));
}

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(coefficients) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(coefficients) {  // AMB
	invocationExpression->checkNumberOfArgs(1,2);// variable number of args
	intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(coefficients_forC5(a->theRingElem));
	vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
  vector<RingElem> res;
  for (long i=0; i<len(v); ++i)
    res.push_back(CoeffOfTerm_forC5(a->theRingElem, v[i]));
  return Value::from(res);
}

// // variable number of args
// DECLARE_ARITYCHECK_FUNCTION(NewMat) { return (2<=nArg) && (nArg<=3); }
// DECLARE_BUILTIN_FUNCTION(NewMat) {
// 	invocationExpression->checkNumberOfArgs(2,3);
//   if (invocationExpression->args.size()==2)
//     throw RuntimeException("NewMat(NR,NC) not allowed, use NewMat(R:RING,NR:INT,NC:INT) instead", ARG(0).exp);
//   intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
//   intrusive_ptr<INT> nr = runtimeEnv->evalArgAs<INT>(ARG(1));
//   intrusive_ptr<INT> nc = runtimeEnv->evalArgAs<INT>(ARG(2));
//   return Value::from(ZeroMat(R->theRing, ConvertTo<long>(nr->theBigInt), ConvertTo<long>(nc->theBigInt)));
// }

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(jacobian) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(jacobian) {  // AMB
	invocationExpression->checkNumberOfArgs(1,2);// variable number of args
	vector<RingElem> v0 = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  if (invocationExpression->args.size()==1)
  {
    if (v0.size()==0) throw RuntimeException("Empty list", ARG(0).exp);
    return new MAT(jacobian(v0, indets(owner(v0[0]))));
  }
	vector<RingElem> v1 = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
  return new MAT(jacobian(v0, v1));
}

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(CompleteToOrd) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(CompleteToOrd) {  // AMB
	invocationExpression->checkNumberOfArgs(1,2);// variable number of args
	intrusive_ptr<const MAT> M0 = runtimeEnv->evalArgAs<const MAT>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(NewMatCompleteOrd(M0->theMatrix));
	intrusive_ptr<const MAT> M1 = runtimeEnv->evalArgAs<const MAT>(ARG(1));
  return Value::from(NewMatCompleteOrd(M0->theMatrix, M1->theMatrix));
}


DECLARE_STD_BUILTIN_FUNCTION(HomogElimMat, 2) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
	vector<BigInt> v1 = runtimeEnv->evalArgAsListOf<INT>(ARG(1));
  return Value::from(HomogElimMat_forC5(M->theMatrix, v1));
}
END_STD_BUILTIN_FUNCTION


// DECLARE_ARITYCHECK_FUNCTION(GrammSchmidtRows) { return (1<=nArg) && (nArg<=2); }
// DECLARE_BUILTIN_FUNCTION(GrammSchmidtRows) {  // AMB
// 	invocationExpression->checkNumberOfArgs(1,2);
// 	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
//   if (invocationExpression->args.size()==1)
//   {
//     GrammSchmidtRows(M->theMatrix);
//     return VoidValue::theInstance;
//   }
// 	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
//   long n;
//   if (!IsConvertible(n, N->theBigInt))
//     throw RuntimeException("invalid row index", ARG(1).exp);
//   GrammSchmidtRows(M->theMatrix, n);
//   return VoidValue::theInstance;
// }


DECLARE_STD_BUILTIN_FUNCTION(append, 2) {
	intrusive_ptr<RightValue> elem = runtimeEnv->evalArgAs<RightValue>(ARG(1));
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	list->addValue(elem);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(remove, 2) {
	const Argument indexArg = ARG(1);
	intrusive_ptr<RightValue> shouldBeAnIndex = runtimeEnv->evalArgAs<RightValue>(indexArg);
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	LIST::ContainerType::size_type index = list->checkIndex(shouldBeAnIndex, indexArg.exp);
	list->removeValue(index);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(insert, 3) {
	const Argument &indexArg = ARG(1);
	intrusive_ptr<RightValue> shouldBeAnIndex = runtimeEnv->evalArgAs<RightValue>(indexArg);
	intrusive_ptr<RightValue> elem = runtimeEnv->evalArgAs<RightValue>(ARG(2));
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	LIST::ContainerType::size_type index = list->checkIndex(shouldBeAnIndex, indexArg.exp);
	list->insertValue(index, elem);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(div, 2) {
	intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
	BigInt result;
	div(result, a->theBigInt, b->theBigInt);
	return new INT(result);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(mod, 2) {
	intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
	BigInt result;
	mod(result, a->theBigInt, b->theBigInt);
	return new INT(result);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NumDigits, 2) {
	intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
	long base;
	if (!IsConvertible(base, b->theBigInt) || base<2 || base>36)
		throw RuntimeException("Base must be in the range 2..36", ARG(1).exp);
	return new INT(NumDigits(a->theBigInt, base));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ContFracToRat, 1) { // JAA
  const vector<BigInt> CFQuots = runtimeEnv->evalArgAsListOf<INT>(ARG(0));
    ContFracApproximant ans;
    const long n = len(CFQuots);
    for (long i=0; i < n; ++i)
    {
      if (i > 0 && sign(CFQuots[i])<=0)
        throw RuntimeException("Cont frac quotients must  be positive", ARG(0).exp);
      ans.myAppendQuot(CFQuots[i]);
    }
    return Value::from(ans.myRational());
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CRT, 4) {
	intrusive_ptr<INT> R1 = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> M1 = runtimeEnv->evalArgAs<INT>(ARG(1));
	intrusive_ptr<INT> R2 = runtimeEnv->evalArgAs<INT>(ARG(2));
	intrusive_ptr<INT> M2 = runtimeEnv->evalArgAs<INT>(ARG(3));
        CRTMill CRT;
        CRT.myAddInfo(R1->theBigInt, M1->theBigInt);
        CRT.myAddInfo(R2->theBigInt, M2->theBigInt);

        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setField("residue", Value::from(residue(CRT)));
        ans->setField("modulus", Value::from(modulus(CRT)));
        return ans;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(RatReconstructWithBounds, 5) {
	intrusive_ptr<INT> E = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> P = runtimeEnv->evalArgAs<INT>(ARG(1));
	intrusive_ptr<INT> Q = runtimeEnv->evalArgAs<INT>(ARG(2));
	const vector<BigInt> ResList = runtimeEnv->evalArgAsListOf<INT>(ARG(3));
	const vector<BigInt> ModList = runtimeEnv->evalArgAsListOf<INT>(ARG(4));

        const long NumRes = len(ResList);
        vector<long> res(NumRes);
        for (long i=0; i < NumRes; ++i)
          res[i] = ConvertTo<long>(ResList[i]);
        const long NumMod = len(ModList);
        vector<long> mod(NumMod);
        for (long i=0; i < NumMod; ++i)
          mod[i] = ConvertTo<long>(ModList[i]);
        const long e = ConvertTo<long>(E->theBigInt);
        const BigRat result = RatReconstructWithBounds(e, P->theBigInt, Q->theBigInt, res, mod);
        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setField("failed", Value::from(den(result) == 0));
        if (den(result) != 0)
          ans->setField("ReconstructedRat", Value::from(result));
        return ans;
}
END_STD_BUILTIN_FUNCTION

DECLARE_ARITYCHECK_FUNCTION(RatReconstructByContFrac) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(RatReconstructByContFrac) {
  invocationExpression->checkNumberOfArgs(2,3);
	intrusive_ptr<INT> X = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> M = runtimeEnv->evalArgAs<INT>(ARG(1));

        BigInt threshold; // Determine threshold: 0 means use default value
        if (invocationExpression->args.size()==3)
        {
          intrusive_ptr<INT> thresh = runtimeEnv->evalArgAs<INT>(ARG(2));
//          if (thresh->theBigInt < 0) throw RuntimeException("Threshold must be >= 0", ARG(2).exp);
          threshold = thresh->theBigInt;
        }

        RatReconstructByContFrac reconstructor(threshold);
        reconstructor.myAddInfo(X->theBigInt, M->theBigInt);

        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setField("failed", Value::from(!IsConvincing(reconstructor)));
        if (IsConvincing(reconstructor))
          ans->setField("ReconstructedRat", Value::from(ReconstructedRat(reconstructor)));
        return ans;
}

DECLARE_ARITYCHECK_FUNCTION(RatReconstructByLattice) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(RatReconstructByLattice) {
  invocationExpression->checkNumberOfArgs(2,3);
	intrusive_ptr<INT> X = runtimeEnv->evalArgAs<INT>(ARG(0));
	intrusive_ptr<INT> M = runtimeEnv->evalArgAs<INT>(ARG(1));

        BigInt SafetyFactor; // Determine threshold: 0 means use default value
        if (invocationExpression->args.size()==3)
        {
          intrusive_ptr<INT> safety = runtimeEnv->evalArgAs<INT>(ARG(2));
//          if (safety->theBigInt < 0) throw RuntimeException("SafetyFactor must be >= 0", ARG(2).exp);
          SafetyFactor = safety->theBigInt;
        }

        RatReconstructByLattice reconstructor(SafetyFactor);
        reconstructor.myAddInfo(X->theBigInt, M->theBigInt);

        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setField("failed", Value::from(!IsConvincing(reconstructor)));
        if (IsConvincing(reconstructor))
          ans->setField("ReconstructedRat", Value::from(ReconstructedRat(reconstructor)));
        return ans;
}

DECLARE_STD_BUILTIN_FUNCTION(TmpMantissaAndExponent2, 1)
{
  intrusive_ptr<RINGELEM> tmp = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  const RingElem& x = tmp->theRingElem;

  if (!IsRingTwinFloat(owner(x))) throw RuntimeException("Arg must belong to a Twinfloat ring", ARG(0).exp);
  MantExp2 ME = MantissaAndExponent2(x);
  ME.myMantissa *= ME.mySign;

  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setField("mantissa", Value::from(ME.myMantissa));
  ans->setField("exponent", Value::from(ME.myExponent));
  ans->setField("NumDigits", Value::from(ME.myNumDigits));
  return ans;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(MantissaAndExponent2, 2)
{
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, Ndigits->theBigInt) || n < 1)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  MantExp2 ME = MantissaAndExponent2(q, n);
  ME.myMantissa *= ME.mySign;

  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setField("mantissa", Value::from(ME.myMantissa));
  ans->setField("exponent", Value::from(ME.myExponent));
  ans->setField("NumDigits", Value::from(ME.myNumDigits));
  return ans;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(MantissaAndExponent10, 2)
{
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, Ndigits->theBigInt) || n < 1)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  MantExp10 ME = MantissaAndExponent10(q, n);
  ME.myMantissa *= ME.mySign;

  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setField("mantissa", Value::from(ME.myMantissa));
  ans->setField("exponent", Value::from(ME.myExponent));
  ans->setField("NumDigits", Value::from(ME.myNumDigits));
  return ans;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(FloatApprox, 2)
{
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  intrusive_ptr<INT> Nbits = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, Nbits->theBigInt) || n < 2)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(FloatApprox(q, n));
}
END_STD_BUILTIN_FUNCTION

DECLARE_ARITYCHECK_FUNCTION(FloatStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(FloatStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  if (invocationExpression->args.size()==1)
    return Value::from(FloatStr(q));
  long n;
  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (!IsConvertible(n, Ndigits->theBigInt) || n < 1)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(FloatStr(q, n));
//???  return Value::from(ScientificStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}

DECLARE_ARITYCHECK_FUNCTION(ScientificStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(ScientificStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  if (invocationExpression->args.size()==1)
    return Value::from(ScientificStr(q));
  long n;
  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (!IsConvertible(n, Ndigits->theBigInt) || n < 1)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(ScientificStr(q, n));
//???  return Value::from(ScientificStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}

DECLARE_ARITYCHECK_FUNCTION(DecimalStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(DecimalStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }

  if (invocationExpression->args.size()==1)
    return Value::from(DecimalStr(q));
  long n;
  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (!IsConvertible(n, Ndigits->theBigInt) || n < 1)
    throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(DecimalStr(q, n));
//???  return Value::from(DecimalStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}

DECLARE_STD_BUILTIN_FUNCTION(sign, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return INT::fromInt(sign(RefTo<BigInt>(v)));
  case 2: return INT::fromInt(sign(RefTo<BigRat>(v)));
  case 3: return INT::fromInt(sign(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsZero, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5orT6orT7<INT, RAT, RINGELEM, MODULEELEM, IDEAL, MAT, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsZero(RefTo<BigInt>(v)));
  case 2: return Value::from(IsZero(RefTo<BigRat>(v)));
  case 3: return Value::from(IsZero(RefTo<RingElem>(v)));
  case 4: return Value::from(IsZero(RefTo<ModuleElem>(v)));
  case 5: return Value::from(IsZero(RefTo<ideal>(v)));
  case 6: return Value::from(IsZero(RefTo<matrix>(v)));
  case 7: return Value::from(IsZero(RefTo<module>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsOne, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsOne(RefTo<BigInt>(v)));
  case 2: return Value::from(IsOne(RefTo<BigRat>(v)));
  case 3: return Value::from(IsOne(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsMinusOne, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsMinusOne(RefTo<BigInt>(v)));
  case 2: return Value::from(IsMinusOne(RefTo<BigRat>(v)));
  case 3: return Value::from(IsMinusOne(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NumTerms, 1) { // AMB
	return Value::from(BigInt(NumTerms(runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem)));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(binomial, 2) {
	int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, INT>(ARG(0), which);
  intrusive_ptr<INT> n = runtimeEnv->evalArgAs<INT>(ARG(1));
  switch (which) {
  case 1: return Value::from(binomial(RefTo<RingElem>(x), n->theBigInt));
  case 2: return Value::from(binomial(RefTo<BigInt>(x), n->theBigInt));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(zero, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RING, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(zero(RefTo<ring>(x)));
  case 2: return Value::from(zero(RefTo<module>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(gens, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(gens(RefTo<ideal>(x)));
  case 2: return Value::from(gens(RefTo<module>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsContained, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsContained(RefTo<ideal>(x), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  case 2: return Value::from(IsContained(RefTo<module>(x), runtimeEnv->evalArgAs<MODULE>(ARG(1))->theModule));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsElem, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsElem(RefTo<RingElem>(x), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  case 2: return Value::from(IsElem(RefTo<ModuleElem>(x), runtimeEnv->evalArgAs<MODULE>(ARG(1))->theModule));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GBasis, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(GBasis(RefTo<ideal>(x)));
  case 2: return Value::from(TidyGens(RefTo<module>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(MinGens, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(MinGens(RefTo<ideal>(x)));
  case 2: return Value::from(MinGens(RefTo<module>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(syz, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<IDEAL, MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: CoCoA_ERROR(ERR::NYI, "minimal syzigies");
    return Value::from(NewFreeModule(RingZZ(),1));
  case 2: CoCoA_ERROR(ERR::NYI, "minimal syzigies");
    return Value::from(NewFreeModule(RingZZ(),1));
  case 3: {
//     if ()
//     {
//       vector<ModuleElem> x = runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(0));
//       return Value::from(SyzOfGens(submodule(x)));
//     }
//     else
      vector<RingElem> x = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
      return Value::from(SyzOfGens(ideal(x)));
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(SyzOfGens) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(SyzOfGens) {
	invocationExpression->checkNumberOfArgs(1,2);
  long n = invocationExpression->args.size();
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(n-1), which);
  if (n==1)
    switch (which) {
    case 1: return Value::from(SyzOfGens(RefTo<ideal>(x)));
    case 2: return Value::from(SyzOfGens(RefTo<module>(x)));
    default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
    }
  intrusive_ptr<MODULE> M = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  switch (which) {
  case 1: return Value::from(SyzOfGens(M->theModule, RefTo<ideal>(x)));
  default:return Value::from(SyzOfGens(M->theModule, RefTo<module>(x)));
  }
}
// END_STD_BUILTIN_FUNCTION no: variable number of args

DECLARE_STD_BUILTIN_FUNCTION(RingElem, 2) {
	intrusive_ptr<const RING> R(runtimeEnv->evalArgAs<RING>(ARG(0)));
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<INT, RAT, STRING, RINGELEM, LIST>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(RingElem(R->theRing, RefTo<BigInt>(v)));
  case 2: return Value::from(RingElem(R->theRing, RefTo<BigRat>(v)));
  case 3: return Value::from(RingElem(R->theRing, symbol(RefTo<string>(v))));
  case 4: return Value::from(RingElem(R->theRing, RefTo<RingElem>(v)));
  case 5: {// Expecting list of [STRING, INT, INT, ..., INT] specifying a symbol
    const intrusive_ptr<const LIST> l = dynamic_pointer_cast<const LIST>(v);
    LIST::ContainerType::size_type ListLen = l->size();
    string SymbolHead;
    if (ListLen == 0) throw RuntimeException("Non-empty list expected", ARG(1).exp);
    if ( const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(l->getValue(0)))
      SymbolHead = s->theString;
    else
      throw RuntimeException("String expected as first entry in list", ARG(1).exp);
    std::vector<BigInt> indices;
    for(LIST::ContainerType::size_type a=1; a<ListLen; ++a)
    {
      if (const boost::intrusive_ptr<INT> elem = boost::dynamic_pointer_cast<INT>(l->getValue(a)))
        indices.push_back(theValue(elem));
      else
        throw RuntimeException("Integer expected for symbol index", ARG(1).exp);
    }
    return Value::from(RingElem(R->theRing, symbol(SymbolHead, VectorLong(indices, "RingElem"))));
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(homog, 2) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM, MODULEELEM, IDEAL, MODULE, LIST>(ARG(0), which);
	intrusive_ptr<const RINGELEM> h(runtimeEnv->evalArgAs<RINGELEM>(ARG(1)));
  switch (which) {
  case 1: return Value::from(homog(RefTo<RingElem>(v), h->theRingElem));
  case 2: return Value::from(homog(RefTo<ModuleElem>(v), h->theRingElem));
  case 3: return Value::from(homog(RefTo<ideal>(v), h->theRingElem));
  case 4: CoCoA_ERROR(ERR::NYI, "homog(module)");
//return Value::from(homog(intrusive_ptr_cast<MODULE>(v)->theModule, h->theRingElem));
  case 5: {
    vector<RingElem> v1 = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
    intrusive_ptr<LIST> returnValue = new LIST();
    for(long a=0; a<len(v1); ++a)
      returnValue->addValue(Value::from(homog(v1[a], h->theRingElem)));
    return returnValue;
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(abs, 1) {
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(abs(RefTo<BigInt>(v)));
  case 2: return Value::from(abs(RefTo<BigRat>(v)));
  case 3: return Value::from(abs(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(floor, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
	case 2: return Value::from(floor(RefTo<BigRat>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ceil, 1) {
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
	case 2: return Value::from(ceil(RefTo<BigRat>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(round, 1) {
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
	case 2: return Value::from(round(RefTo<BigRat>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(fields, 1) {
	intrusive_ptr<RECORD> r = runtimeEnv->evalArgAs<RECORD>(ARG(0));
	return r->fieldNames();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ConcatList, 1) {
	intrusive_ptr<LIST> listOfList = runtimeEnv->evalArgAs<LIST>(ARG(0));
	intrusive_ptr<LIST> result(new LIST);
	const LIST::ContainerType::size_type nLists = listOfList->size();
	for(LIST::ContainerType::size_type a=0; a<nLists; ++a) {
		intrusive_ptr<LIST> l = dynamic_pointer_cast<LIST>(listOfList->getValue(a));
		if (!l)
			throw RuntimeException("The argument is not a list of lists", ARG(0).exp);
		const LIST::ContainerType::size_type len = l->size();
		for(LIST::ContainerType::size_type b=0; b<len; ++b)
			result->addValue(l->getValue(b));
	}
	return result;
}
END_STD_BUILTIN_FUNCTION


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(gcd) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(gcd) {
  invocationExpression->checkNumberOfArgs(1,2);
  if (invocationExpression->args.size()==1)
  {
    int which;
    intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<LIST, INT, RINGELEM>(ARG(0), which);
    switch (which) {
    case 1: { // LIST
      intrusive_ptr<LIST> l = runtimeEnv->evalArgAs<LIST>(ARG(0));
      if (l->size() == 0) return INT::zero;
      vector<BigInt> v1;
      if (IsVectorBigInt(v1, l))  return Value::from(gcd_forC5(v1));
      return Value::from(gcd_forC5(runtimeEnv->evalArgAsListOfRingElem(ARG(0))));
    }
    case 2: return v;
    case 3: return v;
    default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
    }
  }
  int which0;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  int which1;
  intrusive_ptr<RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  if (which0==1 && which1==1)
    return Value::from(gcd(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
  if (which0==1)
  {
    RingElem x1 = RefTo<RingElem>(v1);
    return Value::from(gcd(RingElem(owner(x1), RefTo<BigInt>(v0)), x1));
  }
  if (which1==1)
  {
    RingElem x0 = RefTo<RingElem>(v0);
    return Value::from(gcd(x0, RingElem(owner(x0), RefTo<BigInt>(v1))));
  }
  return Value::from(gcd(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
}


DECLARE_ARITYCHECK_FUNCTION(lcm) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(lcm)
{
  invocationExpression->checkNumberOfArgs(1,2);
  if (invocationExpression->args.size()==1)
  {
    int which;
    intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<LIST, INT, RINGELEM>(ARG(0), which);
    switch (which) {
    case 1: { // LIST
      intrusive_ptr<LIST> l = runtimeEnv->evalArgAs<LIST>(ARG(0));
      if (l->size() == 0) return INT::one;
      vector<BigInt> v1;
      if (IsVectorBigInt(v1, l))  return Value::from(lcm_forC5(v1));
      return Value::from(lcm_forC5(runtimeEnv->evalArgAsListOfRingElem(ARG(0))));
    }
    case 2:  return v;
    case 3:  return v;
    default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
    }
  }
  int which0, which1;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  intrusive_ptr<RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  if (which0==1 && which1==1)
    return Value::from(lcm(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
  if (which0==1)
  {
    RingElem x1 = RefTo<RingElem>(v1);
    return Value::from(lcm(RingElem(owner(x1), RefTo<BigInt>(v0)), x1));
  }
  if (which1==1)
  {
    RingElem x0 = RefTo<RingElem>(v0);
    return Value::from(lcm(x0, RingElem(owner(x0), RefTo<BigInt>(v1))));
  }
  return Value::from(lcm(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
}


DECLARE_STD_BUILTIN_FUNCTION(num, 1) { // AMB (changed)
	const Argument &arg = ARG(0);
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(arg, which);
  switch (which) {
  case 1: return v;
  case 2: return Value::from(num(RefTo<BigRat>(v)));
  case 3: return Value::from(num(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(den, 1) {
	const Argument &arg = ARG(0);
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(arg, which);
  switch (which) {
  case 1: return Value::from(BigInt(1));
  case 2: return Value::from(den(RefTo<BigRat>(v)));
  case 3: return Value::from(den(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SimplestRatBetween, 2) { //JAA 2012-12-11
	const Argument &arg0 = ARG(0);
	const Argument &arg1 = ARG(1);
        BigRat A, B;
        int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg0, which);
  switch (which) {
  case 1: A = RefTo<BigInt>(v); break;
  case 2: A = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
	v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg1, which);
  switch (which) {
  case 1: B = RefTo<BigInt>(v); break;
  case 2: B = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
  return new RAT(SimplestBigRatBetween(A, B));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(sprint, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<OutputStringStreamValue> outs(new OutputStringStreamValue);
	outs->print(runtimeEnv, runtimeEnv->evalArgAs<RightValue>(arg), arg.exp);
	return outs->close(runtimeEnv, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GetErrMesg, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<ERROR> err = runtimeEnv->evalArgAs<ERROR>(arg);
	return new STRING(err->message);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(close, 1) {
	const Argument &arg = ARG(0);
	return runtimeEnv->evalArgAs<OSTREAM>(arg)->close(runtimeEnv, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenOString, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	return new OutputStringStreamValue();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GetEnv, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> varName = runtimeEnv->evalArgAs<STRING>(arg);
	char *value = ::getenv(varName->theString.c_str());
	if (!value)
		return STRING::empty;
	return new STRING(string(value));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenOFile, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> filename = runtimeEnv->evalArgAs<STRING>(arg);
	return new OutputFileStreamValue(filename, arg.exp);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ContentsOfFile, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> filename = runtimeEnv->evalArgAs<STRING>(arg);
	filtering_istream input;
	input.push(newline_filter(newline::posix));
	const file_source fs(filename->theString.c_str());
	if (!fs.is_open())
		throw RuntimeException("Cannot open file \""+filename->theString+"\" for reading.", invocationExpression);
	input.push(fs);
	string wholeFile;
	for(;;) {
		string line;
		getline(input, line);
		if (input.bad())
			throw RuntimeException("Cannot read from \""+filename->theString+"\".", invocationExpression);
		wholeFile += line;
		if (input.eof())
			break;
		wholeFile += '\n';
	}
	return new STRING(wholeFile);
}
END_STD_BUILTIN_FUNCTION

const string invalidTagError("The tag must be a valid identifier");

DECLARE_STD_BUILTIN_FUNCTION(tagged, 2) {
	const Argument &argValue = ARG(0);
	const Argument &argTag = ARG(1);
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAs<RightValue>(argValue);
	const string tag = runtimeEnv->evalArgAs<STRING>(argTag)->theString;
	// TODO: must be fixed
	//if (!LexerNS::Lexer::isValidIdentifier(tag))
	//	throw RuntimeException(invalidTagError, argTag.exp);
	return new TaggedValue(v->untagged(), invocationExpression->packageName+"."+tag);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(untagged, 1) {
	return runtimeEnv->evalArgAs<RightValue>(ARG(0))->untagged();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(tag, 1) {
	const Argument &argValue = ARG(0);
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAs<RightValue>(argValue);
	if (intrusive_ptr<TaggedValue> taggedValue = dynamic_pointer_cast<TaggedValue>(v)) {
		const string::size_type l = invocationExpression->packageName.length()+1;
		if (taggedValue->tag.substr(0, l)==invocationExpression->packageName+".")
			return new STRING(taggedValue->tag.substr(l, string::npos));
		return new STRING(taggedValue->tag);
	}
	return STRING::empty;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(type, 1) {
	return runtimeEnv->evalArgAs<RightValue>(ARG(0))->getType();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PkgName, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	return new STRING(invocationExpression->packageName);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(TAGGED, 1) {
	const Argument &arg = ARG(0);
	const string tag = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	if (!LexerNS::Lexer::isValidIdentifier(tag))
		throw RuntimeException(invalidTagError, arg.exp);
	return TYPE::tagType(tag);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(types, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	throw RuntimeException("Cannot return a list of all types because they're infinite (CoCoA 4 lies about this, I don't). You might want to use CurrentTypes() to get a list of currently instantiated types", invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CurrentTypes, 0) {
	return runtimeEnv->currentTypes();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(error, 1) {
	const Argument &arg = ARG(0);
	const string message = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	throw RuntimeException(message, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(Error, 1) {
	const Argument &arg = ARG(0);
	const string message = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	throw RuntimeException(message, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ascii, 1) {
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<STRING, INT, LIST>(ARG(0), which);
  switch (which) {
  case 1: {
    intrusive_ptr<LIST> l(new LIST());
		const string s(RefTo<string>(v));
		BOOST_FOREACH(const char c, s)
      l->addValue(new INT(static_cast<int>(c)));
		return l; }
	case 2: {
		long n;
		if (!IsConvertible(n, RefTo<BigInt>(v)) || n<0 || n>255)
			throw RuntimeException("When the argument is an INT, it must be in the range [0..255]", ARG(0).exp);
		return Value::from(string(1, static_cast<char>(n))); }
	case 3: {
		string retValue;
    intrusive_ptr<LIST> list = dynamic_pointer_cast<LIST>(v);
		LIST::ContainerType::size_type size = list->size();
		for(LIST::ContainerType::size_type a=0; a<size; ++a) {
			intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(list->getValue(a));
			long n;
			if (!N || !IsConvertible(n, N->theBigInt) || n<0 || n>255)
				throw RuntimeException("When the argument is a LIST, it must be a list of INT only in the range [0..255]", ARG(0).exp);
			retValue += static_cast<char>(n);
		}
		return Value::from(retValue);
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ILogBase, 2) { // AMB
  intrusive_ptr<INT> base = runtimeEnv->evalArgAs<INT>(ARG(1));
  int which;
	intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(ILogBase(RefTo<BigInt>(v0), base->theBigInt));
  case 2: return Value::from(ILogBase(RefTo<BigRat>(v0), base->theBigInt));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(log, 1) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  vector<BigInt> expv;
  if (NumIndets(owner(f->theRingElem))==1)
    expv.push_back(BigInt(deg(f->theRingElem)));
  else 
    BigExponents(expv, LPP(f->theRingElem));
	return new LIST(expv);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(wdeg, 1) { // AMB
  int which;
  intrusive_ptr<const RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(wdeg(RefTo<RingElem>(x)));
  case 2: return Value::from(wdeg(RefTo<ModuleElem>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsHomog, 1) { // AMB
	int which;
	intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM, MODULEELEM, IDEAL, MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsHomog(RefTo<RingElem>(v)));
  case 2: return Value::from(IsHomog(RefTo<ModuleElem>(v)));
  case 3: return Value::from(IsHomog(RefTo<ideal>(v)));
  case 4: return Value::from(IsHomog(RefTo<module>(v)));
  case 5: {
    if (const intrusive_ptr<const LIST> l = intrusive_ptr_cast<const LIST>(v))
    {
      LIST::ContainerType::size_type size = l->size();
      for(LIST::ContainerType::size_type a=0; a<size; ++a)
      {
        if (const boost::intrusive_ptr<const RINGELEM> f = boost::dynamic_pointer_cast<const RINGELEM>(l->getValue(a)))
        {  if (!IsHomog(f->theRingElem)) return Value::from(false);}
        else throw RuntimeException("RingElem expected", ARG(1).exp);
      }
      return Value::from(true);
    }
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsFactorClosed, 1)
{ // JAA
  const vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  if (v.empty()) CoCoA_ERROR(ERR::Empty, "IsFactorClosed");
  const long n = len(v);
  if (n == 1 && IsOne(v[0])) return Value::from(true);
  const ring R = owner(v[0]);
  for (long i=1; i < n; ++i)
    if (owner(v[i]) != R) CoCoA_ERROR(ERR::MixedRings, "IsFactorClosed");
  if (!IsSparsePolyRing(R)) CoCoA_ERROR(ERR::NotSparsePolyRing, "IsFactorClosed");
  for (long i=0; i < n; ++i)
    if (!IsMonomial(v[i])) CoCoA_ERROR(ERR::NotMonomialGens, "IsFactorClosed");
  vector<PPMonoidElem> SetPP; SetPP.reserve(n);
  for (long i=0; i < n; ++i)
    SetPP.push_back(LPP(v[i]));
  return Value::from(IsFactorClosed(SetPP));
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(LPP, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LPP_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LPP_forC5(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LC, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LC(RefTo<RingElem>(v)));
  case 2: return Value::from(LC(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LM, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LM_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LM_forC5(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LT, 1) { // AMB
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4<RINGELEM, MODULEELEM, IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LT_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LT_forC5(RefTo<ModuleElem>(v)));
  case 3: return Value::from(LT(RefTo<ideal>(v)));
  case 4: return Value::from(LT(RefTo<module>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LF, 1) { // AMB
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4<RINGELEM, MODULEELEM, IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LF(RefTo<RingElem>(v)));
  case 2: CoCoA_ERROR(ERR::NYI,"LF for MODULEELEM");
    return Value::from(false); // just to keep the compiler quiet
    //return Value::from(LF(RefTo<ModuleElem>(v)));
  case 3: return Value::from(LF(RefTo<ideal>(v)));
  case 4: CoCoA_ERROR(ERR::NYI,"LF for MODULE");
    //return Value::from(LF(RefTo<module>(v)));
    return Value::from(false); // just to keep the compiler quiet
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RingOf, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM,IDEAL,MAT,MODULE,MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(owner(RefTo<RingElem>(v)));
  case 2: return Value::from(RingOf(RefTo<ideal>(v)));
  case 3: return Value::from(RingOf(RefTo<matrix>(v)));
  case 4: return Value::from(RingOf(RefTo<module>(v)));
  case 5: return Value::from(RingOf(owner(RefTo<ModuleElem>(v))));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ModuleOf, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<MODULEELEM, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(owner(RefTo<ModuleElem>(v)));
  case 2: return Value::from(AmbientFreeModule(RefTo<module>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


//---- RINGHOM -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(InducedHom, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(1));
  if (IsFractionField(R->theRing))
    return Value::from(InducedHom(FractionField(R->theRing),phi->theRingHom));
  if (IsQuotientRing(R->theRing))
    return Value::from(InducedHom(QuotientRing(R->theRing), phi->theRingHom));
  throw RuntimeException("FractionField or QuotientRing expected", ARG(0).exp);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(PolyAlgebraHom, 3) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(1));
  vector<RingElem> IndetImages = runtimeEnv->evalArgAsListOfRingElem(ARG(2), R->theRing);
	return Value::from(PolyAlgebraHom(P->theRing, R->theRing, IndetImages));
}
END_STD_BUILTIN_FUNCTION 


DECLARE_STD_BUILTIN_FUNCTION(PolyRingHom, 4) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(1));
  intrusive_ptr<RINGHOM> CoeffHom = runtimeEnv->evalArgAs<RINGHOM>(ARG(2));
  vector<RingElem> IndetImages = runtimeEnv->evalArgAsListOfRingElem(ARG(3), R->theRing);
	return new RINGHOM(PolyRingHom(P->theRing, R->theRing, CoeffHom->theRingHom, IndetImages));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(apply, 2) {
	intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(0));
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<RINGELEM, MAT, LIST>(ARG(1), which);
  switch (which) {
  case 1:	return Value::from(apply(phi->theRingHom, RefTo<RingElem>(v)));
  case 2: return Value::from(apply(phi->theRingHom, RefTo<matrix>(v)));
  case 3: {
    vector<RingElem> w = runtimeEnv->evalArgAsListOfRingElem(ARG(1), domain(phi->theRingHom));
		return Value::from(apply(phi->theRingHom, w));
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(indets) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(indets) { // AMB+JAA
  invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(indets((runtimeEnv->evalArgAs<RING>(ARG(0)))->theRing));
  return Value::from(indets((runtimeEnv->evalArgAs<RING>(ARG(0)))->theRing,
                            runtimeEnv->evalArgAs<STRING>(ARG(1))->theString));
}


DECLARE_STD_BUILTIN_FUNCTION(indet, 2) { // AMB
	intrusive_ptr<RING> a = runtimeEnv->evalArgAs<RING>(ARG(0));
	intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n;
  if (!IsConvertible(n, b->theBigInt))
    throw RuntimeException("invalid indet index", ARG(1).exp);
	return Value::from(indet(a->theRing, n-1));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetIndex, 1) { // AMB
  long i;
  if (!IsIndet(i, runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  return Value::from(i+1);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetName, 1) { // AMB
	intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  long i;
  if (!IsIndet(i,a->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  const symbol s=IndetSymbol(owner(a->theRingElem), i);
	return Value::from(head(s));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetSubscripts, 1) { // AMB
	intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  long i;
  if (!IsIndet(i,a->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  const symbol s=IndetSymbol(owner(a->theRingElem), i);
  vector<long> iss;
  for (long n=0; n<NumSubscripts(s); ++n) iss.push_back(subscript(s,n));
	return Value::from(iss);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ContentWRT, 2) { // AMB
	intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  int which;
	intrusive_ptr<RightValue> idt = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(ContentWRT(f->theRingElem, RefTo<RingElem>(idt)));
  case 2: {
//     vector<long> v1 = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadIndetIndex, "ContentWRT");
//     return Value::from(ContentWRT(f->theRingElem, v1));
    vector<RingElem> v1 = runtimeEnv->evalArgAsListOfRingElem(ARG(1), owner(f->theRingElem));
    return Value::from(ContentWRT_forC5(f->theRingElem, v1));
  }
  default: throw RuntimeException(ERRORMissingCode(idt),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CoefficientsWRT, 2) { // AMB
	intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  int which;
	intrusive_ptr<RightValue> idt = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(1), which);
  std::vector<CoeffPP> CoeffWRT;
  switch (which) {
  case 1: CoeffWRT=CoefficientsWRT(f->theRingElem, RefTo<RingElem>(idt)); break;
  case 2: {
//     vector<long> v1 = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadIndetIndex, "ContentWRT");
//     CoeffWRT = CoefficientsWRT(f->theRingElem, v1);
    vector<RingElem> v1 = runtimeEnv->evalArgAsListOfRingElem(ARG(1), owner(f->theRingElem));
    CoeffWRT = CoefficientsWRT_forC5(f->theRingElem, v1);
    break;
  }
  default: throw RuntimeException(ERRORMissingCode(idt),invocationExpression);
  }
  // create return value: list of ...
  const SparsePolyRing P = owner(f->theRingElem);
	intrusive_ptr<LIST> returnValue(new LIST);
  BOOST_FOREACH(const CoeffPP cpp, CoeffWRT)
  {
    // create return value: ... record
    intrusive_ptr<RECORD> cpp5(new RECORD);
    cpp5->setField("coeff", Value::from(cpp.myCoeff));
    cpp5->setField("PP",    Value::from(monomial(P,1,cpp.myPP)));
    returnValue->addValue(cpp5);
	}
  return returnValue;
}
END_STD_BUILTIN_FUNCTION



//---- MAT -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(SetEntry, 4) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
	intrusive_ptr<INT> J = runtimeEnv->evalArgAs<INT>(ARG(2));
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(3), which);
  RingElem x(RingOf(M->theMatrix));
  switch (which) {
  case 1: x = RefTo<BigInt>(v); break;
  case 2: x = RefTo<BigRat>(v); break;
  case 3: x = RefTo<RingElem>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v), invocationExpression);
  }
	SetEntry_forC5(M->theMatrix, I->theBigInt, J->theBigInt, x);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetRow, 3) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(2), RingOf(M->theMatrix));
	SetRow_forC5(M->theMatrix, I->theBigInt, v);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetCol, 3) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(2), RingOf(M->theMatrix));
	SetCol_forC5(M->theMatrix, I->theBigInt, v);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SwapRows, 3) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	intrusive_ptr<INT> I1 = runtimeEnv->evalArgAs<INT>(ARG(1));
	intrusive_ptr<INT> I2 = runtimeEnv->evalArgAs<INT>(ARG(2));
	SwapRows_forC5(M->theMatrix, I1->theBigInt, I2->theBigInt);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SwapCols, 3) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	intrusive_ptr<INT> I1 = runtimeEnv->evalArgAs<INT>(ARG(1));
	intrusive_ptr<INT> I2 = runtimeEnv->evalArgAs<INT>(ARG(2));
	SwapCols_forC5(M->theMatrix, I1->theBigInt, I2->theBigInt);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(AssignZero, 1) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
	AssignZero(M->theMatrix);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsZeroRow, 2) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  const ErrorInfo ErrMesg(ERR::BadRowIndex, "IsZeroRow");
  return Value::from(IsZeroRow(M->theMatrix, ConvertTo<long>(N->theBigInt, ErrMesg)-1));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsZeroCol, 2) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  const ErrorInfo ErrMesg(ERR::BadColIndex, "IsZeroCol");
  return Value::from(IsZeroCol(M->theMatrix, ConvertTo<long>(N->theBigInt, ErrMesg)-1));
}
END_STD_BUILTIN_FUNCTION

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(IsPositiveGrading) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(IsPositiveGrading) { // AMB variable number of args
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(IsPositiveGrading(M->theMatrix)); // 1 arg: done
  intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  long GrDim = NumRows(M->theMatrix);
  if (!IsConvertible(GrDim, N->theBigInt))
    throw RuntimeException("invalid row index", ARG(1).exp);
  return Value::from(IsPositiveGrading(M->theMatrix, GrDim));
}
//END_STD_BUILTIN_FUNCTION // no: variable number of args

DECLARE_STD_BUILTIN_FUNCTION(ElimMat, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, MAT>(ARG(0), which);
  vector<BigInt> elim = runtimeEnv->evalArgAsListOf<INT>(ARG(1));
  switch (which) {
  case 1: return Value::from(ElimMat_forC5(RefTo<BigInt>(x), elim));
  case 2: return Value::from(ElimMat_forC5(RefTo<matrix>(x), elim));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(submat, 3) { // AMB
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  vector<long> rows = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadRowIndex, "submat");
  vector<long> cols = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(2)), ERR::BadColIndex, "submat");
	return new MAT(NewDenseMat(submat(M->theMatrix, rows, cols)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NR, 2) { // AMB
	intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
	vector<RingElem> v = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(1));
	return Value::from(NR(f->theRingElem,v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(BinSequentialToric, 2) { // AMB
  intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
  vector<BigInt> indices = runtimeEnv->evalArgAsListOf<INT>(ARG(1));
	return Value::from(SequentialToric_C(I->theIdeal, VectorLong(indices, "BinSequentialToric")));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(MatSequentialToric, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(1));
	return Value::from(SequentialToric_C(R->theRing, M->theMatrix));
}
END_STD_BUILTIN_FUNCTION

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(ideal) { return 1<=nArg; }
DECLARE_BUILTIN_FUNCTION(ideal) {
	const int nArgs = invocationExpression->args.size();
	if (nArgs==0)
		throw RuntimeException("Wrong number of arguments; found: 0, expecting: at least 1", invocationExpression);
  int which;
	intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<RING, LIST, RINGELEM>(ARG(0), which);
  if (nArgs==2 && which==1)
    return Value::from(ideal(RefTo<ring>(x), runtimeEnv->evalArgAsListOfRingElem(ARG(1), RefTo<ring>(x))));
  if (nArgs==1 && which==2)
    return Value::from(ideal(runtimeEnv->evalArgAsListOfRingElem(ARG(0))));
	return Value::from(ideal(runtimeEnv->evalAllArgsAsListOf<RINGELEM>(invocationExpression)));
}

//---- MODULE -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(ModuleElem, 2) { // AMB
  intrusive_ptr<MODULE> M = runtimeEnv->evalArgAs<MODULE>(ARG(0));
	vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(1), RingOf(M->theModule));
  return Value::from(NewFreeModuleElem(M->theModule, v));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NewFreeModule, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
	int which;
	intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, MAT>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(NewFreeModule(R->theRing, ConvertTo<long>(RefTo<BigInt>(x))));
  case 2: return Value::from(NewFreeModule_forC5(R->theRing, RefTo<matrix>(x)));
  default: throw RuntimeException(ERRORMissingCode(x), invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(shifts, 1) { // AMB
	intrusive_ptr<MODULE> v = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(shifts(v->theModule));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewFreeModuleForSyz, 1) { // AMB
	int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<IDEAL, MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(NewFreeModuleForSyz(gens(RefTo<ideal>(v))));
  case 2: return Value::from(NewFreeModuleForSyz(gens(RefTo<module>(v))));
  case 3: {
    vector<RingElem> x = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
    return Value::from(NewFreeModuleForSyz(x));
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(submodule) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(submodule) {  // AMB
	invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(submodule(RefTo<module>(v0), runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(1))));
  case 2: return Value::from(submodule(runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(0))));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}


DECLARE_STD_BUILTIN_FUNCTION(GensAsRows, 1) { // AMB
  intrusive_ptr<MODULE> F = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(GensAsRows(F->theModule));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GensAsCols, 1) { // AMB
  intrusive_ptr<MODULE> F = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(GensAsCols(F->theModule));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NumCompts, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<MODULE, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(NumCompts_forC5(RefTo<module>(x)));
  case 2: return Value::from(NumCompts(RefTo<ModuleElem>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


//---- IDEAL -------------------------------------------------------

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(janet) { return 1<=nArg; }
DECLARE_BUILTIN_FUNCTION(janet) {
	const int nArgs = invocationExpression->args.size();
	if (nArgs==0)
		throw RuntimeException("Wrong number of arguments; found: 0, expecting: at least 1", invocationExpression);
	intrusive_ptr<RightValue> x = runtimeEnv->evalArgAs<RightValue>(ARG(0));
  if (invocationExpression->args.size()==2)
    if (intrusive_ptr<RING> R = dynamic_pointer_cast<RING>(x))
      return Value::from(ExtendedJanetBasis(runtimeEnv->evalArgAsListOfRingElem(ARG(1), R->theRing)));
	return Value::from(ExtendedJanetBasis(runtimeEnv->evalAllArgsAsListOf<RINGELEM>(invocationExpression)));
}

//---- RING -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(NewRingTwinFloat, 1) { // AMB
	intrusive_ptr<INT> prec = runtimeEnv->evalArgAs<INT>(ARG(0));
  long d;
  if (!IsConvertible(d, prec->theBigInt))
		throw RuntimeException("invalid precision", ARG(0).exp);
	return Value::from(NewRingTwinFloat(d));
}
END_STD_BUILTIN_FUNCTION


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(NewPolyRing) { return (nArg==2) || (nArg==4); }
DECLARE_BUILTIN_FUNCTION(NewPolyRing) {
  //	invocationExpression->checkNumberOfArgs(2,4);
	const int nArgs = invocationExpression->args.size();
	if (nArgs!=2 && nArgs!=4)
		throw RuntimeException("Wrong number of arguments; found: "+boost::lexical_cast<std::string>(nArgs)+", expecting: 2 or 4", invocationExpression);
	intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  //	vector<string> v = runtimeEnv->evalArgAsListOf<STRING>(ARG(1));
	intrusive_ptr<LIST> l = runtimeEnv->evalArgAs<LIST>(ARG(1));
	LIST::ContainerType::size_type size = l->size();
  vector<symbol> v;
  for (unsigned long i=0; i<size; ++i)
  {
		if (const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(l->getValue(i)))
      v.push_back(symbol(s->theString));
    else if (const boost::intrusive_ptr<RECORD> r = boost::dynamic_pointer_cast<RECORD>(l->getValue(i)))
    {
      const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(r->getField("head"));
      const boost::intrusive_ptr<LIST> inds = boost::dynamic_pointer_cast<LIST>(r->getField("indices"));
      vector<long> indices;
      long tmp;
      LIST::ContainerType::size_type NumIndices = inds->size();
      for (unsigned long j=0; j<NumIndices; ++j)
        if (IsConvertible(tmp, (boost::dynamic_pointer_cast<INT>(inds->getValue(j)))->theBigInt))
          indices.push_back(tmp);
        else
          throw RuntimeException("All indices must fit into a machine-integer", ARG(1).exp);
      v.push_back(symbol(s->theString, indices));
    }
    else throw RuntimeException("List of symbols must have strings or records", ARG(1).exp);
  }
  if (nArgs==2)
    return Value::from(NewPolyRing(R->theRing, v));
	intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(2));
	intrusive_ptr<INT> GradingDim = runtimeEnv->evalArgAs<INT>(ARG(3));
  long d;
  if (!IsConvertible(d, GradingDim->theBigInt))
		throw RuntimeException("invalid GradingDim", ARG(3).exp);
  const PPOrdering PPO = NewMatrixOrdering(v.size(), d, M->theMatrix);
	return new RING(NewPolyRing(R->theRing, NewPPMonoid(v, PPO)));
}
  //END_STD_BUILTIN_FUNCTION // no: variable number of args

//---- IDEAL -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(elim, 2) { // AMB
	intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
  ideal J = I->theIdeal;
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<LIST, RINGELEM>(ARG(0), which);
  vector<RingElem> ElimIndets;
  switch (which) {
  case 1: ElimIndets = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(0));break;
  case 2: ElimIndets.push_back(RefTo<RingElem>(x));break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
  MakeUnique(J)->myElim(ElimIndets); // I->myElim(v1) ~~~> const
	return Value::from(J);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(minimalized, 1) { // AMB
  int which;
	intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<MODULE, IDEAL>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(minimalized(RefTo<module>(x)));
  case 2: {
    ideal J = RefTo<ideal>(x);
    MinGens(J); // so GBasis and such are stored in original ideal
    MakeUnique(J)->myMinimalize();
    return Value::from(J);
  }
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(saturate, 2) { // AMB
	intrusive_ptr<IDEAL> I0 = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
	intrusive_ptr<IDEAL> I1 = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
  ideal J = I0->theIdeal;
  MakeUnique(J)->mySaturate(I1->theIdeal);
	return Value::from(J);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IdealOfPoints, 2) { // JAA 2013-01-19
	intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
	intrusive_ptr<MAT> pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
  return Value::from(IdealOfPoints(P->theRing, pts->theMatrix));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(TmpNBM, 3) { // AMB
	intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
	intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
	intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
  vector<RingElem> QB;
  vector<RingElem> BB;
  vector<RingElem> AV;
  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("QuotientBasis", Value::from(QB));
  rec->setField("AlmostVanishing", Value::from(AV));
  if (BB.empty())
    rec->setField("StableBBasisFound", Value::from(false));
  else
  {
    rec->setField("StableBBasisFound", Value::from(true));
    rec->setField("BBasis", Value::from(BB));
  }
  return rec;
}
END_STD_BUILTIN_FUNCTION

// DECLARE_STD_BUILTIN_FUNCTION(TmpSOI, 3) { // AMB
// 	intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
// 	intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
// 	intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
//   vector<RingElem> QB;
//   vector<RingElem> BB;
//   vector<RingElem> AV;
//   SOI_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
//   intrusive_ptr<RECORD> rec(new RECORD);
//   rec->setField("QuotientBasis", Value::from(QB));
//   rec->setField("AlmostVanishing", Value::from(AV));
//   if (BB.empty())
//     rec->setField("StableBBasisFound", Value::from(false));
//   else
//   {
//     rec->setField("StableBBasisFound", Value::from(true));
//     rec->setField("BBasis", Value::from(BB));
//   }
//   return rec;
// }
// END_STD_BUILTIN_FUNCTION


// DECLARE_STD_BUILTIN_FUNCTION(ClosePassingPoly, 3) { // AMB
// 	intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
// 	intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
// 	intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
// 	intrusive_ptr<RAT> MaxTol = runtimeEnv->evalArgAs<RAT>(ARG(3));
//   ClosePassingPoly_forC5(P->theRing, Pts->theMatrix, Tols->theMatrix);
//   return Value::from(ClosePassingPoly_forC5(P->theRing, Pts->theMatrix, Tols->theMatrix);
// }
// END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(PreprocessPts, 2) { // JAA
	intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("auto", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setField("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsGrid, 2) { // JAA
	intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("grid", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setField("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsAggr, 2) { // JAA
	intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("aggr", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setField("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsSubdiv, 2) { // JAA
	intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
	intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("subdiv", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setField("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION


//---- general ----

DECLARE_STD_BUILTIN_FUNCTION(TopLevelFunctions, 0) { // GL
	return runtimeEnv->topLevelFunctions();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CocoaPackagePath, 0) { // GL
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
	return Value::from(packageDir);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(VersionInfo, 0) { // AMB
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setField("CoCoALibVersion", Value::from(BuildInfo::version()));
  rec->setField("CoCoAVersion", Value::from(CoCoAVersion()));
  rec->setField("CompilationDate", Value::from(CompilationDate()));
  rec->setField("CompilationFlags", Value::from(BuildInfo::CompilationFlags()));
  rec->setField("CompilationDefines", Value::from(BuildInfo::CompilationDefines()));
  rec->setField("Compiler", Value::from(BuildInfo::compiler()));
  rec->setField("MachineIntNumBits", new INT(std::numeric_limits<unsigned int>::digits));
  rec->setField("MachineLongNumBits", new INT(std::numeric_limits<unsigned long>::digits));

  return rec;
}
END_STD_BUILTIN_FUNCTION

//---- manual CoCoAHelp ----

DECLARE_STD_BUILTIN_FUNCTION(ReloadMan, 0) { // AMB
	try {
    ostringstream os;
    OnlineHelp::ReloadMan(os);
    runtimeEnv->getOutputStream()->print(os.str())->flush();
	} catch (const std::exception& err) {
		throw RuntimeException(err.what(), invocationExpression->getBegin(), invocationExpression->getEnd());
	}
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


  //----------------------------------------------------------------------
void RuntimeEnvironment::initBuiltInFunctions()
{
  BOOST_FOREACH(NameFunPair &p, builtIns())
		this->setTopLevelVar(p.first, p.second, VariableSlot::VSF_SystemProtected);
}
  //----------------------------------------------------------------------

} // namespace InterpreterNS
} // namespace CoCoA
