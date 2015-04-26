//   Copyright (c) 2010-2012 Giovanni Lagorio and John Abbott, Anna M. Bigatti
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

#ifndef BUILTINFUNCTIONS_H_
#define BUILTINFUNCTIONS_H_

#include <iostream>
#include <sstream>
#include <exception>
#include <cassert>
#include <vector>
#include <cstdlib>
#include <memory>
#include <boost/scope_exit.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>
#if 0
#include <zip.h>
#endif

#include "CoCoA/library.H"
#include "CoCoALibSupplement.H"
#include "Interpreter.H"
#include "CompilationDate.H"
#include "OnlineHelp.H"
#include "VersionInfo.H"

// using namespace std;
// using namespace boost;
// using namespace boost::iostreams;
// using namespace CoCoA::AST;
// using namespace CoCoA::LexerNS;
// using namespace CoCoA::ParserNS;

namespace CoCoA {
namespace InterpreterNS {

template <> BigInt inline theValue<INT>(boost::intrusive_ptr<INT> x) { return x->theBigInt; }
template <> const BigInt inline theValue<const INT>(boost::intrusive_ptr<const INT> x) { return x->theBigInt; }
template <> BigRat inline theValue<RAT>(boost::intrusive_ptr<RAT> x) { return x->theBigRat; }
template <> const BigRat inline theValue<const RAT>(boost::intrusive_ptr<const RAT> x) { return x->theBigRat; }
template <> RingElem inline theValue<RINGELEM>(boost::intrusive_ptr<RINGELEM> x) { return x->theRingElem; }
template <> const RingElem inline theValue<const RINGELEM>(boost::intrusive_ptr<const RINGELEM> x) { return x->theRingElem; }
template <> ring inline theValue<RING>(boost::intrusive_ptr<RING> x) { return x->theRing; }
template <> const ring inline theValue<const RING>(boost::intrusive_ptr<const RING> x) { return x->theRing; }
template <> RingHom inline theValue<RINGHOM>(boost::intrusive_ptr<RINGHOM> x) { return x->theRingHom; }
template <> const RingHom inline theValue<const RINGHOM>(boost::intrusive_ptr<const RINGHOM> x) { return x->theRingHom; }
template <> matrix inline theValue<MAT>(boost::intrusive_ptr<MAT> x) { return x->theMatrix; }
template <> const matrix inline theValue<const MAT>(boost::intrusive_ptr<const MAT> x) { return x->theMatrix; }
template <> std::string inline theValue<STRING>(boost::intrusive_ptr<STRING> x) { return x->theString; }
template <> const std::string inline theValue<const STRING>(boost::intrusive_ptr<const STRING> x) { return x->theString; }
template <> std::string inline theValue<TYPE>(boost::intrusive_ptr<TYPE> x) { return x->name; }
template <> const std::string inline theValue<const TYPE>(boost::intrusive_ptr<const TYPE> x) { return x->name; }
template <> bool inline theValue<BOOL>(boost::intrusive_ptr<BOOL> x) { return x->theBool; }
template <> bool inline theValue<const BOOL>(boost::intrusive_ptr<const BOOL> x) { return x->theBool; }
template <> ideal inline theValue<IDEAL>(boost::intrusive_ptr<IDEAL> x) { return x->theIdeal; }
template <> const ideal inline theValue<const IDEAL>(boost::intrusive_ptr<const IDEAL> x) { return x->theIdeal; }
template <> module inline theValue<MODULE>(boost::intrusive_ptr<MODULE> x) { return x->theModule; }
template <> const module inline theValue<const MODULE>(boost::intrusive_ptr<const MODULE> x) { return x->theModule; }
template <> ModuleElem inline theValue<MODULEELEM>(boost::intrusive_ptr<MODULEELEM> x) { return x->theModuleElem; }
template <> const ModuleElem inline theValue<const MODULEELEM>(boost::intrusive_ptr<const MODULEELEM> x) { return x->theModuleElem; }

  //---------------------------------------
template <typename T> const T& RefTo(boost::intrusive_ptr<const RightValue>);  

template <> const BigInt& RefTo<BigInt>(boost::intrusive_ptr<const RightValue>);
template <> const BigRat& RefTo<BigRat>(boost::intrusive_ptr<const RightValue>);
template <> const ModuleElem& RefTo<ModuleElem>(boost::intrusive_ptr<const RightValue>);
template <> const RingElem& RefTo<RingElem>(boost::intrusive_ptr<const RightValue>);
template <> const ideal& RefTo<ideal>(boost::intrusive_ptr<const RightValue>);
template <> const matrix& RefTo<matrix>(boost::intrusive_ptr<const RightValue>);
template <> const module& RefTo<module>(boost::intrusive_ptr<const RightValue>);
template <> const ring& RefTo<ring>(boost::intrusive_ptr<const RightValue>);
template <> const std::string& RefTo<std::string>(boost::intrusive_ptr<const RightValue>);
  //--
template <> inline const BigInt& RefTo<BigInt>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const INT>(x)->theBigInt; }
template <> inline const BigRat& RefTo<BigRat>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RAT>(x)->theBigRat; }
template <> inline const ModuleElem& RefTo<ModuleElem>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MODULEELEM>(x)->theModuleElem; }
template <> inline const RingElem& RefTo<RingElem>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RINGELEM>(x)->theRingElem; }
template <> inline const ideal& RefTo<ideal>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const IDEAL>(x)->theIdeal; }
template <> inline const matrix& RefTo<matrix>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MAT>(x)->theMatrix; }
template <> inline const module& RefTo<module>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MODULE>(x)->theModule; }
template <> inline const ring& RefTo<ring>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RING>(x)->theRing; }
template <> inline const std::string& RefTo<std::string>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const STRING>(x)->theString; }
  //-----------------------
inline std::string ERRORMissingCode(const boost::intrusive_ptr<const RightValue>& x)
{ return "!!!CODE MISSING!!! for "+x->getType()->name; }



  //------------------------------------------------------------
  // The fn  builtIns  behaves like a global variable;  it is filled
  // by the ctor for AddBuiltIn, and the filled vector is used in the
  // the fn initBuiltInFunctions -- see BuiltInFunctions.C

typedef std::pair<std::string, boost::intrusive_ptr<BuiltInFunction> > NameFunPair;

std::vector<NameFunPair>& builtIns(); // declared in BuiltInFunctions.C

struct AddBuiltIn
{
  inline AddBuiltIn(const std::string &fname, boost::intrusive_ptr<BuiltInFunction> fun)
  {
    builtIns().push_back(std::make_pair(fname, fun));
  }
};

  //------------------------------------------------------------

#define ARG(N) invocationExpression->args[N]

#define DECLARE_ARITYCHECK_FUNCTION(FNAME) bool builtin_check_arity_##FNAME(int nArg)

#define DECLARE_BUILTIN_FUNCTION(FNAME) \
		intrusive_ptr<RightValue> builtin_##FNAME(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment * const runtimeEnv); \
		DECLARE_ARITYCHECK_FUNCTION(FNAME) ;\
		AddBuiltIn add_##FNAME(#FNAME, new BuiltInFunction(builtin_##FNAME, builtin_check_arity_##FNAME)); \
		intrusive_ptr<RightValue> builtin_##FNAME(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment * const runtimeEnv)

#define DECLARE_STD_BUILTIN_FUNCTION(FNAME, NARG) \
	DECLARE_ARITYCHECK_FUNCTION(FNAME) { return nArg==(NARG); } \
	DECLARE_BUILTIN_FUNCTION(FNAME) { \
	invocationExpression->checkNumberOfArgs(NARG);

#define END_STD_BUILTIN_FUNCTION }

#define DECLARE_MISSING_EXTLIB(FnName, ExtLibName)\
DECLARE_ARITYCHECK_FUNCTION(FnName) { return true; }\
DECLARE_BUILTIN_FUNCTION(FnName) \
{\
  CoCoA_ERROR(#FnName " undefined: CoCoA was built without ExternalLib " #ExtLibName, "BuiltInFunction"); \
  return VoidValue::theInstance; \
}

} // namespace AST
} // namespace CoCoA

#endif /* BUILTINFUNCTIONS_H_ */
