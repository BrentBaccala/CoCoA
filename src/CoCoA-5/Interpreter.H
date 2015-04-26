//   Copyright (c) 2010 Giovanni Lagorio (lagorio@disi.unige.it)
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

#ifndef INTERPRETER_H_
#define INTERPRETER_H_

#include "CoCoA/matrix.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/factorization.H"
#include "CoCoA/ideal.H"
#include "CoCoA/symbol.H"
#include "CoCoA/TmpJBEnv.H"

#include <ostream>
#include <map>
#include <vector>
#include <set>
#include <cassert>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cstdint.hpp>
#ifdef C5IDE
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
class QTreeWidget;
class QTreeWidgetItem;
#endif // #ifdef C5IDE
#include <sys/types.h>
#include <sys/stat.h>
#include "Parser.H"

namespace CoCoA {

  class ErrorInfo; // forward declaration -- defined in error.H

#ifdef C5IDE
namespace IDE {
	class Debugger;
	class Console;
	int launchTheIDE(CoCoA::LexerNS::WarningSeverity warningLevel, bool warnAboutCocoa5, const std::vector<std::string> &packageList, bool fullCoCoALibError);
	bool loadPackages(Console *console);
} // namespace IDE
#endif // #ifdef C5IDE

namespace InterpreterNS {

class RuntimeEnvironment;
class Frame;

#ifdef ULONG_FRAMEIDENTIFIERS
typedef unsigned long FrameIdentifier;
#else // #ifdef ULONG_FRAMEIDENTIFIERS
typedef uint_fast64_t FrameIdentifier; // if you get a compilation error here, then compile defining the macro ULONG_FRAMEIDENTIFIERS
#endif // #ifdef ULONG_FRAMEIDENTIFIERS

class FramePointer {
	Frame *frame;
	FrameIdentifier id;
public:
	FramePointer(Frame *const frame=0);
	Frame *toCheckedPointer() const;
	Frame *toNonNullCheckedPointer(boost::intrusive_ptr<const AST::Identifier> expId) const;
	inline void reset() { this->frame=0; }
};

class Value;

class VariableSlot {
public:
	enum Flags {
		VSF_None,
		VSF_Protected=1,
		VSF_SystemProtected = VSF_Protected|2,
		VSF_Args = VSF_SystemProtected|4,
		VSF_Package = 8,
		VSF_IterationVariable = 16,
		VSF_Exported = 32
	};
	boost::intrusive_ptr<Value> value;
private:
	Flags flags;
	friend class RuntimeEnvironment;
	std::string protectionReason;
public:
	inline Flags getFlags() const { return this->flags; }
	inline void resetFlags() { this->flags = VSF_None; }
	inline bool isPackage() const { return (this->flags & VSF_Package)==VSF_Package; }
	inline bool isProtected() const { return (this->flags & VSF_Protected)==VSF_Protected; }
	inline bool isArgs() const { return (this->flags & VSF_Args)==VSF_Args; }
	inline bool isSystemProtected() const { return (this->flags & VSF_SystemProtected)==VSF_SystemProtected; }
	inline bool isIterationVariable() const { return (this->flags & VSF_IterationVariable)==VSF_IterationVariable; }
	inline bool hasBeenExported() const { return (this->flags & VSF_Exported)==VSF_Exported; }
	inline void setHasBeenExported() { this->flags = static_cast<Flags>(this->flags|VSF_Exported); }
	void protect(std::string reason) {
		this->setProtectionReason(reason);
		this->flags = static_cast<Flags>(this->flags|VSF_Protected);
	}
	void unprotect() {
		this->flags = static_cast<Flags>(this->flags&~VSF_SystemProtected);
	}
	const std::string& getProtectionReason() const {
		return this->protectionReason;
	}
	void setProtectionReason(const std::string &reason) {
		this->protectionReason = reason;
	}
	explicit VariableSlot(boost::intrusive_ptr<Value> value, Flags flags) :
			value(value),
			flags(flags)
	{}
};

class UserDefinedFunction;

class Frame {
public: // all members are public
	FrameIdentifier id;
	FramePointer accessLink;
	std::vector<VariableSlot> varSlots;
	boost::intrusive_ptr<const AST::InvocationExpression> invocationExp;
	boost::intrusive_ptr<const AST::ParsedObject> block;
	boost::intrusive_ptr<const UserDefinedFunction> userdefinedFun;
	Frame() : id(1) {}
#ifdef C5IDE
	bool singleStepOnPop;
#endif // #ifdef C5IDE
};

class DeadEnviromentException : public RuntimeException {
	static std::string errorMessage(const std::string &id) {
		return "The environment for \""+id+"\" does not exist anymore";
	}
public:
	DeadEnviromentException(boost::intrusive_ptr<const AST::Identifier> id) :
		RuntimeException(errorMessage(id->identifier), id)
	{}
	DeadEnviromentException(const LexerNS::Token &tokId) :
		RuntimeException(errorMessage(tokId.lexeme()), tokId)
	{}
	~DeadEnviromentException() throw () {}
};

class WrongNumberOfArgumentsException : public RuntimeException {
	static std::string formatExpectedArgs(int found, int min, int max);
public:
	WrongNumberOfArgumentsException(boost::intrusive_ptr<const AST::ParsedObject> po, int found, int minExpecting, int maxExpecting) :
		RuntimeException(formatExpectedArgs(found, minExpecting, maxExpecting), po)
	{}
	WrongNumberOfArgumentsException(boost::intrusive_ptr<const AST::ParsedObject> po, int found, int expecting) :
		RuntimeException(formatExpectedArgs(found, expecting, expecting), po)
	{}
	~WrongNumberOfArgumentsException() throw () {}
};

class CorrespondingParameterNotReferenceException : public RuntimeException {
public:
	CorrespondingParameterNotReferenceException(boost::intrusive_ptr<const AST::ParsedObject> po) :
		RuntimeException("The corresponding parameter is not declared as a reference", po)
	{}
	~CorrespondingParameterNotReferenceException() throw () {}
};

class FieldNotFoundException : public RuntimeException {
	static std::string errorMessage(const std::string &fieldName, boost::intrusive_ptr<const RECORD> rv);
public:
	const std::string fieldName;
	FieldNotFoundException(const std::string &fieldName, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to, boost::intrusive_ptr<const RECORD> rv) :
		RuntimeException(errorMessage(fieldName, rv), from, to),
		fieldName(fieldName)
	{}
	FieldNotFoundException(const std::string &fieldName, boost::intrusive_ptr<const AST::ParsedObject> po, boost::intrusive_ptr<const RECORD> rv) :
		RuntimeException(errorMessage(fieldName, rv), po),
		fieldName(fieldName)
	{}
	FieldNotFoundException(const std::string &fieldName, const LexerNS::Token &token, boost::intrusive_ptr<const RECORD> rv) :
		RuntimeException(errorMessage(fieldName, rv), token),
		fieldName(fieldName)
	{}
	~FieldNotFoundException() throw () {}
};

class ProtectedVariableException : public RuntimeException {
public:
	const std::string protectedVarName;
	ProtectedVariableException(const std::string &protectedVarName, const std::string &errorMessage, boost::intrusive_ptr<const AST::ParsedObject> po) :
		RuntimeException(errorMessage, po),
		protectedVarName(protectedVarName)
	{}
	~ProtectedVariableException() throw () {}
};

class DivisionByZeroException : public RuntimeException {
public:
	DivisionByZeroException(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
		RuntimeException("Division by zero", from, to)
	{}
	~DivisionByZeroException() throw () {}
};

class RECORD;
class TYPE;

class Return {
	explicit Return(boost::intrusive_ptr<RightValue> value) :
		value(value)
	{}
	friend class AST::ReturnStatement;
public:
	const boost::intrusive_ptr<RightValue> value;
};

class Break {
	explicit Break(const std::string &label) : label(label) {}
	friend class AST::BreakStatement;
public:
	const std::string label;
};

class Continue {
	explicit Continue(const std::string &label) : label(label) {}
	friend class AST::ContinueStatement;
public:
	const std::string label;
};

class TaggedValue : public RightValue {
	boost::intrusive_ptr<TYPE> type;
public:
	const boost::intrusive_ptr<RightValue> value;
	const std::string tag;
	TaggedValue(boost::intrusive_ptr<RightValue> value, const std::string &tag);
	boost::intrusive_ptr<RightValue> untagged() { return this->value; }
	boost::intrusive_ptr<const RightValue> untagged() const { return this->value; }
	boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
	boost::intrusive_ptr<RightValue> clone();
	std::ostream &dumpAsString(std::ostream &out) const;
	bool needsToBeCopiedBeforeChanges() const;
	boost::intrusive_ptr<TYPE> getType() const;
	virtual ~TaggedValue() throw () {}
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
};

class ERROR : public ImmutableRightValue {
public:
	const std::string message;
	ERROR(std::string message) :
		message(message)
	{}
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<TYPE> getType() const;
	virtual ~ERROR() throw () {}
};

class OSTREAM : public ImmutableRightValue {
	bool canBeClosed;
	bool isClosed;
protected:
	OSTREAM(bool canBeClosed, bool isClosed) :
		canBeClosed(canBeClosed),
		isClosed(isClosed)
	{}
public:
	std::ostream &dumpAsString(std::ostream &out) const;
	virtual boost::intrusive_ptr<OSTREAM> print(const std::string &s) = 0;
	virtual boost::intrusive_ptr<OSTREAM> print(boost::intrusive_ptr<const RightValue> v) = 0;
	inline boost::intrusive_ptr<OSTREAM> println(const std::string &s) {
		return this->print(s)->newline();
	}
	inline boost::intrusive_ptr<OSTREAM> println(boost::intrusive_ptr<const RightValue> v) {
		return this->print(v)->newline();
	}
	virtual boost::intrusive_ptr<OSTREAM> newline() = 0;
	boost::intrusive_ptr<TYPE> getType() const;
	void print(RuntimeEnvironment *runtimeEnv, boost::intrusive_ptr<RightValue> v,  boost::intrusive_ptr<const AST::Expression> sourceExp);
	virtual boost::intrusive_ptr<RightValue> close(RuntimeEnvironment *runtimeEnv, boost::intrusive_ptr<const AST::Expression> sourceExp);
	virtual void flush() = 0;
	static boost::intrusive_ptr<TYPE> type;
	virtual ~OSTREAM() throw () {}
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
};

class CppOSTREAM : public OSTREAM {
public:
	CppOSTREAM(std::ostream &out, bool canBeClosed = false, bool isClosed = false) :
		OSTREAM(canBeClosed, isClosed),
		out(out)
	{}
	std::ostream &out;
	using OSTREAM::print;
	using OSTREAM::println;
	boost::intrusive_ptr<OSTREAM> print(const std::string &s);
	boost::intrusive_ptr<OSTREAM> print(boost::intrusive_ptr<const RightValue> v);
	boost::intrusive_ptr<OSTREAM> newline();
	void flush();
	virtual ~CppOSTREAM() throw () {}
};

//JAA 20140901 struct ZZModValue { static boost::intrusive_ptr<TYPE> type; };
struct RatFunValue { static boost::intrusive_ptr<TYPE> type; };
  //struct MODULEELEM { static boost::intrusive_ptr<TYPE> type; };
  //struct MODULE { static boost::intrusive_ptr<TYPE> type; };


class RINGELEM : public ImmutableRightValue {
public:
	RINGELEM(const RingElem &theRingElem) : theRingElem(theRingElem) {}
	const RingElem theRingElem;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
};

class RING : public ImmutableRightValue {
public:
	typedef std::pair<symbol, boost::intrusive_ptr<RINGELEM> > SymbolPair;
	const ring theRing;
	RING(const ring &theRing) :
		theRing(theRing)
	{}
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<TYPE> getType() const;
	void injectIndeterminates(RuntimeEnvironment *runtimeEnv, const boost::intrusive_ptr<const AST::RingDefinition> ringDefinition, const std::string &ringName);
	void removeInjectedIndeterminates(RuntimeEnvironment *runtimeEnv);
	inline int nIndeterminates() { return this->allSymbolValues().size(); }
	std::vector<SymbolPair> allSymbolValues();
private:
	std::vector<int> injectedSlots;
	static void createAllSymbolValues(const ring &r, std::vector<SymbolPair> &allSymbolValues);
};

class OutputStringStreamValue : public CppOSTREAM {
public:
	std::ostringstream ss;
	OutputStringStreamValue() :
		CppOSTREAM(ss, true, false)
	{}
	boost::intrusive_ptr<RightValue> close(RuntimeEnvironment *runtimeEnv, boost::intrusive_ptr<const AST::Expression> sourceExp);
	virtual ~OutputStringStreamValue() throw () {}
};

class OutputFileStreamValue : public CppOSTREAM {
public:
	std::ofstream stream;
	OutputFileStreamValue(boost::intrusive_ptr<STRING> filename, boost::intrusive_ptr<const AST::Expression> sourceExp);
	boost::intrusive_ptr<RightValue> close(RuntimeEnvironment *runtimeEnv, boost::intrusive_ptr<const AST::Expression> sourceExp);
	virtual ~OutputFileStreamValue() throw () {}
};

class VariableName : public LeftValue {
public:
	RuntimeEnvironment *runtimeEnvironment;
	FramePointer framePointer;
	int slotIndex;
	const boost::intrusive_ptr<const AST::Identifier> expId;
	VariableName(RuntimeEnvironment *runtimeEnvironment, Frame *frame, int slotIndex, boost::intrusive_ptr<const AST::Identifier> expId) :
		runtimeEnvironment(runtimeEnvironment),
		framePointer(frame),
		slotIndex(slotIndex),
		expId(expId)
	{
		assert(this->runtimeEnvironment);
		assert(this->framePointer.toCheckedPointer());
		assert(this->expId);
		assert(this->slotIndex>=AST::StaticEnv::TOP_LEVEL);
	}
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> asRightValue();
	void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv);
	boost::intrusive_ptr<VariableName> getBase();
	bool assignmentNeedsOwnership() const;
	virtual ~VariableName() throw () {}
};

class FieldAccess : public LeftValue {
	boost::intrusive_ptr<RECORD> targetAsRECORD();
public:
	const boost::intrusive_ptr<LeftValue> targetLV;
	const std::string fieldName;
	const boost::intrusive_ptr<const AST::FieldAccessExpression> fieldAccessExp;
	FieldAccess(const boost::intrusive_ptr<LeftValue> targetLV, const std::string &fieldName, boost::intrusive_ptr<const AST::FieldAccessExpression> fieldAccessExp) :
		targetLV(targetLV),
		fieldName(fieldName),
		fieldAccessExp(fieldAccessExp)
	{
		this->targetAsRECORD();
	}
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> asRightValue();
	void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv);
	boost::intrusive_ptr<VariableName> getBase();
	bool assignmentNeedsOwnership() const;
	virtual ~FieldAccess() throw () {}
};

class LIST : public RightValue {
public:
	typedef std::vector<boost::intrusive_ptr<RightValue> > ContainerType;
private:
	ContainerType container;
public:
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> clone();
	bool needsToBeCopiedBeforeChanges() const;
	boost::intrusive_ptr<RightValue> getValue(ContainerType::size_type index) const;
	void setValue(ContainerType::size_type index, boost::intrusive_ptr<RightValue> newValue);
	void addValue(boost::intrusive_ptr<RightValue> newValue);
	void insertValue(ContainerType::size_type position, boost::intrusive_ptr<RightValue> newValue);
	void removeValue(ContainerType::size_type position);
	ContainerType::size_type size() const;
	LIST() {}
	template <typename T> explicit LIST(const std::vector<T> &v);
	explicit LIST(const ContainerType::size_type reserveSize) { this->container.reserve(reserveSize); }
	ContainerType::size_type checkIndex(const boost::intrusive_ptr<const RightValue> v, const boost::intrusive_ptr<const AST::Expression> originalExp) const;
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
};

class MAT : public RightValue {
	MAT(const MAT &that) : RightValue(), theMatrix(that.theMatrix) {}
public:
	matrix theMatrix;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> clone();
	bool needsToBeCopiedBeforeChanges() const;
	std::size_t numRows() const { return NumRows(this->theMatrix); }
	std::size_t numColumns() const { return NumCols(this->theMatrix); }
	MAT(const matrix &theMatrix) : theMatrix(theMatrix) {}
	MAT(ConstMatrixView theMatrix) : theMatrix(NewDenseMat(theMatrix)) {}
	MAT(boost::intrusive_ptr<RING> ring, std::size_t nRows, std::size_t nColumns) : theMatrix(NewDenseMat(ring->theRing, nRows, nColumns)) {}
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
};

class IDEAL : public ImmutableRightValue {
public:
	ideal theIdeal;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	IDEAL(const ideal &theIdeal) : theIdeal(theIdeal) {}
	boost::intrusive_ptr<TYPE> getType() const;
};

class MODULE : public ImmutableRightValue {
public:
	module theModule;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	MODULE(const module &theModule) : theModule(theModule) {}
	boost::intrusive_ptr<TYPE> getType() const;
};

class MODULEELEM : public ImmutableRightValue {
public:
	MODULEELEM(const ModuleElem &theModuleElem) : theModuleElem(theModuleElem) {}
	ModuleElem theModuleElem;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
};

class JBMillValue : public ImmutableRightValue {
public:
	JBMill theJBMill;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	JBMillValue(const JBMill &theJBMill) : theJBMill(theJBMill) {}
	boost::intrusive_ptr<TYPE> getType() const;
};



class MatrixRowValue : public ImmutableRightValue {
public:
	const boost::intrusive_ptr<MAT> matrix;
	const std::size_t nRow;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	MatrixRowValue(boost::intrusive_ptr<MAT> matrix, std::size_t nRow) : matrix(matrix), nRow(nRow) {}
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
};

class IntMapValue : public ImmutableRightValue {
public:
	typedef long KeyType;
	typedef std::map<KeyType, boost::intrusive_ptr<RightValue> > MapType;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> getValue(KeyType index) const {
		MapType::const_iterator it = this->map.find(index);
		if (it==this->map.end())
			return 0;
		return it->second;
	}
	inline void addValue(KeyType index, boost::intrusive_ptr<RightValue> v) {
		assert(this->map.find(index)==this->map.end());
		this->map.insert(std::make_pair(index, v));
	}
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
private:
	MapType map;
};

class IndexedAccess : public LeftValue {
public:
	const boost::intrusive_ptr<LeftValue> targetLV;
	boost::intrusive_ptr<VariableName> getBase();
	bool assignmentNeedsOwnership() const;
	const LexerNS::CharPointer targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd;
	explicit IndexedAccess(const boost::intrusive_ptr<LeftValue> targetLV, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd) :
		targetLV(targetLV),
		targetExpBegin(targetExpBegin),
		targetExpEnd(targetExpEnd),
		indexExpBegin(indexExpBegin),
		indexExpEnd(indexExpEnd)
	{
	}
	virtual ~IndexedAccess() throw () {}
};

class StringIndexedAccess : public IndexedAccess {
	boost::intrusive_ptr<RECORD> targetAsRECORD();
public:
	const std::string fieldName;
	StringIndexedAccess(const boost::intrusive_ptr<LeftValue> targetLV, const std::string &fieldName,  const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd) :
		IndexedAccess(targetLV, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd),
		fieldName(fieldName)
	{
		this->targetAsRECORD();
		if (!LexerNS::Lexer::isValidIdentifier(fieldName))
			throw RuntimeException("Invalid identifier", indexExpBegin, indexExpEnd);
	}
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> asRightValue();
	void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv);
	virtual ~StringIndexedAccess() throw () {}
};

class IntegerIndexedAccess : public IndexedAccess {
	boost::intrusive_ptr<RightValue> targetAsListOrMatrixOrMatrixRow(int &which);
	LIST::ContainerType::size_type getIndexFromBigInt(const BigInt &N, boost::intrusive_ptr<const LIST> list) const;
	std::string::size_type getIndexFromBigInt(const BigInt &N, boost::intrusive_ptr<const STRING> str) const;
	std::size_t getIndexFromBigInt(const BigInt &N, boost::intrusive_ptr<const MatrixRowValue> mrv) const;
public:
	const BigInt index;
	IntegerIndexedAccess(const boost::intrusive_ptr<LeftValue> targetLV, const BigInt &index,  const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> asRightValue();
	void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv);
	virtual ~IntegerIndexedAccess() throw () {}
};

class ReferenceVariable : public LeftValue {
	const boost::intrusive_ptr<const AST::Identifier> expId;
	FramePointer framePointer;
	int slotIndex;
	void throwIncompatibleChanges() const;
	VariableSlot *referencedSlot() const;
	inline boost::intrusive_ptr<LeftValue> referencedLV() const { return intrusive_ptr_cast<LeftValue>(this->referencedSlot()->value); }
public:
	ReferenceVariable(const boost::intrusive_ptr<const AST::Identifier> expId, FramePointer framePointer, int slotIndex) :
		expId(expId),
		framePointer(framePointer),
		slotIndex(slotIndex)
	{
		assert(this->expId);
		assert(this->framePointer.toCheckedPointer());
		assert(this->slotIndex>=0);
		assert(static_cast<std::vector<VariableSlot>::size_type>(this->slotIndex) < this->framePointer.toCheckedPointer()->varSlots.size());
	}
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> asRightValue();
	void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv);
	boost::intrusive_ptr<VariableName> getBase();
	bool assignmentNeedsOwnership() const;
	virtual ~ReferenceVariable() throw () {}
};

class BOOL : public ImmutableRightValue {
	explicit BOOL(const bool theBool) : theBool(theBool) {}
public:
	static boost::intrusive_ptr<BOOL> trueValue, falseValue;
	inline static boost::intrusive_ptr<BOOL> fromBool(const bool b) {
		return b ? trueValue : falseValue;
	}
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	const bool theBool;
	boost::intrusive_ptr<TYPE> getType() const;
	virtual ~BOOL() throw () {}
};

class FUNCTION : public ImmutableRightValue {
public:
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	virtual bool canBeCalledWith(int nArg) = 0;
	virtual boost::intrusive_ptr<RightValue> eval(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression, RuntimeEnvironment * const) const = 0;
	boost::intrusive_ptr<TYPE> getType() const;
	virtual ~FUNCTION() throw () {}
};

class BuiltInFunction : public FUNCTION {
	typedef boost::intrusive_ptr<RightValue> (*FunType)(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression, RuntimeEnvironment * const);
	typedef bool (*ArityCheckType)(int arity);
	FunType function;
	ArityCheckType arityCheck;
public:
	BuiltInFunction(FunType function, ArityCheckType arityCheck) :
		function(function),
		arityCheck(arityCheck)
	{}
	bool canBeCalledWith(int nArg);
	boost::intrusive_ptr<RightValue> eval(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression, RuntimeEnvironment * const) const;
	void describe(boost::intrusive_ptr<OSTREAM> out) const;
	virtual ~BuiltInFunction() throw () {}
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
};

class RINGHOM : public FUNCTION {
public:
	const RingHom theRingHom;
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	RINGHOM(const RingHom &theRingHom) : theRingHom(theRingHom)
	{}
	bool canBeCalledWith(int nArg);
	boost::intrusive_ptr<RightValue> eval(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression, RuntimeEnvironment * const) const;
	void describe(boost::intrusive_ptr<OSTREAM> out) const;
	boost::intrusive_ptr<TYPE> getType() const;
	virtual ~RINGHOM() throw () {}
};

class UserDefinedFunction : public FUNCTION {
public:
	const boost::intrusive_ptr<const AST::FunctionDeclaration> fnDecl;
	const FramePointer fp;
	std::vector<boost::intrusive_ptr<RightValue> > capturedValues;
	UserDefinedFunction(RuntimeEnvironment *runtimeEnv, const boost::intrusive_ptr<const AST::FunctionDeclaration> fnDecl);
	boost::intrusive_ptr<RightValue> eval(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression, RuntimeEnvironment * const) const;
	bool canBeCalledWith(int nArg);
	void describe(boost::intrusive_ptr<OSTREAM> out) const;
	virtual ~UserDefinedFunction() throw () {}
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
};

class STRING : public ImmutableRightValue {
public:
	explicit STRING(const std::string &theString) : theString(theString) {}
	static boost::intrusive_ptr<TYPE> type;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string theString;
	static boost::intrusive_ptr<STRING> empty;
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByBigInt() { return true; }
};

class RECORD : public RightValue {
public:
	static boost::intrusive_ptr<TYPE> type;
	typedef std::map<std::string, boost::intrusive_ptr<RightValue> > MapType;
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<RightValue> clone();
	bool needsToBeCopiedBeforeChanges() const;
	boost::intrusive_ptr<RightValue> getField(const std::string &fieldName) const;
	virtual void setField(const std::string &fieldName, boost::intrusive_ptr<RightValue> newValue);
	MapType::size_type numberOfFields() const;
	MapType::const_iterator begin() const;
	MapType::const_iterator end() const;
	void collectSimilarlyNamedFields(const std::string &id, std::set<std::string> &set) const;
  std::vector<std::string> myFieldNamesStrings() const;
	boost::intrusive_ptr<LIST> fieldNames() const;
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<Value> indexedByString(boost::intrusive_ptr<STRING> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	bool canBeIndexedByString() { return true; }
#ifdef C5IDE
	void initTreeWidget(QTreeWidgetItem *) const;
#endif // #ifdef C5IDE
private:
	MapType fields;
};

class PackageValue : public ImmutableRightValue {
	friend void AST::PackageStatement::implExecute(InterpreterNS::RuntimeEnvironment *);
	PackageValue(const boost::intrusive_ptr<const AST::PackageStatement> pkgDecl);
	std::vector<int> ownedTopLevelIndexes;
public:
	const boost::intrusive_ptr<const AST::PackageStatement> pkgDecl;
	static boost::intrusive_ptr<TYPE> type;
	const std::string pkgName;
	const std::string prefix; // that is, the name plus the "."
	std::ostream &dumpAsString(std::ostream &out) const { assert(false); return out; }
	boost::intrusive_ptr<RightValue> clone() { assert(false); return this; }
	bool needsToBeCopiedBeforeChanges() const { return false; }
	boost::intrusive_ptr<TYPE> getType() const;
	boost::intrusive_ptr<VariableName> toVariableName(RuntimeEnvironment *runtimeEnv, const std::string &memberName, const LexerNS::Token &tokMemberName);
	void describe(boost::intrusive_ptr<OSTREAM> out) const;
	virtual ~PackageValue() throw () {}
};

std::ostream &operator<<(std::ostream &out, const Value * const value);

class IndexOutOfRangeException : public RuntimeException {
	typedef LIST::ContainerType::size_type size_type;
	static std::string errorMessage(const BigInt &foundIndex, size_type numOfElements) {
		return numOfElements==0 ?
				"Cannot index an empty list" :
				"The index "+ boost::lexical_cast<std::string>(foundIndex) + " is not in the range 1.." + boost::lexical_cast<std::string>(numOfElements);
	}
public:
	const BigInt foundIndex;
	const size_type numOfElements;
	IndexOutOfRangeException(BigInt foundIndex, size_type numOfElements, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
		RuntimeException(errorMessage(foundIndex, numOfElements), from, to),
		foundIndex(foundIndex),
		numOfElements(numOfElements)
	{}
	IndexOutOfRangeException(BigInt foundIndex, size_type numOfElements, boost::intrusive_ptr<const AST::ParsedObject> po) :
		RuntimeException(errorMessage(foundIndex, numOfElements), po),
		foundIndex(foundIndex),
		numOfElements(numOfElements)
	{}
	IndexOutOfRangeException(BigInt foundIndex, size_type numOfElements, const LexerNS::Token &token) :
		RuntimeException(errorMessage(foundIndex, numOfElements), token),
		foundIndex(foundIndex),
		numOfElements(numOfElements)
	{}
	~IndexOutOfRangeException() throw () {}
};

class NonIntegerIndexableException : public RuntimeException {
public:
	NonIntegerIndexableException(boost::intrusive_ptr<const RightValue> value, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
		RuntimeException(value->getType()->name+" cannot be indexed by integers", from, to)
	{}
	~NonIntegerIndexableException() throw () {}
};

class NonStringIndexableException : public RuntimeException {
public:
	NonStringIndexableException(boost::intrusive_ptr<const RightValue> value, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
		RuntimeException(value->getType()->name+" cannot be indexed by strings", from, to)
	{}
	~NonStringIndexableException() throw () {}
};

class InterruptException : public RuntimeException {
public:
	explicit InterruptException(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
		RuntimeException("\n\n\n*** Interrupted ***\n\n\n", from, to)
	{}
	~InterruptException() throw () {}
};

class CannotUseEllipsisInThisContextException : public RuntimeException {
public:
	explicit CannotUseEllipsisInThisContextException(const boost::intrusive_ptr<const AST::ParsedObject> po) :
		RuntimeException("Ellipsis can be only used inside functions/procedures with a variable number of arguments", po)
	{}
	~CannotUseEllipsisInThisContextException() throw () {}
};

class Ciao {};


class BinaryOpDispatcher : boost::noncopyable {
	static BinaryOpDispatcher uselessObjectToForceMapInitialization;
public:
	static inline bool isValidIndex(int i) { return i>=0 && i<TYPE::NUM_OF_DISPATCHABLE_TYPES; }
	BinaryOpDispatcher();
};

class Interpreter;

template <typename T> struct CoCoALibType {};
template <> struct CoCoALibType<RAT> { typedef BigRat type; };
template <> struct CoCoALibType<const RAT> { typedef const BigRat type; };
template <> struct CoCoALibType<INT> { typedef BigInt type; };
template <> struct CoCoALibType<const INT> { typedef const BigInt type; };
template <> struct CoCoALibType<RINGELEM> { typedef RingElem type; };
template <> struct CoCoALibType<const RINGELEM> { typedef const RingElem type; };
template <> struct CoCoALibType<RING> { typedef ring type; };
template <> struct CoCoALibType<const RING> { typedef const ring type; };
template <> struct CoCoALibType<RINGHOM> { typedef RingHom type; };
template <> struct CoCoALibType<const RINGHOM> { typedef const RingHom type; };
template <> struct CoCoALibType<MAT> { typedef matrix type; };
template <> struct CoCoALibType<const MAT> { typedef const matrix type; };
template <> struct CoCoALibType<STRING> { typedef std::string type; };
template <> struct CoCoALibType<const STRING> { typedef const std::string type; };
template <> struct CoCoALibType<TYPE> { typedef std::string type; };
template <> struct CoCoALibType<const TYPE> { typedef const std::string type; };
template <> struct CoCoALibType<BOOL> { typedef bool type; };
template <> struct CoCoALibType<const BOOL> { typedef bool type; };
template <> struct CoCoALibType<IDEAL> { typedef ideal type; };
template <> struct CoCoALibType<const IDEAL> { typedef const ideal type; };
template <> struct CoCoALibType<MODULE> { typedef module type; };
template <> struct CoCoALibType<const MODULE> { typedef const module type; };
template <> struct CoCoALibType<MODULEELEM> { typedef ModuleElem type; };
template <> struct CoCoALibType<const MODULEELEM> { typedef const ModuleElem type; };
template <> struct CoCoALibType<JBMillValue> { typedef JBMill type; };
template <> struct CoCoALibType<const JBMillValue> { typedef const JBMill type; };



template <typename T> typename CoCoALibType<T>::type theValue(boost::intrusive_ptr<T> x);
template <> BigInt theValue<INT>(boost::intrusive_ptr<INT> x);
template <> const BigInt theValue<const INT>(boost::intrusive_ptr<const INT> x);
template <> BigRat theValue<RAT>(boost::intrusive_ptr<RAT> x);
template <> const BigRat theValue<const RAT>(boost::intrusive_ptr<const RAT> x);
template <> RingElem theValue<RINGELEM>(boost::intrusive_ptr<RINGELEM> x);
template <> const RingElem theValue<const RINGELEM>(boost::intrusive_ptr<const RINGELEM> x);
template <> ring theValue<RING>(boost::intrusive_ptr<RING> x);
template <> const ring theValue<const RING>(boost::intrusive_ptr<const RING> x);
template <> RingHom theValue<RINGHOM>(boost::intrusive_ptr<RINGHOM> x);
template <> const RingHom theValue<const RINGHOM>(boost::intrusive_ptr<const RINGHOM> x);
template <> matrix theValue<MAT>(boost::intrusive_ptr<MAT> x);
template <> const matrix theValue<const MAT>(boost::intrusive_ptr<const MAT> x);
template <> std::string theValue<STRING>(boost::intrusive_ptr<STRING> x);
template <> const std::string theValue<const STRING>(boost::intrusive_ptr<const STRING> x);
template <> std::string theValue<TYPE>(boost::intrusive_ptr<TYPE> x);
template <> const std::string theValue<const TYPE>(boost::intrusive_ptr<const TYPE> x);
template <> bool theValue<BOOL>(boost::intrusive_ptr<BOOL> x);
template <> bool theValue<const BOOL>(boost::intrusive_ptr<const BOOL> x);
template <> ideal theValue<IDEAL>(boost::intrusive_ptr<IDEAL> x);
template <> const ideal theValue<const IDEAL>(boost::intrusive_ptr<const IDEAL> x);
template <> module theValue<MODULE>(boost::intrusive_ptr<MODULE> x);
template <> const module theValue<const MODULE>(boost::intrusive_ptr<const MODULE> x);
template <> ModuleElem theValue<MODULEELEM>(boost::intrusive_ptr<MODULEELEM> x);
template <> const ModuleElem theValue<const MODULEELEM>(boost::intrusive_ptr<const MODULEELEM> x);
template <> JBMill theValue<JBMillValue>(boost::intrusive_ptr<JBMillValue> x);
template <> const JBMill theValue<const JBMillValue>(boost::intrusive_ptr<const JBMillValue> x);



class RuntimeEnvironment {
  std::vector<Frame> frames; // JAA
	Frame *currentFrame;
	FrameIdentifier nextFrameId;
	std::map<std::string, int> topLevelIdentifiers;
	void initBuiltInFunctions();
	void initMaps();
	void registerType(boost::intrusive_ptr<TYPE> type);
	boost::intrusive_ptr<OSTREAM> standardOutput;
	std::set<std::string> registeredTypeNames;
	template <class T> static boost::intrusive_ptr<T> evalArgAsCheckType(const AST::Argument &arg, const boost::intrusive_ptr<RightValue> v);
public:
  long ResizeStack(long NewSize, const boost::intrusive_ptr<const AST::Expression> OrigExp);
	enum EvalKind { EVAL_BY_REF, EVAL_BY_VALUE };
	bool isRegisteredType(const std::string &typeName) { return this->registeredTypeNames.find(typeName)!=this->registeredTypeNames.end(); }
	const boost::intrusive_ptr<AST::StaticEnv> topLevelStaticEnv;
	const boost::intrusive_ptr<AST::AliasEnv> topLevelAliases;
	const int itSlot, currentRingSlot;
	Interpreter * const interpreter;
  RuntimeEnvironment(Interpreter * const interpreter,  boost::intrusive_ptr<OSTREAM> standardOutput, long MaxStackSize);
	int setTopLevelVar(const std::string &id, boost::intrusive_ptr<Value> v = VoidValue::theInstance, VariableSlot::Flags flags = VariableSlot::VSF_None);
	int slotFor(const std::string &id);
	inline Frame *getCurrentFrame() { return this->currentFrame; }
  /*JAA*/  inline Frame *getTopLevelFrame() { return &(this->frames[0]); }
  /*JAA*/  inline const Frame *getTopLevelFrame() const { return &(this->frames[0]); }
	Frame *getCurrentFunctionFrame();
	Frame *getCurrentFunctionFrameOrNull();
	Frame *pushFrame(const FramePointer &accessLink, boost::intrusive_ptr<const AST::InvocationExpression> invocationExp, boost::intrusive_ptr<const AST::ParsedObject> block, boost::intrusive_ptr<const UserDefinedFunction> userdefinedFun);
	Frame *pushIterationFrame(boost::intrusive_ptr<const AST::ParsedObject> block);
	void popFrame();
	std::vector<SnapshotFrame> takeSnapshot();
        bool CheckArity(int arity, std::map<std::string, int>::const_iterator pos) const; // added JAA 20140903
	void collectSimilarlyNamedIdentifiers(const std::string &id, int arity, std::vector<std::string> &NearMatches, bool &thereIsAnExactMatch) const;
	boost::intrusive_ptr<const AST::StaticEnv> findStaticEnv();
	boost::intrusive_ptr<Value> evalArg(const AST::Argument &arg, EvalKind evalKind);
	template <typename T> boost::intrusive_ptr<T> evalArgAs(const AST::Argument &arg);
	template <typename T> std::vector<typename CoCoALibType<T>::type > evalArgAsListOf(const AST::Argument &arg);
	template <typename T> std::vector<typename CoCoALibType<T>::type > evalAllArgsAsListOf(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression);
	std::vector<RingElem> evalArgAsListOfRingElem(const AST::Argument &arg, const ring &ring);
	std::vector<RingElem> evalArgAsListOfRingElem(const AST::Argument &arg);
	template <typename T1, typename T2> boost::intrusive_ptr<RightValue> evalArgAsT1orT2(const AST::Argument &arg, bool &isT1); // not implemented
	template <typename T1, typename T2> boost::intrusive_ptr<RightValue> evalArgAsT1orT2(const AST::Argument &arg, int &which);
	template <typename T1, typename T2, typename T3> boost::intrusive_ptr<RightValue> evalArgAsT1orT2orT3(const AST::Argument &arg, int &which);
	template <typename T1, typename T2, typename T3, typename T4> boost::intrusive_ptr<RightValue> evalArgAsT1orT2orT3orT4(const AST::Argument &arg, int &which);
	template <typename T1, typename T2, typename T3, typename T4, typename T5> boost::intrusive_ptr<RightValue> evalArgAsT1orT2orT3orT4orT5(const AST::Argument &arg, int &which);
	template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6> boost::intrusive_ptr<RightValue> evalArgAsT1orT2orT3orT4orT5orT6(const AST::Argument &arg, int &which);
	template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7> boost::intrusive_ptr<RightValue> evalArgAsT1orT2orT3orT4orT5orT6orT7(const AST::Argument &arg, int &which);
	template <typename T> boost::intrusive_ptr<T> obtainUnshared(const AST::Argument &arg);
	boost::intrusive_ptr<OSTREAM> getOutputStream();
	typedef boost::intrusive_ptr<RightValue> (*BinaryOpFunc)(RuntimeEnvironment *, const boost::intrusive_ptr<const RightValue> leftOp, const boost::intrusive_ptr<const RightValue> rightOp,  const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition);
	typedef BinaryOpFunc DispatchMapType[TYPE::NUM_OF_DISPATCHABLE_TYPES][TYPE::NUM_OF_DISPATCHABLE_TYPES];
	static DispatchMapType opPlusMap;
	static DispatchMapType opMinusMap;
	static DispatchMapType opStarMap;
	static DispatchMapType opSlashMap;
	static DispatchMapType opLessThanMap;
	static DispatchMapType opLessOrEqualMap;
	static DispatchMapType opEqualMap;
	static DispatchMapType opNotEqualMap;
	static DispatchMapType opGreaterOrEqualMap;
	static DispatchMapType opGreaterThanMap;
	static DispatchMapType opPowerMap;
	static DispatchMapType opModMap;
	static DispatchMapType opColonMap;
	boost::intrusive_ptr<RightValue> binaryOperatorDispatch(boost::intrusive_ptr<const RightValue> leftOp, boost::intrusive_ptr<const RightValue> rightOp, DispatchMapType map, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition);
	boost::intrusive_ptr<RING> currentRing;
	boost::intrusive_ptr<LIST> topLevelFunctions();
	boost::intrusive_ptr<LIST> currentTypes();
	void announceCLE(const ErrorInfo& err);
#ifdef C5IDE
	friend class IDE::Debugger;
	void aboutToExecute(boost::intrusive_ptr<RuntimeEnvironment> runtimeEnv, boost::intrusive_ptr<AST::Statement> stmt);
	void aboutToEvaluate(boost::intrusive_ptr<RuntimeEnvironment> runtimeEnv, boost::intrusive_ptr<AST::Expression> exp);
#endif // #ifdef C5IDE
};

template <typename T>
boost::intrusive_ptr<T> RuntimeEnvironment::obtainUnshared(const AST::Argument &arg) {
	boost::intrusive_ptr<LeftValue> lv = intrusive_ptr_cast<LeftValue>(this->evalArg(arg, EVAL_BY_REF));
	boost::intrusive_ptr<RightValue> rv = lv->asRightValue();
	boost::intrusive_ptr<T> value = boost::dynamic_pointer_cast<T>(rv);
	if (!value)
		throw WrongTypeException(
				T::type->name,
				rv->getType()->name,
				arg.exp);
	rv = 0;
	value = 0;
	lv->obtainOwnership();
	value = intrusive_ptr_cast<T>(lv->asRightValue());
	return value;
}

template <typename T> inline boost::intrusive_ptr<T> RuntimeEnvironment::evalArgAsCheckType(const AST::Argument &arg, const boost::intrusive_ptr<RightValue> v) {
	if (const boost::intrusive_ptr<T> x = boost::dynamic_pointer_cast<T>(v))
		return x;
	throw WrongTypeException(T::type->name, v->getType()->name, arg.exp);
}

template <> inline boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsCheckType(const AST::Argument &, const boost::intrusive_ptr<RightValue> v) {
	return v;
}

template <typename T> inline boost::intrusive_ptr<T> RuntimeEnvironment::evalArgAs(const AST::Argument &arg) {
	return evalArgAsCheckType<T>(arg, intrusive_ptr_cast<RightValue>(this->evalArg(arg, EVAL_BY_VALUE)));
}

template <typename T> std::vector<typename CoCoALibType<T>::type > RuntimeEnvironment::evalArgAsListOf(const AST::Argument &arg) {
	const boost::intrusive_ptr<LIST> list = this->evalArgAs<LIST>(arg);
	const LIST::ContainerType::size_type size = list->size();
	std::vector<typename CoCoALibType<T>::type > result;
	for(LIST::ContainerType::size_type a=0; a<size; ++a) {
		if (const boost::intrusive_ptr<T> elem = boost::dynamic_pointer_cast<T>(list->getValue(a)))
			result.push_back(theValue(elem));
		else
			throw RuntimeException("List of "+T::type->name+" is required", arg.exp);
	}
	return result;
}

template <typename T> std::vector<typename CoCoALibType<T>::type > RuntimeEnvironment::evalAllArgsAsListOf(const boost::intrusive_ptr<const AST::InvocationExpression> invocationExpression) {
	const size_t nArgs = invocationExpression->args.size();
	std::vector<typename CoCoALibType<T>::type > result;
	if (nArgs==0)
		return result;
	const boost::intrusive_ptr<RightValue> firstArg = this->evalArgAs<RightValue>(invocationExpression->args.front());
	if (nArgs==1) {
		if (const boost::intrusive_ptr<LIST> list = boost::dynamic_pointer_cast<LIST>(firstArg)) {
			const LIST::ContainerType::size_type size = list->size();
			for(LIST::ContainerType::size_type a=0; a<size; ++a) {
				if (const boost::intrusive_ptr<T> elem = boost::dynamic_pointer_cast<T>(list->getValue(a)))
					result.push_back(theValue(elem));
				else
					throw RuntimeException("List of "+T::type->name+" is required", invocationExpression->args.front().exp);
			}
			return result;
		}
	}
	if (const boost::intrusive_ptr<T> v = boost::dynamic_pointer_cast<T>(firstArg))
		result.push_back(theValue(v));
	else
		throw WrongTypeException(T::type->name, firstArg->getType()->name, invocationExpression->args.front().exp);
	for(size_t a=1; a<nArgs; ++a)
		result.push_back(theValue(this->evalArgAs<T>(invocationExpression->args[a])));
	return result;
}

// not implemented to prevent improper use with bool (obsolete)
// template <typename T1, typename T2> boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2(const AST::Argument &arg, bool &isT1)

template <typename T1, typename T2> boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	throw WrongTypeException(T1::type->name+" or "+T2::type->name, v->getType()->name, arg.exp);
}

template <typename T1, typename T2, typename T3> boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2orT3(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	if (const boost::intrusive_ptr<T3> x = boost::dynamic_pointer_cast<T3>(v)) {
		which = 3;
		return x;
	}
	throw WrongTypeException(T1::type->name+", "+T2::type->name+" or "+T3::type->name, v->getType()->name, arg.exp);
}

  template <typename T1, typename T2, typename T3, typename T4> boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2orT3orT4(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	if (const boost::intrusive_ptr<T3> x = boost::dynamic_pointer_cast<T3>(v)) {
		which = 3;
		return x;
	}
	if (const boost::intrusive_ptr<T4> x = boost::dynamic_pointer_cast<T4>(v)) {
		which = 4;
		return x;
	}
	throw WrongTypeException(T1::type->name+", "+T2::type->name+", "+T3::type->name+" or "+T4::type->name, v->getType()->name, arg.exp);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5> boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2orT3orT4orT5(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	if (const boost::intrusive_ptr<T3> x = boost::dynamic_pointer_cast<T3>(v)) {
		which = 3;
		return x;
	}
	if (const boost::intrusive_ptr<T4> x = boost::dynamic_pointer_cast<T4>(v)) {
		which = 4;
		return x;
	}
	if (const boost::intrusive_ptr<T5> x = boost::dynamic_pointer_cast<T5>(v)) {
		which = 5;
		return x;
	}
	throw WrongTypeException(T1::type->name+", "+T2::type->name+", "+T3::type->name+", "+T4::type->name+" or "+T5::type->name, v->getType()->name, arg.exp);
}

  template <typename T1, typename T2, typename T3,
            typename T4, typename T5, typename T6> 
  boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2orT3orT4orT5orT6(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	if (const boost::intrusive_ptr<T3> x = boost::dynamic_pointer_cast<T3>(v)) {
		which = 3;
		return x;
	}
	if (const boost::intrusive_ptr<T4> x = boost::dynamic_pointer_cast<T4>(v)) {
		which = 4;
		return x;
	}
	if (const boost::intrusive_ptr<T5> x = boost::dynamic_pointer_cast<T5>(v)) {
		which = 5;
		return x;
	}
	if (const boost::intrusive_ptr<T6> x = boost::dynamic_pointer_cast<T6>(v)) {
		which = 6;
		return x;
	}
	throw WrongTypeException(T1::type->name+", "+T2::type->name
                           +", "+T3::type->name+", "+T4::type->name
                           +", "+T5::type->name+" or "+T6::type->name,
                           v->getType()->name, arg.exp);
}

  template <typename T1, typename T2, typename T3,
            typename T4, typename T5, typename T6, typename T7> 
  boost::intrusive_ptr<RightValue> RuntimeEnvironment::evalArgAsT1orT2orT3orT4orT5orT6orT7(const AST::Argument &arg, int &which) {
	const boost::intrusive_ptr<RightValue> v = this->evalArgAs<RightValue>(arg);
	if (const boost::intrusive_ptr<T1> x = boost::dynamic_pointer_cast<T1>(v)) {
		which = 1;
		return x;
	}
	if (const boost::intrusive_ptr<T2> x = boost::dynamic_pointer_cast<T2>(v)) {
		which = 2;
		return x;
	}
	if (const boost::intrusive_ptr<T3> x = boost::dynamic_pointer_cast<T3>(v)) {
		which = 3;
		return x;
	}
	if (const boost::intrusive_ptr<T4> x = boost::dynamic_pointer_cast<T4>(v)) {
		which = 4;
		return x;
	}
	if (const boost::intrusive_ptr<T5> x = boost::dynamic_pointer_cast<T5>(v)) {
		which = 5;
		return x;
	}
	if (const boost::intrusive_ptr<T6> x = boost::dynamic_pointer_cast<T6>(v)) {
		which = 6;
		return x;
	}
	if (const boost::intrusive_ptr<T7> x = boost::dynamic_pointer_cast<T7>(v)) {
		which = 7;
		return x;
	}
	throw WrongTypeException(T1::type->name+", "+T2::type->name
                           +", "+T3::type->name+", "+T4::type->name
                           +", "+T5::type->name+", "+T6::type->name
                           +" or "+T7::type->name,
                           v->getType()->name, arg.exp);
}

class VariableNotFoundException : public RuntimeException {
public:
	static std::string errorMessage(const std::string &varName, int arity, boost::intrusive_ptr<const AST::StaticEnv> env);
	const std::string varName;
	VariableNotFoundException(boost::intrusive_ptr<const AST::Identifier> id, RuntimeEnvironment *runtimeEnv) :
		RuntimeException(errorMessage(id->identifier, id->arity, runtimeEnv->findStaticEnv()), id),
		varName(id->identifier)
	{}
	~VariableNotFoundException() throw () {}
};

class Interpreter : public LexerNS::ReferenceCountedObject {
public:
	const bool warnAboutCocoa5;
	const bool fullCoCoALibError;
private:
	boost::intrusive_ptr<LexerNS::LineProvider> lineProvider;
	typedef std::pair<dev_t, ino_t> FilePair;
	std::set<FilePair> sourcedFiles;
	RuntimeEnvironment runtimeEnvironment;
#ifdef C5IDE
	friend class IDE::Debugger;
	boost::mutex mut;
	boost::condition_variable condVar;
	volatile bool mustResume;
#endif // #ifdef C5IDE
public:
	void checkForInterrupts(const boost::intrusive_ptr<const AST::ParsedObject> po);
	void checkForInterrupts(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	volatile bool controlC;
#ifdef C5IDE
	volatile bool singleStepExecution;
	volatile bool doStepOver;
	enum InterpreterStatus {
		IS_WAITING_FOR_COMMAND,
		IS_WAITING_FOR_COMMAND_COMPLETION,
		IS_RUNNING,
		IS_RUNNING_BUILTIN,
		IS_PAUSED,
		IS_ENDED
	};
	boost::intrusive_ptr<const AST::ParsedObject> pausedOn;
	inline InterpreterStatus getStatus() { return this->status; }
	inline void UpdateStatusToWFCC() {
		assert(this->status==IS_WAITING_FOR_COMMAND || this->status==IS_WAITING_FOR_COMMAND_COMPLETION);
		this->status = IS_WAITING_FOR_COMMAND_COMPLETION;
	}
	inline InterpreterStatus UpdateStatusStartingBuiltIn() {
		InterpreterStatus prevStatus = this->status;
		assert(prevStatus==IS_RUNNING || prevStatus==IS_RUNNING_BUILTIN);
		this->status = IS_RUNNING_BUILTIN;
		return prevStatus;
	}
	inline void UpdateStatusEndingBuiltIn(InterpreterStatus prevStatus) {
		assert(prevStatus==IS_RUNNING || prevStatus==IS_RUNNING_BUILTIN);
		assert(this->status==IS_RUNNING_BUILTIN);
		this->status = prevStatus;
	}
	bool pause(boost::intrusive_ptr<const AST::ParsedObject> po);
	void resume();
	int run(IDE::Console *mainWindow);
#else // #ifdef C5IDE
	int run();
#endif // #ifdef C5IDE
	const boost::intrusive_ptr<LexerNS::ErrorReporter> errorReporter;
  Interpreter(bool warnAboutCocoa5, boost::intrusive_ptr<LexerNS::LineProvider> lineProvider, boost::intrusive_ptr<LexerNS::ErrorReporter> errorReporter, boost::intrusive_ptr<OSTREAM> standardOutput, bool fullCoCoALibError, long MaxStackSize);
  void readAndExecute(const std::string &filename, bool calledFromMain, bool immediateExecution, long FromLine=-1, long FromChar=-1, long ToLine=-1, long ToChar=-1);
	void reportWarning(const std::string &msg, const boost::intrusive_ptr<const AST::ParsedObject> po);
	void reportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void reportError(const std::string &msg);
	void reportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	static const std::string IT, CURRENT_RING;
#ifdef C5IDE
private:
	InterpreterStatus status;
#endif // #ifdef C5IDE
};

class ResolvePackageNamesVisitor : public AST::ParsedObjectVisitor {
public:
	const AST::PackageStatement &packageStatement;
	ResolvePackageNamesVisitor(const AST::PackageStatement &packageStatement) :
		packageStatement(packageStatement)
	{}
  using ParsedObjectVisitor::visit; // disables warnings of overloading
	void visit(AST::InvocationExpression &invocationExp);
	~ResolvePackageNamesVisitor() throw () {}
};

class CheckNamesVisitor : public AST::ParsedObjectVisitor {
	const boost::intrusive_ptr<LexerNS::ErrorReporter> errorReporter;
	boost::intrusive_ptr<AST::StaticEnv> env;
	boost::intrusive_ptr<AST::AliasEnv> aliases;
	boost::intrusive_ptr<AST::StaticEnv> addEnvironmentForIterationVar(const std::string &identifier, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	bool insideFunctions;
	bool foundErrors;
	RuntimeEnvironment *runtimeEnvironment;
	bool insidePackage;
	void popEnv() {
		assert(this->env);
		this->env = this->env->parent;
		assert(this->env);
	}
	void reportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
public:
	bool thereAreErrors() { return this->foundErrors; }
	CheckNamesVisitor(RuntimeEnvironment *runtimeEnvironment, const boost::intrusive_ptr<LexerNS::ErrorReporter> errorReporter) :
		errorReporter(errorReporter),
		env(runtimeEnvironment->topLevelStaticEnv),
		aliases(runtimeEnvironment->topLevelAliases),
		insideFunctions(false),
		foundErrors(false),
		runtimeEnvironment(runtimeEnvironment),
		insidePackage(false)
	{
		assert(this->runtimeEnvironment);
		assert(this->errorReporter);
	}
  using ParsedObjectVisitor::visit; // disables warnings of overloading
	void visit(AST::IdInExpSuchThatExp &idInExpSuchThatExp);
	void visit(AST::FieldAccessExpression &fieldAccessExp);
	void visit(AST::AliasStatement &aliasStmt);
	void visit(AST::ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp);
	void visit(AST::Identifier &);
	void visit(AST::ForStatement &forStmt);
	void visit(AST::ForeachStatement &foreachStmt);
	void visit(AST::DefineStatement &defineStmt);
	void visit(AST::FunctionDeclaration &funDecl);
	void visit(AST::PackageStatement &packageStmt);
	void visit(AST::InvocationExpression &invocationExp);
	void visit(AST::AssignmentStatement &assignmentStmt);
	void visit(AST::TryStatement &tryStatement);
	~CheckNamesVisitor() throw () {}
};

inline void Interpreter::checkForInterrupts(const boost::intrusive_ptr<const AST::ParsedObject> po) {
	if (this->controlC) {
		this->controlC = false;
		throw InterruptException(po->getBegin(), po->getEnd());
	}
}

inline void Interpreter::checkForInterrupts(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) {
	if (this->controlC) {
		this->controlC = false;
		throw InterruptException(from, to);
	}
}

inline boost::intrusive_ptr<RINGELEM> Value::from(const RingElem &x) { return new RINGELEM(x); }
inline boost::intrusive_ptr<RAT> Value::from(const BigRat &x) { return new RAT(x); }
inline boost::intrusive_ptr<INT> Value::from(const BigInt &x) { return new INT(x); }
inline boost::intrusive_ptr<INT> Value::from(const long x) { return new INT(x); }
// no Value::from(const size_t n) : CoCoALib never returns unsigned
template<typename T> inline boost::intrusive_ptr<LIST> Value::from(const std::vector<T> &v) { return new LIST(v); }
inline boost::intrusive_ptr<MAT> Value::from(const matrix &x) { return new MAT(x); }
inline boost::intrusive_ptr<MAT> Value::from(ConstMatrixView x) { return new MAT(x); }
inline boost::intrusive_ptr<BOOL> Value::from(const bool x) { return BOOL::fromBool(x); }
inline boost::intrusive_ptr<RINGHOM> Value::from(const RingHom &x) { return new RINGHOM(x); }
inline boost::intrusive_ptr<RING> Value::from(const ring &x) { return new RING(x); }
template <typename T> LIST::LIST(const std::vector<T> &v) {
	BOOST_FOREACH(const T &elem, v) {
		this->addValue(Value::from(elem));
	}
}
inline boost::intrusive_ptr<STRING> Value::from(const std::string &x) { return new STRING(x); }
inline boost::intrusive_ptr<IDEAL> Value::from(const ideal &x) { return new IDEAL(x); }
inline boost::intrusive_ptr<MODULE> Value::from(const module &x) { return new MODULE(x); }
inline boost::intrusive_ptr<MODULEELEM> Value::from(const ModuleElem &x) { return new MODULEELEM(x); }
inline boost::intrusive_ptr<JBMillValue> Value::from(const JBMill &theJBMill) { return new JBMillValue(theJBMill); }

template <> inline boost::intrusive_ptr<RightValue> Value::simplifiedValueFrom<BigRat>(const BigRat &value) {
	return (IsOneDen(value)) ? static_cast<RightValue *>(new INT(num(value))) : static_cast<RightValue *>(new RAT(value));
}

} // namespace InterpreterNS
} // namespace CoCoA

// NOTE: for some strange reason the following global consts are
//       defined in Main.C rather than Interpreter.C.
extern std::string packageDir;

extern const std::string PACKAGES_NOT_FOUND;
extern const std::string PACKAGE_AUTOLOAD_LOADING;
extern const std::string PACKAGE_AUTOLOAD_FAILURE_MESSAGE;
extern const std::string PACKAGE_AUTOLOAD_SKIPPING_PKG_DUE_TO_FAILURE;

#endif /* INTERPRETER_H_ */
