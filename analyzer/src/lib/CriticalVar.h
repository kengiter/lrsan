#ifndef _CRITICAL_VAR_H
#define _CRITICAL_VAR_H

#include <llvm/Analysis/BasicAliasAnalysis.h>

#include "LRSan.h"

class CriticalVarPass : public IterativeModulePass {

	typedef std::pair<TerminatorInst *, unsigned> CFGEdge;
	typedef std::pair<CFGEdge, Value *> EdgeValue;

	enum ErrnoFlag {
		Reserved,
		Must_Return_Errno,
		May_Return_Errno
	};

	enum UseDefFlag {
		Reserverd,
		Used,
		Defined,
	};

	typedef std::map<CFGEdge, ErrnoFlag> EdgeErrnoFlag;
	typedef std::pair<Instruction *, UseDefFlag> InstructionUseDef;

	private:

	unsigned int sc_counter; // counter of security checks
	unsigned int cuc_counter; // counter of check-use chains
	unsigned int lcc_counter; // counter of lacking-check cases
	std::set<Value *> cv_set; // set of critical variables

	std::set<Value *> TrackedAddr;

	// pairs for values and their checks; a value may correspond
	// multiple checks? 
	std::map<Value *, Instruction *>valueCheckPairs;

	// Map function name to function definition.
	NameFuncMap AllFuncs;

	// check if the value is a constant.
	bool isConstant(Value *V);

	// Check if the value is an errno.
	bool isValueErrno(Value *V);

	// Check if the value is a parameter of the function.
	bool isFunctionParameter(Value *V, Function *F);

	// Make the given edge with errno flag.
	void markEdge(CFGEdge &CE, ErrnoFlag flag, EdgeErrnoFlag &errnoEdges);

	// Mark the edge from predBB to succBB with errno flag.
	void markEdgeBlockToBlock(BasicBlock *predBB, BasicBlock *succBB, ErrnoFlag flag,
			EdgeErrnoFlag &errnoEdges);

	// Mark all edges from BB with errno flag.
	void markAllEdgesFromBlock(BasicBlock *BB, ErrnoFlag flag, EdgeErrnoFlag &errnoEdges);

	// Mark all edges to BB with errno flag.
	void markAllEdgesToBlock(BasicBlock *BB, ErrnoFlag flag, EdgeErrnoFlag &errnoEdges);

	// Add edges to analyze for a phinode.
	void addPHINodeEdges(PHINode *PN, std::list<EdgeValue> &EV);

	// Recursively do errno check.
	void recurCheckValueErrno(Function *F, std::list<EdgeValue> &EV,
			std::set<Value *> &PV, EdgeErrnoFlag &errnoEdges);

	// Check if the value must or may be an errno. This function is called recursively,
	void checkValueErrno(Function *F, Value *V, EdgeErrnoFlag &errnoEdges);

	// Check if the returned value must be or may be an errno.
	void checkReturnValue(Function *F, BasicBlock *BB, Value *V, EdgeErrnoFlag &errnoEdges);

	// Dump marked edges.
	void dumpEdges(EdgeErrnoFlag &errnoEdges);

	// Find security checks.
	void findSecurityChecks(Function *F, EdgeErrnoFlag &errnoEdges,
			std::set<Value *> &securityChecks);

	// Find critical variables.
	void findCriticalVariable(Function *, Value *, Value *, 
			std::map<Value *, std::set<Value *>> &, std::set<Value *> &);

	// Identify critical variables.
	void identifyCriticalVariables(Function *, std::set<Value *> &,
			std::map<Value *, std::set<Value *>> &);

	void filterNonGlobalVariables(std::map<Value *, std::set<Value *>> &);

	// Check alias result of two values.
	bool checkAlias(Value *, Value *,
			PointerAnalysisMap &);

	// Get aliased pointers for the pointer.
	void getAliasPointers(Value *, std::set<Value *> &,
			PointerAnalysisMap &);

	// Get called function name.
	StringRef getCalledFuncName(Value *);

	// Track aliased address after modification.
	void trackAliasAddrAfterModification(Value *,
			std::set<Value *> &,
			Instruction *,
			std::set<InstructionUseDef> &,
			PointerAnalysisMap &aliasPtrs);

	// Check whether the value is exploitable.
	bool isExploitable(Value *, Value *, std::map<Value *, bool> &,
			std::set<Value *> &);

	// Find modifications by function calls.
	void findUseDefFromCVAddr(Function *, Value *,
			std::set<Value *> &,
			std::set<InstructionUseDef> &,
			PointerAnalysisMap &);

	// Check if it is possible to use the value stored
	// by the store instruction.
	bool possibleUseStResult(Instruction *, Instruction *);

	// Find uses of the critical variable.
	void findUseDefOfCriticalVar(Value *CV, Value *cvAddr,
			std::set<Value *> &pSet,
			std::set<InstructionUseDef> &udSet,
			PointerAnalysisMap &aliasPtrs);

	// Filter out uses before definition.
	void filterUseBeforeDef(std::set<InstructionUseDef> &);

	// Filter definition before security check.
	void filterDefBeforeCheck(Value *, std::set<InstructionUseDef> &,
			Value *);

	void filterDefFromCheckedValue(Value *, std::set<InstructionUseDef> &,
			Value *);

	void trackCalledFunc(Function *, Value *, std::set<Function *> &,
			std::list<Instruction *> &,
			std::map<Instruction *, std::list<Instruction *>> &);

	Value *trackSrcOfVal(Function *, Value *, std::set<Value *> &);

	void trackSrc(Function *, Instruction *, Value *,
			std::set<Function *> &,
			std::list<Instruction *> &,
			std::map<Instruction *, std::list<Instruction *>> &);

	// Find out where the source comes from for the modification.
	void findSrcForModification(Function *, Instruction *,
			std::list<Instruction *> &,
			std::map<Instruction *, std::list<Instruction *>> &);

	// Print out source code information.
	void printSourceCodeInfo(Instruction *);
	void printSourceCodeInfo(Value *);

	// Find the dominating check
	Instruction * findDominatingCheck(Value *RV);

	// Find the checked variable or function
	Value * findCheckedValue(Instruction *CI);

	// Check if lacking-check issues may occur
	bool findLackingChecks(Value *V, Instruction *CI);

	public:
	CriticalVarPass(GlobalContext *Ctx_)
		: IterativeModulePass(Ctx_, "CriticalVar") { }
	virtual bool doInitialization(llvm::Module *);
	virtual bool doFinalization(llvm::Module *);
	virtual bool doModulePass(llvm::Module *);
};

#endif
