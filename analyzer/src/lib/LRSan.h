#ifndef _LRSAN_GLOBAL_H
#define _LRSAN_GLOBAL_H

#include <llvm/IR/DebugInfo.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Instructions.h>
#include <llvm/ADT/DenseMap.h>
#include <llvm/ADT/SmallPtrSet.h>
#include <llvm/ADT/StringExtras.h>
#include <llvm/Support/Path.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Analysis/AliasAnalysis.h>
#include "llvm/Support/CommandLine.h"
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Common.h"


// 
// typedefs
//
typedef std::vector< std::pair<llvm::Module*, llvm::StringRef> > ModuleList;
// Mapping module to its file name.
typedef std::unordered_map<llvm::Module*, llvm::StringRef> ModuleNameMap;
// The set of all functions.
typedef llvm::SmallPtrSet<llvm::Function*, 8> FuncSet;
// Mapping from function name to function.
typedef std::unordered_map<std::string, llvm::Function*> NameFuncMap;
typedef llvm::SmallPtrSet<llvm::CallInst*, 8> CallInstSet;
typedef llvm::DenseMap<llvm::Function*, CallInstSet> CallerMap;
typedef llvm::DenseMap<llvm::CallInst *, FuncSet> CalleeMap;
// Pointer analysis types.
typedef std::map<llvm::Value *, std::set<llvm::Value *>> PointerAnalysisMap;
typedef std::map<llvm::Function *, PointerAnalysisMap> FuncPointerAnalysisMap;
typedef std::map<llvm::Function *, AAResults *> FuncAAResultsMap;

struct GlobalContext {

  GlobalContext() {
    NumFunctions = 0;
  }

  unsigned NumFunctions;

	// Map global function name to function defination.
	NameFuncMap Funcs;

  // Functions whose addresses are taken.
  FuncSet AddressTakenFuncs;

	// Map a callsite to all potential callee functions.
	CalleeMap Callees;

	// Map a function to all potential caller instructions.
	CallerMap Callers;

  // Indirect call instructions.
  std::vector<CallInst *>IndirectCallInsts;

	// Modules.
  ModuleList Modules;
  ModuleNameMap ModuleMaps;
  std::set<std::string> InvolvedModules;

	std::map<std::string, uint8_t> MemWriteFuncs;
  std::set<std::string> CriticalFuncs;

  // Pinter analysis results.
  FuncPointerAnalysisMap FuncPAResults;
  FuncAAResultsMap FuncAAResults;
};

class IterativeModulePass {
protected:
	GlobalContext *Ctx;
	const char * ID;
public:
	IterativeModulePass(GlobalContext *Ctx_, const char *ID_)
		: Ctx(Ctx_), ID(ID_) { }

	// Run on each module before iterative pass.
	virtual bool doInitialization(llvm::Module *M)
		{ return true; }

	// Run on each module after iterative pass.
	virtual bool doFinalization(llvm::Module *M)
		{ return true; }

	// Iterative pass.
	virtual bool doModulePass(llvm::Module *M)
		{ return false; }

	virtual void run(ModuleList &modules);
};

#endif
