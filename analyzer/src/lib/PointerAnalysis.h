#ifndef _POINTER_ANALYSIS_H
#define _POINTER_ANALYSIS_H

#include "LRSan.h"

class PointerAnalysisPass : public IterativeModulePass {
  typedef std::pair<Value *, MemoryLocation *> AddrMemPair;

  private:
    TargetLibraryInfo *TLI;

    void detectAliasPointers(Function *, AAResults &,
                             PointerAnalysisMap &);

  public:
    PointerAnalysisPass(GlobalContext *Ctx_)
      : IterativeModulePass(Ctx_, "PointerAnalysis") { }
    virtual bool doInitialization(llvm::Module *);
    virtual bool doFinalization(llvm::Module *);
    virtual bool doModulePass(llvm::Module *);
};

#endif
