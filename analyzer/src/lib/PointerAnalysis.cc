#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/InstIterator.h>
#include <llvm/IR/LegacyPassManager.h>

#include "PointerAnalysis.h"

/// Alias types used to do pointer analysis.
#define MUST_ALIAS

bool PointerAnalysisPass::doInitialization(Module *M) {
    return false;
}

bool PointerAnalysisPass::doFinalization(Module *M) {
    return false;
}

/// Detect aliased pointers in this function.
void PointerAnalysisPass::detectAliasPointers(Function *F,
                                      AAResults &AAR,
                           PointerAnalysisMap &aliasPtrs) {
  std::set<AddrMemPair> addrSet;
  Value *Addr1, *Addr2;
  MemoryLocation *MemLoc1, *MemLoc2;
  std::set<Value*> DstSet, LPSet;

  DstSet.clear();
  LPSet.clear();
  addrSet.clear();
  aliasPtrs.clear();

  // Scan instructions to extract all pointers.
  for (inst_iterator i = inst_begin(F), ei = inst_end(F);
       i != ei; ++i) {
    Instruction *iInst = dyn_cast<Instruction>(&*i);

    LoadInst *LI = dyn_cast<LoadInst>(iInst);
    StoreInst *SI = dyn_cast<StoreInst>(iInst);
    CallInst *CI = dyn_cast<CallInst>(iInst);

    if (LI) {
      MemLoc1 = new MemoryLocation();
      *MemLoc1 = MemoryLocation::get(LI);
      addrSet.insert(std::make_pair(LI->getPointerOperand(),
                                    MemLoc1));
      LPSet.insert(LI->getPointerOperand());
    } else if (SI) {
      MemLoc1 = new MemoryLocation();
      *MemLoc1 = MemoryLocation::get(SI);
      addrSet.insert(std::make_pair(SI->getPointerOperand(),
                                    MemLoc1));
    } else if (CI) {
      ImmutableCallSite CS(CI);
      for (unsigned j = 0, ej = CI->getNumArgOperands();
           j < ej; ++j) {
        Value *Arg = CI->getArgOperand(j);

        if (!Arg->getType()->isPointerTy())
          continue;

        MemLoc1 = new MemoryLocation();
        *MemLoc1 = MemoryLocation::getForArgument(CS, j, *TLI);
        addrSet.insert(std::make_pair(Arg, MemLoc1));

        Function *CF = CI->getCalledFunction();
        if (CF && CF->getName() == "copy_from_user" && j < 2)  
          DstSet.insert(Arg);
      }
    }
  }

  for (std::set<AddrMemPair>::iterator it1 = addrSet.begin(),
       it1e = addrSet.end(); it1 != it1e; ++it1) {
    Addr1 = it1->first;
    MemLoc1 = it1->second;

    for (std::set<AddrMemPair>::iterator it2 = addrSet.begin(),
         it2e = addrSet.end(); it2 != it2e; ++it2) {
      if (it2 == it1)
        continue;

      Addr2 = it2->first;
      MemLoc2 = it2->second;

      if (Addr1 == Addr2)
        continue;

      AliasResult AResult = AAR.alias(*MemLoc1, *MemLoc2);

#ifdef MUST_ALIAS
      if (AResult != MustAlias && AResult != PartialAlias) {
#else
			if (AResult == NoAlias) {
#endif
        bool flag = true;

        if (AResult == MayAlias) {
          CallInst *CI;
          Function *CF;
          StringRef CFName;

          CI = dyn_cast<CallInst>(Addr1);
          if (CI) {
            CF = CI->getCalledFunction();
            if (CF) {
              CFName = CF->getName();
              if (CFName.contains("kmalloc"))  
                flag = false;
            }
          }
          CI = dyn_cast<CallInst>(Addr2);
          if (CI) {
            CF = CI->getCalledFunction();
            if (CF) {
              CFName = CF->getName();
              if (CFName.contains("kmalloc"))
                flag = false;
            }
          }
          //Hack for copy_from_user - dst and SCheck
          if (DstSet.find(Addr1) != DstSet.end() && 
              LPSet.find(Addr2) != LPSet.end())
            flag = false;
          if (DstSet.find(Addr2) != DstSet.end() &&
              LPSet.find(Addr1) != LPSet.end())
            flag = false;
          if (DstSet.find(Addr1) != DstSet.end() && 
            DstSet.find(Addr2) !=  DstSet.end())
            flag = false;
        }

        if (flag)
          continue;
      }

      auto as = aliasPtrs.find(Addr1);
      if (as == aliasPtrs.end()) {
        std::set<Value *> sv;
        sv.clear();
        sv.insert(Addr2);
        aliasPtrs[Addr1] = sv;
      } else if (as->second.count(Addr2) == 0) {
        as->second.insert(Addr2);
      }
    }
  }
}

bool PointerAnalysisPass::doModulePass(Module *M) {
  // Save TargetLibraryInfo.
  Triple ModuleTriple(M->getTargetTriple());
  TargetLibraryInfoImpl TLII(ModuleTriple);
  TLI = new TargetLibraryInfo(TLII);

  // Run BasicAliasAnalysis pass on each function in this module.
  // XXX: more complicated alias analyses may be required.
  legacy::FunctionPassManager *FPasses = new legacy::FunctionPassManager(M);
  AAResultsWrapperPass *AARPass = new AAResultsWrapperPass();

  FPasses->add(AARPass);

  FPasses->doInitialization();
  for (Function &F : *M) {
    if (F.isDeclaration())
      continue;
    FPasses->run(F);
  }
  FPasses->doFinalization();

  // Basic alias analysis result.
  AAResults &AAR = AARPass->getAAResults();

  for (Module::iterator f = M->begin(), fe = M->end();
       f != fe; ++f) {
    Function *F = &*f;
    PointerAnalysisMap aliasPtrs;

    if (F->empty())
      continue;

    detectAliasPointers(F, AAR, aliasPtrs);

    // Save pointer analysis result.
    Ctx->FuncPAResults[F] = aliasPtrs;
    Ctx->FuncAAResults[F] = &AAR;
  }

  return false;
}
