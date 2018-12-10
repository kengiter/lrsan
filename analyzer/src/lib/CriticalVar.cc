#include <llvm/IR/DebugInfo.h>
#include <llvm/Pass.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/Function.h>
#include <llvm/Support/Debug.h>
#include <llvm/IR/InstIterator.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/Value.h>
#include <llvm/IR/CFG.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/InlineAsm.h>
#include <llvm/ADT/StringExtras.h>
#include <llvm/Analysis/CallGraph.h>

#include "CriticalVar.h"
#include "Config.h"

// don't check critical function
#define VARIABLE_ONLY 1

#define ENABLE_SOURCE_TRACKING 0

// Number of addresses to track for each identified
// modification of a critical variable
#define NUM_ADDRESSES_TO_TRACK 50

// Number of functions to track for each address
#define NUM_FUNCTIONS_TO_TRACK 50

#define DETECT_MEMFUNC_MODIFICATIONS true

#define DETECT_STORE_MODIFICATIONS true

//#define DEBUG_PRINT 1


using namespace llvm;

/// Check if the value is a constant.
bool CriticalVarPass::isConstant(Value *V) {
  // Invalid input.
  if (!V)
    return false;

  // The value is a constant.
  Constant *Ct = dyn_cast<Constant>(V);
  if (Ct)
    return true;

  return false;
}

/// Check if the value is an errno.
bool CriticalVarPass::isValueErrno(Value *V) {
  // Invalid input.
  if (!V)
    return false;

  // The value is a constant integer.
  ConstantInt *CI = dyn_cast<ConstantInt>(V);
  if (CI && (CI->getType()->getBitWidth() == 32 ||
      CI->getType()->getBitWidth() == 64)) {
    const APInt &value = CI->getValue();
    // The value is an errno (negative or positive).
    if (is_errno(-value) || is_errno(value))
      return true;
  }

  // The value is a constant expression.
  ConstantExpr *CE = dyn_cast<ConstantExpr>(V);
  if (CE) {
    for (unsigned i = 0, e = CE->getNumOperands();
         i != e; ++i) {
      if (isValueErrno(CE->getOperand(i)))
        return true;
    }
  }

  return false;
}

/// Check if the value is a parameter of the function.
/// The current analysis is limited to a scope of a function.
bool CriticalVarPass::isFunctionParameter(Value *V, Function *F)
{
  for (Function::arg_iterator arg = F->arg_begin();
       arg != F->arg_end(); arg++) {
    if (arg == V)
      return true;
  }
  return false;
}

/// Dump the marked CFG edges.
void CriticalVarPass::dumpEdges(EdgeErrnoFlag &errnoEdges) {
  for (auto it = errnoEdges.begin(); it != errnoEdges.end(); ++it) {
    CFGEdge edge = it->first;
    ErrnoFlag flag = it->second;
    TerminatorInst *TI = edge.first;

    if (NULL == TI) {
      OP << "    An errno ";
      OP << (flag == Must_Return_Errno ? "must" : "may");
      OP << " be returned.\n";
      continue;
    }

    BasicBlock *predBB = TI->getParent();
    BasicBlock *succBB = TI->getSuccessor(edge.second);

    OP << "    ";
    predBB->printAsOperand(OP, false);
    OP << " -> ";
    succBB->printAsOperand(OP, false);
    OP << ": ";
    OP << (flag == Must_Return_Errno ? "Must" : "May");
    OP << '\n';
  }
}

/// Make the given edge with errno flag.
void CriticalVarPass::markEdge(CFGEdge &CE, ErrnoFlag flag, EdgeErrnoFlag &errnoEdges) {
	errnoEdges[CE] = flag;
}

/// Mark the edge from predBB to succBB with errno flag.
void CriticalVarPass::markEdgeBlockToBlock(BasicBlock *predBB, BasicBlock *succBB,
                                           ErrnoFlag flag, EdgeErrnoFlag &errnoEdges) {
  TerminatorInst *TI;

  // Invalid inputs.
  if (!predBB || !succBB)
    return;

  TI = predBB->getTerminator();
  for (unsigned i = 0, e = TI->getNumSuccessors(); i != e; ++i) {
    BasicBlock *BB = TI->getSuccessor(i);
    if (BB != succBB)
      continue;

    // Mark the edge
    CFGEdge CE = std::make_pair(TI, i);
    markEdge(CE, flag, errnoEdges);
  }
}

/// Recursively mark all edges from this BB with errno flag.
void CriticalVarPass::markAllEdgesFromBlock(BasicBlock *BB, ErrnoFlag flag,
                                            EdgeErrnoFlag &errnoEdges) {
  TerminatorInst *TI;

  // Invalid input.
  if (!BB)
    return;

	std::set<BasicBlock *> PB;
  std::list<BasicBlock *> EB;
  PB.clear();
  EB.clear();

  EB.push_back(BB);
	while (!EB.empty()) {
		
		BasicBlock *TB = EB.front();
    EB.pop_front();
    if (PB.count(TB) != 0)
      continue;
    PB.insert(TB);

		// Iterate on each successor basic block.
		TI = TB->getTerminator();
		for (unsigned i = 0, e = TI->getNumSuccessors(); i != e; ++i) {
			CFGEdge CE = std::make_pair(TI, i);
			// FIXME: this is a backward tracking, so do not overwrite
			// previously marked flags
			if (errnoEdges.count(CE) > 0 && flag != errnoEdges[CE])
				markEdge(CE, May_Return_Errno, errnoEdges);
			else
				markEdge(CE, flag, errnoEdges);


			EB.push_back(TI->getSuccessor(i));
		}
  }
}

/// Mark all edges to this BB with errno flag.
void CriticalVarPass::markAllEdgesToBlock(BasicBlock *BB, ErrnoFlag flag,
                                          EdgeErrnoFlag &errnoEdges) {
  // Invalid input.
  if (!BB)
    return;

	std::set<BasicBlock *> PB;
  std::list<BasicBlock *> EB;
  PB.clear();
  EB.clear();

  EB.push_back(BB);
	while (!EB.empty()) {
		
		BasicBlock *TB = EB.front();
    EB.pop_front();
    if (PB.count(TB) != 0)
      continue;
    PB.insert(TB);
		// Iterate on each predecessor basic block.
		for (BasicBlock *predBB : predecessors(TB)) {
			TerminatorInst *TI = predBB->getTerminator();

			for (unsigned i = 0, e = TI->getNumSuccessors(); i != e; ++i) {
				if (TB != TI->getSuccessor(i))
					continue;

				CFGEdge CE = std::make_pair(TI, i);
				if (errnoEdges.count(CE) > 0 && flag != errnoEdges[CE])
					markEdge(CE, May_Return_Errno, errnoEdges);
				else
					markEdge(CE, flag, errnoEdges);
				break;
			}
			EB.push_back(predBB);
		}
	}
}

/// Add edges to analyze for a phinode.
void CriticalVarPass::addPHINodeEdges(PHINode *PN, std::list<EdgeValue> &EV) {
  BasicBlock *BB = PN->getParent();

  for (unsigned i = 0, e = PN->getNumIncomingValues(); i != e; ++i) {
    Value *IV = PN->getIncomingValue(i);
    BasicBlock *inBB = PN->getIncomingBlock(i);
    TerminatorInst *TI = inBB->getTerminator();

    for (unsigned j = 0, f = TI->getNumSuccessors(); j != f; ++j) {
      if (BB != TI->getSuccessor(j))
        continue;
      EV.push_back(std::make_pair(std::make_pair(TI, j), IV));
    }
  }
}

/// Recursively do errno check for this value.
void CriticalVarPass::recurCheckValueErrno(Function *F, std::list<EdgeValue> &EV,
                                           std::set<Value *> &PV, 
                                           EdgeErrnoFlag &errnoEdges) {
  EdgeValue E = EV.front();
  CFGEdge CE = E.first;
  Value *V = E.second;

  EV.pop_front();

  // Invalid input.
  if (!V)
    return;

  // Check if the value was checked before.
  if (PV.count(V) != 0)
    return;
  PV.insert(V);

  // The value is a load. Let's find out the previous stores.
  if (auto LI = dyn_cast<LoadInst>(V)) {
    Value *LPO = LI->getPointerOperand();

    list<pair<CFGEdge, BasicBlock *>> RList;
    set<BasicBlock *> PSet;

    RList.clear();
    PSet.clear();
    RList.push_back(make_pair(CE, LI->getParent()));

    while (!RList.empty()) {
      auto REL = RList.front();
      RList.pop_front();

      CFGEdge RCE = REL.first;
      BasicBlock *BB = REL.second;
      if (PSet.find(BB) != PSet.end())
        continue;
      PSet.insert(BB);

      BasicBlock::reverse_iterator rit;
      if (BB == LI->getParent())
        rit = ++LI->getIterator().getReverse();
      else
        rit = BB->rbegin();
      BasicBlock::reverse_iterator rie = BB->rend();

      bool found = false;
      for (; rit != rie; ++rit) {
        auto SI = dyn_cast<StoreInst>(&*rit);
        if (SI && SI->getPointerOperand() == LPO) {
          Value *SVO = SI->getValueOperand();
          if (isConstant(SVO)) {
            if (isValueErrno(SVO)) {
							markAllEdgesFromBlock(BB, Must_Return_Errno, errnoEdges);
							markAllEdgesToBlock(BB, Must_Return_Errno, errnoEdges);
						}
						else {
							// likely not error code
							markAllEdgesFromBlock(BB, May_Return_Errno, errnoEdges);
							markAllEdgesToBlock(BB, May_Return_Errno, errnoEdges);
						}
          } else {
            EV.push_back(make_pair(RCE, SVO));
          }
          found = true;
          break;
        }
      }
      if (found)
        continue;

      pred_iterator pi = pred_begin(BB), pe = pred_end(BB);
      for (; pi != pe; ++pi) {
        BasicBlock *PBB = *pi;
        TerminatorInst *TI = PBB->getTerminator();
        unsigned num = TI->getNumSuccessors();
        for (unsigned i = 0; i < num; ++i) {
          if (BB != TI->getSuccessor(i))
            continue;
          RList.push_back(make_pair(make_pair(TI, i), PBB));
        }
      }
    }

    return;
  }

  // The value is a phinode.
  PHINode *PN = dyn_cast<PHINode>(V);
  if (PN) {
    BasicBlock *BB = PN->getParent();

    // Check each incoming value.
    for (unsigned i = 0, e = PN->getNumIncomingValues(); i != e; ++i) {
      Value *IV = PN->getIncomingValue(i);
      BasicBlock *inBB = PN->getIncomingBlock(i);

      // The incoming value is a constant.
      if (isConstant(IV)) {
        if (isValueErrno(IV)) {
          markEdgeBlockToBlock(inBB, BB, Must_Return_Errno, errnoEdges);
        }
      } else {
        // Add the incoming value and the corresponding edge to the list.
        TerminatorInst *TI = inBB->getTerminator();
        for (unsigned j = 0, f = TI->getNumSuccessors(); j != f; ++j) {
          if (BB != TI->getSuccessor(j))
            continue;
          EV.push_back(std::make_pair(std::make_pair(TI, j), IV));
        }
      }
    }

    return;
  }

  // The value is a select instruction.
  SelectInst *SI = dyn_cast<SelectInst>(V);
  if (SI) {
    Value *Cond = SI->getCondition();
    Value *SV;
    bool flag1 = false;
    bool flag2 = false;

    SV = SI->getTrueValue();
    if (isConstant(SV)) {
      if(isValueErrno(SV))
        flag1 = true;
    } else {
      EV.push_back(std::make_pair(CE, SV));
    }

    SV = SI->getFalseValue();
    if (isConstant(SV)) {
      if (isValueErrno(SV))
        flag2 = true;
    } else {
      EV.push_back(std::make_pair(CE, SV));
    }

    if (flag1 && flag2) {
      markEdge(CE, Must_Return_Errno, errnoEdges);
		}
    else if (flag1 || flag2) {
      markEdge(CE, May_Return_Errno, errnoEdges);
		}
    return;
  }

  // The value is a getelementptr instruction.
  GetElementPtrInst *GEP = dyn_cast<GetElementPtrInst>(V);
  if (GEP) {
    Value *PO = GEP->getPointerOperand();

    if (isConstant(PO)) {
      if (isValueErrno(PO)) {
        markEdge(CE, Must_Return_Errno, errnoEdges);
			}
    } else
      EV.push_back(std::make_pair(CE, PO));
    return;
  }

  // The value is an unary instruction.
  UnaryInstruction *UI = dyn_cast<UnaryInstruction>(V);
  if (UI) {
    Value *UO = UI->getOperand(0);
    if (isConstant(UO)) {
      if (isValueErrno(UO)) {
        markEdge(CE, Must_Return_Errno, errnoEdges);
			}
    }
    else
      EV.push_back(std::make_pair(CE, UO));
    return;
  }

  // The value is a binary operator.
  BinaryOperator *BO = dyn_cast<BinaryOperator>(V);
  if (BO) {
    for (unsigned i = 0, e = BO->getNumOperands();
         i != e; ++i) {
      Value *Opd = BO->getOperand(i);
      if (isConstant(Opd)) {
        if (isValueErrno(Opd)) {
          markEdge(CE, May_Return_Errno, errnoEdges);
          return;
        }
      } else
        EV.push_back(std::make_pair(CE, Opd));
    }
    return;
  }

  // The value is an insertvalue instruction.
  InsertValueInst *IV = dyn_cast<InsertValueInst>(V);
  if (IV) {
    Value *AO = IV->getAggregateOperand();
    Value *IVO = IV->getInsertedValueOperand();
    bool flag1 = false;
    bool flag2 = false;

    if (isConstant(AO)) {
      if (isValueErrno(AO))
        flag1 = true;
    } else
      EV.push_back(std::make_pair(CE, AO));

    if (isConstant(IVO)) {
      if (isValueErrno(IVO))
        flag2 = true;
    } else
      EV.push_back(std::make_pair(CE, IVO));

    if (flag1 && flag2) {
      markEdge(CE, Must_Return_Errno, errnoEdges);
		}
    else if (flag1 || flag2) {
      markEdge(CE, May_Return_Errno, errnoEdges);
		}

    return;
  }

  // The value is a icmp instruction. Skip it.
  ICmpInst *ICI = dyn_cast<ICmpInst>(V);
  if (ICI)
    return;

  // The value is a call instruction. 
  CallInst *CaI = dyn_cast<CallInst>(V);
  if (CaI) {
		// FIXME: should be handled
		return;
	}

  // The value is a parameter of the fucntion. Skip it.
  if (isFunctionParameter(V, F))
    return;

  // TODO: support more LLVM IR types.
#ifdef DEBUG_PRINT
  OP << "== Warning: unsupported LLVM IR:"
     << *V << '\n';
#endif
}



/// Check if the value must or may be an errno.
void CriticalVarPass::checkValueErrno(Function *F, Value *V, EdgeErrnoFlag &errnoEdges) {
  std::list<EdgeValue> EV;
  std::set<Value *> PV;

  // Invalid input.
  if (!V)
    return;

  EV.clear();
  PV.clear();

  // Construct the list.
  EV.push_back(std::make_pair(std::make_pair((TerminatorInst *)NULL, 0), V));

#ifdef DEBUG_PRINT
  if (EV.empty())
    OP << "== Warning: unsupported LLVM IR:"
       << *V << '\n';
#endif

  // Iterate each value in the list.
  while (!EV.empty()) {
    // Recursively do errno check for this value.
    recurCheckValueErrno(F, EV, PV, errnoEdges);
  }
}

/// Check if the returned value must be or may be an errno.
/// And mark the traversed edges in the CFG.
void CriticalVarPass::checkReturnValue(Function *F, BasicBlock *BB, Value *V,
                                       EdgeErrnoFlag &errnoEdges) {
  std::list<std::pair<CFGEdge, Value *>> EV;

  // Invalid input.
  if (!V)
    return;
  
  // The returned value is not a constant. Further check is required.
  checkValueErrno(F, V, errnoEdges);
}

/// Traverse the CFG and find security checks.
void CriticalVarPass::findSecurityChecks(Function *F,
                                         EdgeErrnoFlag &errnoEdges,
                                         std::set<Value *> &securityChecks) {
  // Find blocks that contain security checks.
  for (Function::iterator b = F->begin(), e = F->end();
       b != e; ++b) {
    BasicBlock *BB = &*b;
    TerminatorInst *TI = BB->getTerminator();
    bool SecureCond1 = false; /* One path must return an errno. */
    bool SecureCond2 = false; /* One path must have possibility not to return any errno. */

		int NumSucc = TI->getNumSuccessors();
    if (NumSucc < 2)
      continue;
		
		// Query marked edge graph to determine security checks
		int errFlag = 0;
		int NumMayErr = 0, NumMustErr = 0;
		for (unsigned i = 0; i < NumSucc; ++i) {
			errFlag = errnoEdges[std::make_pair(TI, i)];
			if (errFlag == Must_Return_Errno)
				++NumMustErr;
			else
				++NumMayErr;
		}

		// not a security check
		if (!(NumMustErr && NumMayErr)) {
			continue;
		}


    // This BB has security check, try to find it.
    Value *Cond = NULL;
    BranchInst *BI;
    SwitchInst *SI;

    BI = dyn_cast<BranchInst>(TI);
    if (BI)
      Cond = BI->getCondition();
    else {
      SI = dyn_cast<SwitchInst>(TI);
      if (SI)
        Cond = SI->getCondition();
    }

    if (!Cond)
      OP << "Warnning: cannot find the check.\n";
    else {
      //OP << "Condition of the check is: " << *Cond << "\n";
      securityChecks.insert(Cond);
    }
  }
}

/// Backtrack to find critical variables/functions.
void CriticalVarPass::findCriticalVariable(Function *F, Value *SCheck, Value *V,
                                  std::map<Value *, std::set<Value *>> &CheckToVars,
                                          std::set<Value *> &PSet) {
  if (isConstant(V))
    return;

  if (PSet.count(V) != 0)
    return;
  PSet.insert(V);

  if (isFunctionParameter(V, F)) {
#ifdef DEBUG_PRINT
    OP << "== A critical variable is identified (function parameter):\n";
    OP << "\033[31m" << *V << "\033[0m" << "\n";
#endif
    CheckToVars[SCheck].insert(V);
    return;
  }

  CallInst *CI = dyn_cast<CallInst>(V);
  if (CI) {
#ifndef VARIABLE_ONLY
    CheckToVars[SCheck].insert(CI);
#endif
    return;
  }

  LoadInst *LI = dyn_cast<LoadInst>(V);
  if (LI) {
    // XXX: we may need to save the memory address for future analysis.
    // XXX: we may also need to save the type of the loaded value, pointer or non-pointer.
    CheckToVars[SCheck].insert(LI);
    return;
  }

  SelectInst *SI = dyn_cast<SelectInst>(V);
  if (SI) {
    findCriticalVariable(F, SCheck, SI->getTrueValue(), CheckToVars, PSet);
    findCriticalVariable(F, SCheck, SI->getFalseValue(), CheckToVars, PSet);
    return;
  }

  GetElementPtrInst *GEP = dyn_cast<GetElementPtrInst>(V);
  if (GEP)
    return findCriticalVariable(F, SCheck, GEP->getPointerOperand(), CheckToVars, PSet);

  PHINode *PN = dyn_cast<PHINode>(V);
  if (PN) {
    for (unsigned i = 0, e = PN->getNumIncomingValues(); i != e; ++i) {
      Value *IV = PN->getIncomingValue(i);

      findCriticalVariable(F, SCheck, IV, CheckToVars, PSet);
    }
    return;
  }

  ICmpInst *ICmp = dyn_cast<ICmpInst>(V);
  if (ICmp) {
    for (unsigned i = 0; i < ICmp->getNumOperands(); ++i) {
      Value *Opd = ICmp->getOperand(i);

      findCriticalVariable(F, SCheck, Opd, CheckToVars, PSet);
    }
    return;
  }

  BinaryOperator *BO = dyn_cast<BinaryOperator>(V);
  if (BO) {
    for (unsigned i = 0, e = BO->getNumOperands();
         i != e; ++i) {
      Value *Opd = BO->getOperand(i);
      if (isConstant(Opd))
        continue;

      findCriticalVariable(F, SCheck, Opd, CheckToVars, PSet);
    }
    return;
  }

  UnaryInstruction *UI = dyn_cast<UnaryInstruction>(V);
  if (UI) {
    Value *UO = UI->getOperand(0);
    findCriticalVariable(F, SCheck, UO, CheckToVars, PSet);
    return;
  }

  // TODO: more analyses are required.
#ifdef DEBUG_PRINT
  OP << "== Warning: unsupported LLVM IR when find critical variables/functions:";
  OP << "\033[31m" << *V << "\033[0m\n";
#endif
}

/// Identify critical variables/functions used in each security check.
void CriticalVarPass::identifyCriticalVariables(Function *F,
                                        std::set<Value *> &securityChecks,
                               std::map<Value *, std::set<Value *>> &CheckToVars) {
  for (std::set<Value *>::iterator it = securityChecks.begin();
       it != securityChecks.end(); ++it) {
    Value *Cond = *it;
    std::set<Value *> PSet;

    PSet.clear();

    CheckToVars[Cond].clear();

    findCriticalVariable(F, Cond, Cond, CheckToVars, PSet);
  }
}

void CriticalVarPass::filterNonGlobalVariables(
		std::map<Value *, std::set<Value *>> &CheckToVars) {
  for (std::map<Value *, 
			std::set<Value *>>::iterator it = CheckToVars.begin();
      it != CheckToVars.end(); ++it) {
    
		std::set<Value *> cvSet = it->second;

    for (std::set<Value *>::iterator cit = cvSet.begin(); cit != cvSet.end(); cit++) {
      Value *CV = *cit;

      if (dyn_cast<GlobalVariable>(CV)) {
        cvSet.erase(cit);
        continue;
      }

      LoadInst *LI = dyn_cast<LoadInst>(CV);
      if (LI) {
        Value *Addr = LI->getPointerOperand();

        if (dyn_cast<GlobalVariable>(Addr)) {
          cvSet.erase(cit);
          continue;
        }
      }
    }
  }
}

/// Check alias result of two values.
/// True: alias, False: not alias.
bool CriticalVarPass::checkAlias(Value *Addr1, Value *Addr2,
                            PointerAnalysisMap &aliasPtrs) {

  if (Addr1 == Addr2)
    return true;

  auto it = aliasPtrs.find(Addr1);
  if (it != aliasPtrs.end()) {
    if (it->second.count(Addr2) != 0)
      return true;
  }

  // May not need to do this further check.
  it = aliasPtrs.find(Addr2);
  if (it != aliasPtrs.end()) {
    if (it->second.count(Addr1) != 0)
      return true;
  }

  return false;
}

/// Get aliased pointers for this pointer.
void CriticalVarPass::getAliasPointers(Value *Addr,
                           std::set<Value *> &aliasAddr,
                          PointerAnalysisMap &aliasPtrs) {
  aliasAddr.clear();
  aliasAddr.insert(Addr);

  auto it = aliasPtrs.find(Addr);
  if (it == aliasPtrs.end())
    return;

  for (auto itt : it->second)
    aliasAddr.insert(itt);
}

/// Get called function name of V.
StringRef CriticalVarPass::getCalledFuncName(Value *V) {
  assert(V);

  InlineAsm *IA = dyn_cast<InlineAsm>(V);
  if (IA)
    return StringRef(IA->getAsmString());

  User *UV = dyn_cast<User>(V);
  if (UV) {
    if (UV->getNumOperands() == 0)
      return UV->getName();
    Value *VUV = UV->getOperand(0);
    return VUV->getName();
  }

  return V->getName();
}

/// Check if it is possible for an instruction Inst to use
/// the value stored by the store instruction St.
/// If yes, there exists a path on the CFG from St to Inst.
/// Return value: true -> possible, false -> impossible.
bool CriticalVarPass::possibleUseStResult(Instruction *Inst,
                                          Instruction *St) {
  BasicBlock *InstBB = Inst->getParent();
  BasicBlock *StBB = St->getParent();
  std::set<BasicBlock *> PB;
  std::list<BasicBlock *> EB;

  if (!InstBB || !StBB)
    return false;

  // Inst and St are in the same basic block.
  // Iterate over each instruction see which
  // one is seen firstly.
  if (InstBB == StBB) {
    for (BasicBlock::iterator I = InstBB->begin(),
         IE = InstBB->end(); I != IE; ++I) {
        Instruction *II = dyn_cast<Instruction>(&*I);
        if (II == St)
          return true;
        if (II == Inst)
          return false;
    }
    return false;
  }

  // See if there exists a path from StBB to LdBB.
  PB.clear();
  EB.clear();

  EB.push_back(StBB);
  while (!EB.empty()) {
    BasicBlock *TB = EB.front();
    TerminatorInst *TI;

    EB.pop_front();
    if (PB.count(TB) != 0)
      continue;
    PB.insert(TB);

    if (TB == InstBB)
      return true;

    TI = TB->getTerminator();
    for (unsigned i = 0, e = TI->getNumSuccessors();
         i != e; ++i) {
      EB.push_back(TI->getSuccessor(i));
    }
  }

  return false;
}

void CriticalVarPass::trackAliasAddrAfterModification(Value *Addr,
                                          std::set<Value *> &pSet,
                                              Instruction *StInst,
                                  std::set<InstructionUseDef> &udSet,
                                    PointerAnalysisMap &aliasPtrs) {
  std::set<Value *> aliasAddr;

  getAliasPointers(Addr, aliasAddr, aliasPtrs);
  for (auto it : aliasAddr) {
    for (User *UA : it->users()) {
      Instruction *UAInst = dyn_cast<Instruction>(UA);

      if (!UAInst)
        continue;

      // Check if it is possible for this instruciton to
      // use the value stored by the store instruction.
      // It is possible that the instruction may happen before
      // the store instruction on the CFG.
      if (!possibleUseStResult(UAInst, StInst))
        continue;

      LoadInst *UALoad = dyn_cast<LoadInst>(UA);
      if (UALoad) {
        findUseDefOfCriticalVar(dyn_cast<Value>(UA),
                                UALoad->getPointerOperand(),
                                pSet, udSet, aliasPtrs);
        continue;
      }

#if DETECT_MEMFUNC_MODIFICATIONS
      CallInst *UACall = dyn_cast<CallInst>(UA);
      if (UACall) {
        Value *CaV = UACall->getCalledValue();
        StringRef FuncName = getCalledFuncName(CaV);
        auto FIter = Ctx->MemWriteFuncs.find(FuncName.str());
        if (FIter == Ctx->MemWriteFuncs.end())
          continue;
        if (UACall->getNumArgOperands() < FIter->second)
          continue;
        if (UACall->getArgOperand(FIter->second) != it)
          continue;
        udSet.insert(std::make_pair(UACall, Defined));
      }
#endif
    }
  }
}

/// Check whether a value is exploitable. Intuitively,
/// the following cases are unexploitable and will be
/// filtered out:
/// 1. The operation result of the orignal value of a
///      critical variable and a constant.
/// 2. The return value of a function call (???).
/// 3. More cases to be added here ...
///
/// Return value:
/// false: the value is unexploitable.
/// true: the value may be exploitable.
bool CriticalVarPass::isExploitable(Value *V, Value *Addr,
                          std::map<Value *, bool> &Expt,
                                 std::set<Value *> &PV) {
  if (Expt.count(V) > 0)
    return Expt[V];

  // This value was checked before.
  if (PV.count(V) > 0) {
    Expt[V] = true;
    return Expt[V];
  }
  PV.insert(V);

  if (isConstant(V)) {
    Expt[V] = false;
    return Expt[V];
  }

  CallInst *CI = dyn_cast<CallInst>(V);
  if (CI) {
    Expt[V] = false;
    return Expt[V];
  }

  LoadInst *LI = dyn_cast<LoadInst>(V);
  if (LI) {
    if (LI->getPointerOperand() == Addr)
      Expt[V] = false;
    else
      Expt[V] = true;
    return Expt[V];
  }

  UnaryInstruction *UI = dyn_cast<UnaryInstruction>(V);
  if (UI) {
    Expt[V] = isExploitable(UI->getOperand(0), Addr,
                            Expt, PV);
    return Expt[V];
  }

  PHINode *PN = dyn_cast<PHINode>(V);
  if (PN) {
    bool flag = false;
    for (unsigned i = 0, e = PN->getNumIncomingValues();
         i != e; ++i) {
      if (isExploitable(PN->getIncomingValue(i),
                               Addr, Expt, PV)) {
        flag = true;
        break;
      }
    }
    Expt[V] = flag;
    return Expt[V];
  }

  BinaryOperator *BO = dyn_cast<BinaryOperator>(V);
  if (BO) {
    if (BO->getNumOperands() != 2)
      Expt[V] = true;
    else {
      Value *Opd0 = BO->getOperand(0);
      Value *Opd1 = BO->getOperand(1);
      Expt[V] = (isExploitable(Opd0, Addr, Expt, PV) ||
                isExploitable(Opd1, Addr, Expt, PV));
    }
    return Expt[V];
  }

  Expt[V] = true;

  return Expt[V];
}

/// Find modifications by function calls.
void CriticalVarPass::findUseDefFromCVAddr(Function *F,
                                         Value *cvAddr,
                               std::set<Value *> &pSet,
                    std::set<InstructionUseDef> &udSet,
                         PointerAnalysisMap &aliasPtrs) {
  for (inst_iterator i = inst_begin(F), ei = inst_end(F);
       i != ei; ++i) {
    Instruction *Inst = dyn_cast<Instruction>(&*i);
    StringRef FName;
    Value *CaV;
    unsigned ej;

    StoreInst *SI = dyn_cast<StoreInst>(Inst);
    if (SI) {
      CaV = SI->getPointerOperand();
      if (!checkAlias(cvAddr, CaV, aliasPtrs))
        continue;

      // This store may modify a critical variable.
      #if DETECT_STORE_MODIFICATIONS
      std::map<Value *, bool> Expt;
      std::set<Value *> PV;
      Expt.clear();
      PV.clear();
      if (isExploitable(SI->getValueOperand(),
                        SI->getPointerOperand(),
                        Expt, PV))
        udSet.insert(std::make_pair(Inst, Defined));
      #endif

      // We need to continue to track address aliased with this address.
      trackAliasAddrAfterModification(CaV, pSet, Inst, udSet, aliasPtrs);

      // We also need to track the value stored.
      // TODO: further analysis is required to find the source of this value.
      findUseDefOfCriticalVar(SI->getValueOperand(), NULL, pSet, udSet, aliasPtrs);

      continue;
    }

    CallInst *CI = dyn_cast<CallInst>(Inst);
    InvokeInst *II = dyn_cast<InvokeInst>(Inst);
    if (CI) {
      CaV = CI->getCalledValue();
      ej = CI->getNumArgOperands();
    } else if (II) {
      CaV = II->getCalledValue();
      ej = II->getNumArgOperands();
    } else
      continue;

    FName = getCalledFuncName(CaV);
    if (FName.contains("llvm") && FName.contains("lifetime"))
      continue;

    for (unsigned j = 0; j != ej; ++j) {
      Value *Arg;

      if (CI)
        Arg = CI->getArgOperand(j);
      else
        Arg = II->getArgOperand(j);

      if (!checkAlias(cvAddr, Arg, aliasPtrs))
        continue;

#if DETECT_MEMFUNC_MODIFICATIONS
			StringRef FuncName = getCalledFuncName(CaV);
			auto FIter = Ctx->MemWriteFuncs.find(FuncName.str());
			if (FIter != Ctx->MemWriteFuncs.end() && FIter->second == j) {
        udSet.insert(std::make_pair(Inst, Defined));
        break;
      }
#endif

      // We need to continue to track address aliased with this address.
      trackAliasAddrAfterModification(Arg, pSet, Inst, udSet, aliasPtrs);

      break;
    }
  }
}

/// Find all instructions that use the critical variable.
/// Note: a use can be a direct use or an indirect use.
/// XXX: the definition of critical use:
/// 1. A critical variable is used as a memory address 
///    in memory access instructions, e.g., load and store.
/// 2. A critical variable is used as a parameter of a
///    critical function, e.g., kmalloc, memcpy, and memset.
/// 3. A critical variable is used as a function pointer.
/// 4. A critical varialbe is used as an index in GEP.
/// 5. More???
/// XXX: the definition of modification:
/// 1. An updated value of a critical variable is stored 
///    to the address that is aliased with the address of
///    the critical variable. That is a store instruction
///    is required for this type of modification.
/// 2. Stack variable???
void CriticalVarPass::findUseDefOfCriticalVar(Value *CV,
                                           Value *cvAddr,
                                     std::set<Value *> &pSet,
                           std::set<InstructionUseDef> &udSet,
                               PointerAnalysisMap &aliasPtrs) {
  if (!CV || isConstant(CV))
    return;

  if (pSet.count(CV) != 0)
    return;
  pSet.insert(CV);

  for (User *U : CV->users()) {
    if (Instruction *Inst = dyn_cast<Instruction>(U)) {
      /* For the following types of instructions, we continue to 
         find their uses. */
      // XXX: typical use instructions:  BinaryOperator, CallInst,
      // ExtractElementInst, GetElementInst, InsertElementInst,
      // InsertValueInst, PHINode, SelectInst, ShuffleVectorInst,
      // StoreInst, ReturnInst, CastInst, ExtractValueInst, LoadInst.
      //
      // XXX: typical define instructions: StoreInst, 
      // BinaryOperator, CallInst, InvokeInst, InsertElementInst,
      // InsertValueInst, ShuffleVectorInst, and CastInst.

      SelectInst *SI = dyn_cast<SelectInst>(Inst);
      if (SI) {
        Value *V = dyn_cast<Value>(SI);
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      PHINode *PN = dyn_cast<PHINode>(Inst);
      if (PN) {
        Value *V = dyn_cast<Value>(PN);
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      CastInst *CI = dyn_cast<CastInst>(Inst);
      if (CI) {
        Value *V = dyn_cast<Value>(CI);
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      GetElementPtrInst *GEP = dyn_cast<GetElementPtrInst>(Inst);
      if (GEP) {
        // If used as an index, mark it as a critical use.
        for (Value *Index : make_range(GEP->idx_begin(),
                                       GEP->idx_end())) {
          if (Index == CV) {
            //udSet.insert(std::make_pair(Inst, Used));
            break;
          }
        }
        continue;
      }

      BinaryOperator *BO = dyn_cast<BinaryOperator>(Inst);
      if (BO) {
        // The result of BO may be used for modification,
        // we should continue to track.
        Value *V = dyn_cast<Value>(BO);
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      LoadInst *LI = dyn_cast<LoadInst>(Inst);
      if (LI) {
				// used for memory read, a critical use
        udSet.insert(std::make_pair(Inst, Used));
        continue;
      }

      StoreInst *STI = dyn_cast<StoreInst>(Inst);
      if (STI) {
        Value *Addr = STI->getPointerOperand();
        if (Addr == CV) {
          // Used as an address to write, a critical use ???
          udSet.insert(std::make_pair(Inst, Used));
          continue;
        }

        // Used as a value to write to memory address Addr.
        // If Addr is aliased to the memory address of CV,
        // we need mark this store as a modification.
#if DETECT_STORE_MODIFICATIONS
        if (cvAddr) {
          std::map<Value *, bool> Expt;
          std::set<Value *> PV;
          Expt.clear();
          PV.clear();
          if (checkAlias(Addr, cvAddr, aliasPtrs) &&
              isExploitable(STI->getValueOperand(),
                            Addr, Expt, PV))
            udSet.insert(std::make_pair(Inst, Defined));
        }
#endif

        // Try to find all addresses that are aliased with Addr
        // and continue to track all values loaded from the
        // aliased addresses.
        trackAliasAddrAfterModification(Addr, pSet, Inst,
                                        udSet, aliasPtrs);
        continue;
      }

      ExtractElementInst *EEI = dyn_cast<ExtractElementInst>(Inst);
      if (EEI) {
        Value *V = dyn_cast<Value>(EEI);
				// TODO
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      InsertElementInst *IEI = dyn_cast<InsertElementInst>(Inst);
      if (IEI) {
        Value *V = dyn_cast<Value>(IEI);
				// TODO
        findUseDefOfCriticalVar(V, cvAddr, pSet, udSet, aliasPtrs);
        continue;
      }

      InsertValueInst *IVI = dyn_cast<InsertValueInst>(Inst);
      if (IVI) {
				// TODO
        //Value *V = dyn_cast<Value>(IVI);
        //findUseDefOfCriticalVar(V, pSet, udSet);
        continue;
      }

      ShuffleVectorInst *SVI = dyn_cast<ShuffleVectorInst>(Inst);
      if (SVI) {
				// TODO
        //Value *V = dyn_cast<Value>(SVI);
        //findUseDefOfCriticalVar(V, pSet, udSet);
        continue;
      }

      ExtractValueInst *EVI = dyn_cast<ExtractValueInst>(Inst);
      if (EVI) {
				// TODO
        //Value *V = dyn_cast<Value>(EVI);
        //findUseDefOfCriticalVar(V, pSet, udSet);
        continue;
      }

      ReturnInst *RI = dyn_cast<ReturnInst>(Inst);
      if (RI) {
				// TODO: continue tracking caller?
        continue;
      }

      CallInst *CLI = dyn_cast<CallInst>(Inst);
      InvokeInst *II = dyn_cast<InvokeInst>(Inst);
      if (CLI || II) {
        Value *CaV;
        unsigned numArgCall;
				Value *Arg0;

        if (CLI) {
          CaV = CLI->getCalledValue();
          numArgCall = CLI->getNumArgOperands();
          if (numArgCall > 0)
            Arg0 = CLI->getArgOperand(0);
        } else {
          CaV = II->getCalledValue();
          numArgCall = II->getNumArgOperands();
          if (numArgCall > 0)
            Arg0 = II->getArgOperand(0);
        }

        assert(CaV);
        if (CaV == CV) {
          // Used for indirect call, a critical use
          udSet.insert(std::make_pair(Inst, Used));
        } else {
          // Used as a parameter of the function.
          // We need to further check if this function
          // is a predefined critical function.
          StringRef FName = getCalledFuncName(CaV);

#if DETECT_MEMFUNC_MODIFICATIONS
          auto FIter = Ctx->CriticalFuncs.find(FName.str());
          if (FIter != Ctx->CriticalFuncs.end()) {
            udSet.insert(std::make_pair(Inst, Used));
          } else {
						if ((FName.str().compare("get_user") == 0 
								|| FName.str().compare("__get_user") == 0)
								&& Arg0 == CV) {
							udSet.insert(std::make_pair(Inst, Defined));
						}
						else {
							auto Filter = Ctx->MemWriteFuncs.find(FName.str());
							if (Filter != Ctx->MemWriteFuncs.end())
							udSet.insert(std::make_pair(Inst, Used));
						}
          }
#endif

          // Continue to track CV in callee function.
          // Inter-procedural analysis.
          Function *Callee = AllFuncs[FName];
          Argument *ArgCV;
          unsigned argNo;

          if (!Callee)
            continue;

          for (argNo = 0; argNo < numArgCall; argNo++) {
            Value *Arg;
            if (CLI)
              Arg = CLI->getArgOperand(argNo);
            else
              Arg = II->getArgOperand(argNo);
            if (Arg == CV)
              break;
          }

          for (Function::arg_iterator AI = Callee->arg_begin(),
               E = Callee->arg_end(); AI != E; ++AI) {
            ArgCV = &*AI;
            if (ArgCV->getArgNo() == argNo)
              break;
          }

          findUseDefOfCriticalVar(dyn_cast<Value>(ArgCV), NULL,
                          pSet, udSet, Ctx->FuncPAResults[Callee]);
				}
        continue;
      }

      // TODO: more analyses are required.
			//OP << "== Warning: unsupported LLVM IR when handling uses/defines";
			//OP << *Inst << "\n";
    }
  }
}

static void printSourceCode(unsigned lineno, std::string &src_file) {
  std::ifstream sourcefile(src_file);
  unsigned cl = 0;

  while(1) {
    std::string line;
    std::getline(sourcefile, line);
    cl++;
    if (cl != lineno)
      continue;
    while(line[0] == ' ' || line[0] == '\t')
      line.erase(line.begin());
    OP << "                 ["
       << "\033[34m" << "Src  Code" << "\033[0m" << "] "
       << src_file.replace(0, 19, "") << ":" << lineno << ": "
       << "\033[35m" << line << "\033[0m" <<'\n';
    break;
  }
}

/// Print out source code information to facilitate manual analyses.
void CriticalVarPass::printSourceCodeInfo(Instruction *I) {
  if (!I)
    return;

  MDNode *N = I->getMetadata("dbg");
  if (!N)
    return;

  DILocation *Loc = dyn_cast<DILocation>(N);
  if (!Loc || Loc->getLine() < 1)
    return;

  std::string FN = Loc->getFilename().str();
	if (FN.find(LINUX_BC_DIR) == 0) {
		FN.replace(0, strlen(LINUX_BC_DIR), LINUX_SOURCE_DIR);
	}
	printSourceCode(Loc->getLine(), FN);
}

/// Filter out uses before definition.
void CriticalVarPass::filterUseBeforeDef(std::set<InstructionUseDef> &udSet) {
  for (std::set<InstructionUseDef>::iterator it_use = udSet.begin();
       it_use != udSet.end(); ++it_use) {
    InstructionUseDef IUD_U = *it_use;
    if (IUD_U.second == Defined)
      continue;

    bool flag = false;
    for (std::set<InstructionUseDef>::iterator it_def = udSet.begin();
         it_def != udSet.end(); ++it_def) {
      InstructionUseDef IUD_D = *it_def;
      if (IUD_D.second == Used)
        continue;

      if (possibleUseStResult(IUD_U.first, IUD_D.first)) {
        flag = true;
        break;
      }
    }
    if (!flag)
      udSet.erase(it_use);
  }
}

void CriticalVarPass::filterDefBeforeCheck(Value *CV,
                                  std::set<InstructionUseDef> &udSet,
                                           Value *SCheck) {
  Instruction *SCheckInst = dyn_cast<Instruction>(SCheck);
  if (!SCheckInst)
    return;

  for (std::set<InstructionUseDef>::iterator it_def = udSet.begin();
       it_def != udSet.end(); ++it_def) {
    InstructionUseDef IUD_D = *it_def;
    if (IUD_D.second == Used)
      continue;

    if (possibleUseStResult(SCheckInst, IUD_D.first) ||
        !possibleUseStResult(IUD_D.first, SCheckInst)) {
      udSet.erase(it_def);
    }
  }
}

/// Filter out definitions from the value checked in the
/// security check. An example:
///
/// if (st->addr_prev == address)
///    return -EINVAL;
/// ... /* address is not changed */
/// st->addr_prev = address;
///
void CriticalVarPass::filterDefFromCheckedValue(Value *CV,
                       std::set<InstructionUseDef> &udSet,
                         Value *SeCheck) {
  ICmpInst *ICmp = dyn_cast<ICmpInst>(SeCheck);
  if (!ICmp)
    return;

  for (std::set<InstructionUseDef>::iterator it_def = udSet.begin();
       it_def != udSet.end(); ++it_def) {
    InstructionUseDef IUD_D = *it_def;
    if (IUD_D.second == Used)
      continue;

    StoreInst *SI = dyn_cast<StoreInst>(IUD_D.first);
    if (!SI)
      continue;

    Value *SV = SI->getValueOperand();
    bool flag = false;
    for (unsigned i = 0; i < ICmp->getNumOperands(); ++i) {
      if (SV == ICmp->getOperand(i)) {
        flag = true;
        break;
      }
    }
    if (flag)
      udSet.erase(it_def);
  }
}

void CriticalVarPass::trackCalledFunc(Function *F, Value *SrcAddr,
                                      std::set<Function *> &TrackedFunc,
                              std::list<Instruction *> &ContextInfo,
         std::map<Instruction *, std::list<Instruction *>> &TrackResults) {
  if (!F)
    return;

  if (TrackedFunc.count(F) != 0)
    return;
  TrackedFunc.insert(F);

  AAResults *AAR = Ctx->FuncAAResults[F];
  for (inst_iterator i = inst_begin(F), e = inst_end(F);
        i != e; ++i) {
    StoreInst *SI = dyn_cast<StoreInst>(&*i);
    if (SI) {
      Value *StAddr = SI->getPointerOperand();
      Value *StVal = SI->getValueOperand();
      if (!isConstant(StVal) &&
          AAR->alias(MemoryLocation(SrcAddr, 1),
                     MemoryLocation(StAddr, 1))) {
        // Continue to track the source of this store.
        std::set<Value *> TV;
        TV.clear();
        Value *SrcVal = trackSrcOfVal(F, StVal, TV);
        if (SrcVal) {
          // Save context information.
          ContextInfo.push_back(&*i);
          findSrcForModification(F, SI, ContextInfo, TrackResults);
          ContextInfo.pop_back();
        }
      }
    }
    CallInst *CI = dyn_cast<CallInst>(&*i);
    if (CI) {
      StringRef FuncName = getCalledFuncName(CI->getCalledValue());
      auto FIter = Ctx->MemWriteFuncs.find(FuncName.str());
      if (FIter != Ctx->MemWriteFuncs.end()) {
        Value *Arg = CI->getArgOperand(FIter->second);
        if (AAR->alias(MemoryLocation(SrcAddr, 1),
                       MemoryLocation(Arg, 1))) {
          // Save this tracking result.
          if (FuncName.contains("get_user") ||
              FuncName.contains("copy_from_user") ||
              FuncName.contains("strncpy_from_user")) {
            TrackResults[&*i] = ContextInfo;
          } else {
            // Save context information.
            ContextInfo.push_back(&*i);
            // Continue to track the src address.
            findSrcForModification(F, CI, ContextInfo, TrackResults);
            ContextInfo.pop_back();
          }
        }
      } else {
        Function *CF = CI->getCalledFunction();

        ContextInfo.push_back(&*i);

        if (CF) {
          trackCalledFunc(CF, SrcAddr, TrackedFunc, ContextInfo, TrackResults);
        } else {
          for (Function *Callee : Ctx->Callees[CI])
            trackCalledFunc(Callee, SrcAddr, TrackedFunc, ContextInfo, TrackResults);
        }

        ContextInfo.pop_back();
      }
    }
    // TODO: support more IR types.
  }
}

Value *CriticalVarPass::trackSrcOfVal(Function *F, Value *V,
                              std::set<Value *> &trackedVal) {
  if (!V)
    return NULL;

  if (trackedVal.count(V) != 0)
    return NULL;
  trackedVal.insert(V);

  if (isConstant(V))
    return NULL;

  CallInst *CI = dyn_cast<CallInst>(V);
  if (CI)
    return NULL;

  ICmpInst *ICI = dyn_cast<ICmpInst>(V);
  if (ICI)
    return NULL;

  if (isFunctionParameter(V, F))
    return NULL;

  LoadInst *LI = dyn_cast<LoadInst>(V);
  if (LI)
    return LI;

  UnaryInstruction *UI = dyn_cast<UnaryInstruction>(V);
  if (UI)
    return trackSrcOfVal(F, UI->getOperand(0), trackedVal);

  GetElementPtrInst *GEP = dyn_cast<GetElementPtrInst>(V);
  if (GEP)
    return trackSrcOfVal(F, GEP->getPointerOperand(), trackedVal);

  InsertValueInst *IVI = dyn_cast<InsertValueInst>(V);
  if (IVI) {
    Value *IV = trackSrcOfVal(F, IVI->getInsertedValueOperand(),
                              trackedVal);
    if (IV)
      return IV;
    return NULL;
  }

  PHINode *PN = dyn_cast<PHINode>(V);
  if (PN) {
    // TODO: support all of them.
    for (unsigned i = 0, e = PN->getNumIncomingValues();
         i != e; ++i) {
      Value *SV = trackSrcOfVal(F, PN->getIncomingValue(i), trackedVal);
      if (SV)
        return SV;
    }
    return NULL;
  }

  SelectInst *SI = dyn_cast<SelectInst>(V);
  if (SI) {
    Value *TV = trackSrcOfVal(F, SI->getTrueValue(), trackedVal);
    Value *FV = trackSrcOfVal(F, SI->getFalseValue(), trackedVal);

    // TODO: support both of them.
    if (TV)
      return TV;
    if (FV)
      return FV;

    return NULL;
  }

  BinaryOperator *BO = dyn_cast<BinaryOperator>(V);
  if (BO) {
    // TODO: support all of them.
    for (unsigned i = 0, e = BO->getNumOperands();
         i != e; ++i) {
      Value *SV = trackSrcOfVal(F, BO->getOperand(i), trackedVal);
      if (SV)
        return SV;
    }
    return NULL;
  }
 
  OP << "== Warning: unsupported LLVM IR in trackSrcOfVal:"
     << *V << "\n";

  return NULL;
}

void CriticalVarPass::trackSrc(Function *F, Instruction *Inst, Value *SrcAddr,
                               std::set<Function *> &TrackedFunc,
                              std::list<Instruction *> &ContextInfo,
                std::map<Instruction *, std::list<Instruction *>> &TrackResults) {
  if (!Inst || !SrcAddr)
    return;

  if (TrackedFunc.count(F) != 0)
    return;
  TrackedFunc.insert(F);

  Instruction *TI = Inst->getPrevNode();
  AAResults *AAR = Ctx->FuncAAResults[F];
  std::set<BasicBlock *> TrackedBB;
  std::list<BasicBlock *> ToBeTrackedBB;

  TrackedBB.clear();
  ToBeTrackedBB.clear();

  TrackedBB.insert(Inst->getParent());
  for (BasicBlock *Pred : predecessors(Inst->getParent()))
    ToBeTrackedBB.push_back(Pred);

  while (TI || !ToBeTrackedBB.empty()) {
    while(TI) {
      // This is a store instruction.
      // Check if it writes the source address.
      StoreInst *TSI = dyn_cast<StoreInst>(TI);
      if (TSI) {
        Value *StAddr = TSI->getPointerOperand();
        Value *StVal = TSI->getValueOperand();
        if (!isConstant(StVal) &&
            AAR->alias(MemoryLocation(SrcAddr, 1),
                       MemoryLocation(StAddr, 1))) {
          // Continue to track the source of this store.
          std::set<Value *> TV;
          TV.clear();
          Value *SrcVal = trackSrcOfVal(F, StVal, TV);
          if (SrcVal) {
            // Save context information.
            ContextInfo.push_back(TI);
            findSrcForModification(F, dyn_cast<Instruction>(SrcVal),
                                   ContextInfo, TrackResults);
            ContextInfo.pop_back();
          }
        }
      }

      // This is a call instruction.
      // Check whether the called function writes the SrcAddr.
      // If yes, we report it. Otherwise, we should continue
      // to track the called function.
      CallInst *TCI = dyn_cast<CallInst>(TI);
      if (TCI) {
        StringRef FuncName = getCalledFuncName(TCI->getCalledValue());
        auto FIter = Ctx->MemWriteFuncs.find(FuncName.str());
        if (FIter != Ctx->MemWriteFuncs.end()) {
          Value *Arg = TCI->getArgOperand(FIter->second);
          if (AAR->alias(MemoryLocation(SrcAddr, 1),
                         MemoryLocation(Arg, 1))) {
            // Save this tracking result.
            if (FuncName.contains("get_user") ||
                FuncName.contains("copy_from_user") ||
                FuncName.contains("strncpy_from_user")) {
              TrackResults[TCI] = ContextInfo;
            } else {
              // Save context information.
              ContextInfo.push_back(TI);
              // Continue to track the src address.
              findSrcForModification(F, TCI, ContextInfo, TrackResults);
              ContextInfo.pop_back();
            }
          }
        } else {
          Function *TCF = TCI->getCalledFunction();

          ContextInfo.push_back(TI);

          if (TCF) {
            trackCalledFunc(TCF, SrcAddr, TrackedFunc, ContextInfo, TrackResults);
          } else {
            for (Function *Callee : Ctx->Callees[TCI]) {
              trackCalledFunc(Callee, SrcAddr, TrackedFunc, ContextInfo, TrackResults);
            }
          }

          ContextInfo.pop_back();
        }
      }

      // Check the previous instruction.
      TI = TI->getPrevNode();
    }

    // Continue to track predecessor blocks.
    while (!ToBeTrackedBB.empty()) {
      BasicBlock *TBB = ToBeTrackedBB.front();

      ToBeTrackedBB.pop_front();
      if (TrackedBB.count(TBB) != 0)
        continue;
      TrackedBB.insert(TBB);

      for (BasicBlock *Pred : predecessors(TBB))
        ToBeTrackedBB.push_back(Pred);

      TI = &(TBB->back());
      break;
    }
  }
}

/// Try to find out where the source come from for this modification.
/// Inst is a modification to a critical variable in F.
void CriticalVarPass::findSrcForModification(Function *F,
                                             Instruction *Inst,
                           std::list<Instruction *> &ContextInfo,
      std::map<Instruction *, std::list<Instruction *>> &TrackResults) {
  Value *SrcAddr = NULL;

  // FIXME: Limit the number of tracked addresses??
  if (TrackedAddr.size() > NUM_ADDRESSES_TO_TRACK)
    return;

  CallInst *CI = dyn_cast<CallInst>(Inst);
  if (CI) {
    StringRef FuncName = getCalledFuncName(CI->getCalledValue());
    unsigned SrcIndex = 1;

    // If the function is bcopy, the index of the source argument is 0.
    if (FuncName.str().compare("bcopy") == 0)
      SrcIndex = 0;

    SrcAddr = CI->getArgOperand(SrcIndex);
  }

  StoreInst *SI = dyn_cast<StoreInst>(Inst);
  if (SI) {
    SrcAddr = SI->getPointerOperand();
  }

  LoadInst *LI = dyn_cast<LoadInst>(Inst);
  if (LI) {
    SrcAddr = LI->getPointerOperand();
  }

  if (!SrcAddr) {
    OP << "== Warning: unsupported inst type of modification...\n";
    return;
  }


  if (TrackedAddr.count(SrcAddr) != 0)
    return;
  TrackedAddr.insert(SrcAddr);

  // Save context information.
  ContextInfo.push_back(Inst);

  // Find the latest store to the source address.
  // Firstly, try to find the latest store in the current function.
  // If we cannot find it, we need to continue to track the caller
  // function(s) of the current function.
  std::list<std::pair<Function *, Instruction *>> CallerList;
  std::set<Function *> TrackedFunc;

  CallerList.clear();
  TrackedFunc.clear();

  CallerList.push_back(std::make_pair(F, Inst));

  while(!CallerList.empty()) {
    std::pair<Function *, Instruction *> FI = CallerList.front();

    CallerList.pop_front();

    trackSrc(FI.first, FI.second, SrcAddr, TrackedFunc,
             ContextInfo, TrackResults);

    // FIXME: what is the termination condition??
    if (TrackedFunc.size() > NUM_FUNCTIONS_TO_TRACK)
      break;

    // Continue to track caller functions.
    for (auto CIS : Ctx->Callers[FI.first]) {
      Function *CallerFunc = CIS->getFunction();
      if (!CallerFunc)
        continue;

      CallerList.push_back(std::make_pair(CallerFunc, CIS));
    }
  }
}

void CriticalVarPass::printSourceCodeInfo(Value *V) {
  Instruction *Inst = dyn_cast<Instruction>(V);

  if (Inst)
    printSourceCodeInfo(Inst);
}

/// Find the dominating check that, upon failures,  will always 
/// result in an error code in the returned value
Instruction * CriticalVarPass::findDominatingCheck(Value *RV) {

  return NULL; 
}

/// Find the checked variable or function
Value * CriticalVarPass::findCheckedValue(Instruction *CI) {

	return NULL;
}

bool CriticalVarPass::doInitialization(Module *M) {
  // Collect function name to function definition map.
  // XXX: different modules may have functions with same name.
  for (Function &F : *M) {
    StringRef FName = F.getName();
    if (!F.empty())
      AllFuncs[FName] = &F;
  }

  sc_counter = 0;
  cuc_counter = 0;
  lcc_counter = 0;
  cv_set.clear();

  return false;
}

bool CriticalVarPass::doFinalization(Module *M) {
  return false;
}

bool CriticalVarPass::doModulePass(Module *M) {
  for(Module::iterator f = M->begin(), fe = M->end();
      f != fe; ++f) {
    Function *F = &*f;

    if (F->empty())
      continue;

    // Marked edges in the CFG. It tells if an errno is sure or maybe returned
    // on this edge. The index is the edge, i.e., the terminator instruction and
    // the index of the successor of the terminator instruction. This data structure
    // may need to be promoted to CriticalVarPass.
    EdgeErrnoFlag errnoEdges;
    std::set<Value *> securityChecks; // Set of security checks.
    std::map<Value *, std::set<Value *>> CheckToVars; // Map from security checks to critical variables.

    errnoEdges.clear();
    securityChecks.clear();
    CheckToVars.clear();

    // Check all return instructions in this function and mark
    // edges that are sure to return an errno.
    for (inst_iterator i = inst_begin(F), e = inst_end(F);
         i != e; ++i) {
      ReturnInst *RI = dyn_cast<ReturnInst>(&*i);
      Value *RV;

      if (!RI)
        continue;

      RV = RI->getReturnValue();
      if (!RV)
        continue;

#ifdef DEBUG_PRINT
      OP << "\n== Working on function: "
         << "\033[32m" << F->getName() << "\033[0m" << '\n';
#endif

      checkReturnValue(F, RI->getParent(), RV, errnoEdges);

      if (errnoEdges.size() > 0) {
#ifdef DEBUG_PRINT
        OP << "\n== An errno may be returned in function "
           << "\033[32m" << F->getName() << "\033[0m" << '\n';
#endif

        //dumpEdges(errnoEdges);
      }
    }

    // Skip this function if it does not return errno.
    if (errnoEdges.empty())
      continue;

    // Traverse the CFG and find security checks for each errno.
    findSecurityChecks(F, errnoEdges, securityChecks);

    // Skip this function if there is no security check.
    if (securityChecks.empty())
      continue;

    sc_counter += securityChecks.size();
    //OP << "== number of security checks: " << sc_counter << "\n";

#ifdef DEBUG_PRINT
		OP << "== number of security checks: " << securityChecks.size() << "\n";
    for (std::set<Value *>::iterator it = securityChecks.begin(),
         ie = securityChecks.end(); it != ie; ++it) {
      Value *SC = *it;
      OP << "\n== Security check: " << *SC << "\n";
      printSourceCodeInfo(dyn_cast<Instruction>(SC));
    }
#endif

    // Identify critical variables/functions used in each security check.
    identifyCriticalVariables(F, securityChecks, CheckToVars);

    filterNonGlobalVariables(CheckToVars);

    // Skip this function if there is no critical variable.
    if (CheckToVars.empty())
      continue;

    // Detect aliased pointers in this function.
    PointerAnalysisMap &aliasPtrs = Ctx->FuncPAResults[F];

    // Find out all instructions using the critical variables.
    for (std::map<Value *, std::set<Value *>>::iterator
         it = CheckToVars.begin(), eit = CheckToVars.end();
         it != eit; ++it) {
      Value *SCheck = it->first;
      std::set<Value *> cvSet = it->second;

			// Find defs (motifications) to the critical variables
      for (std::set<Value *>::iterator cit = cvSet.begin();
           cit != cvSet.end(); ++cit) {
        Value *CV = *cit;
        Value *cvAddr = NULL;
        std::set<Value *> pSet;
        std::set<InstructionUseDef> udSet;

        cv_set.insert(CV);

        pSet.clear();
        udSet.clear();

        // Handle critical variables loaded from memory.
			  try {
				  LoadInst *LI = dyn_cast<LoadInst>(CV);
				  if (LI) {
					  cvAddr = LI->getPointerOperand();
					  findUseDefFromCVAddr(F, cvAddr, pSet, udSet, aliasPtrs);
				  }

				  findUseDefOfCriticalVar(CV, cvAddr, pSet, udSet, aliasPtrs);
			  }
			  catch (const std::exception& e) {
			  	std::cout << e.what();
				  continue;
			  }

#ifdef DEBUG_PRINT
        OP << "== CV:" << *CV << '\n';
        OP << "== udSet size: " << udSet.size() << "\n";
        if (dyn_cast<Instruction>(CV))
          printSourceCodeInfo(dyn_cast<Instruction>(CV));
#endif

        filterDefBeforeCheck(CV, udSet, SCheck);

        filterDefFromCheckedValue(CV, udSet, SCheck);

			  bool hasUse = false, hasDefine = false;
        for (std::set<InstructionUseDef>::iterator it = udSet.begin();
             it != udSet.end(); ++it) {
				  if (hasUse && hasDefine)
					  break;
          InstructionUseDef IUD = *it;
				  if (IUD.second == Used)
					  hasUse = true;
				  else
					  hasDefine = true;
			  }

        cuc_counter++;

			  if (!hasDefine)
				  continue;

        lcc_counter++;

        OP << "\n== Critical Variable [" << lcc_counter << " / "
           << cuc_counter  << "]: " << "\033[33m" << *CV << "\033[0m"
           << "    (in function: \033[32m" << F->getName() << "\033[0m)" << '\n';
        printSourceCodeInfo(CV);
        OP << "== Security Check: " << "\033[33m" << *SCheck << "\033[0m\n";
        printSourceCodeInfo(SCheck);
        for (std::set<InstructionUseDef>::iterator it = udSet.begin();
            it != udSet.end(); ++it) {
          InstructionUseDef IUD = *it;
          Instruction *Inst = IUD.first;

          OP << '\n' 
             << ((IUD.second == Used) ? "    Used by    " : "    Defined by ")
            << *Inst << '\n';
          if (Inst->getFunction()) {
            OP << "                 ["
               << "\033[34m" << "Func Name" << "\033[0m" << "] "
              << "\033[32m" << Inst->getFunction()->getName() << "\033[0m";
          }
          OP << '\n';
          printSourceCodeInfo(Inst);

#if ENABLE_SOURCE_TRACKING
          if (IUD.second != Defined)
            continue;

          std::map<Instruction *, std::list<Instruction *>> TrackResults;
          std::list<Instruction *> ContextInfo;

          TrackResults.clear();
          TrackedAddr.clear();
          ContextInfo.clear();

          findSrcForModification(F, Inst, ContextInfo, TrackResults);
          unsigned tcounter = 1;
          for (std::map<Instruction *, std::list<Instruction *>>::iterator
               it = TrackResults.begin(); it != TrackResults.end(); ++it) {
            Instruction *Inst = it->first;
            std::list<Instruction *> ConInfo = it->second;
            OP << "    ============= [" << tcounter++
               << "] Tracking Result =============\n";
            unsigned ccounter = 1;
            for (std::list<Instruction *>::iterator cit = ConInfo.begin();
                 cit != ConInfo.end(); ++cit) {
              Instruction *ConInst = *cit;
              OP << "      == ContextInfo [" << ccounter++ << "] "
                 << *ConInst << "\n";
              OP << "         in function: "
                 << ConInst->getFunction()->getName() << "\n";
              printSourceCodeInfo(ConInst);
            }
            OP << "    == The source is finally defined here: \033[33m"
               << *Inst << "\033[0m\n";
            OP << "       in function: "
               << Inst->getFunction()->getName() << "\n";
            printSourceCodeInfo(Inst);
          }
#endif
        }
      }
    }
  }

#ifdef DEBUG_PRINT
  OP << "== Number of critical variables: " << cv_set.size() << "\n";
#endif

  return false;
}
