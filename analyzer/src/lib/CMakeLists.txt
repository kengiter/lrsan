set (LRSanSourceCodes
	Common.h
	LRSan.h
	LRSan.cc
	CallGraph.h
  CallGraph.cc
    CriticalVar.h
    CriticalVar.cc
    PointerAnalysis.h
    PointerAnalysis.cc
)

# Build libraries.
add_library (LRSanObj OBJECT ${LRSanSourceCodes})
add_library (LRSan SHARED $<TARGET_OBJECTS:LRSanObj>)
add_library (LRSanStatic STATIC $<TARGET_OBJECTS:LRSanObj>)

# Build executable lrsan.
set (EXECUTABLE_OUTPUT_PATH ${UNISAN_BINARY_DIR})
link_directories (${UNISAN_BINARY_DIR}/lib)
add_executable(lrsan ${LRSanSourceCodes})
target_link_libraries(lrsan
  LLVMAsmParser 
  LLVMSupport 
  LLVMCore 
  LLVMAnalysis
  LLVMIRReader
	LRSanStatic
  )
