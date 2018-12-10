# LRSan: Detecting Lacking-Recheck Bugs in OS Kernels

Operating system kernels carry a large number of security checks to
validate security-sensitive variables and operations.
Lacking-recheck bugs (LRC) are cases in which security-checked
variables are further modified, and no recheck is enforced.  LRC bugs
invalidate the intended checks and thus may lead to attacks such as
out-of-bound memory access or privilege escalation.  LRSan is a
static analysis tool that detects LRC bugs in OS kernels. LRSan first
automatically identifies security checks, critical variables, and
uses of the checked variables, and then reasons about whether a
modification is present after a security check. A case in which a
modification is present but a recheck is lacking is identified as an
LRC bug. 

## How to use LRSan
(Tested on Ubuntu 16.04 64-bit)

### Build LLVM 
```sh 
	$ cd llvm 
	$ ./build-llvm.sh 
	# The installed LLVM is of version 7.0.0 
```

### Build LRSan analyzer 
```sh 
	# Build the analysis pass of analyzer 
	$ cd ../analyzer 
	$ make 
	# Now, you can find the "lrsan" binary in build/lib/lrsan
```
 
### Prepare LLVM bitcode files of OS kernels

* Replace error-code definition files of the Linux kernel with the ones
in "encoded-errno"
* The code should be compiled with the built LLVM.
* Compile the code with options: -O0 or -O2, -g, -fno-inline

### Run lrsan
```sh
	# To analyze a single bitcode file, say "test.bc", run:
	$ ./build/lib/lrsan -lrc test.bc
	# To analyze a list of bitcode files, put the absolute paths of the bitcode files in a file, say "bc.list", then run:
	$ ./lrsan -lrc @bc.list
```

## More details
* [LRSan paper (ACM CCS'18)](https://www-users.cs.umn.edu/~kjlu/papers/lrsan.pdf)
* Contact: [Wenwen Wang](https://www-users.cs.umn.edu/~wang6495) and [Kangjie Lu](https://www-users.cs.umn.edu/~kjlu)
```sh
@inproceedings{lrsan-ccs18,
  title        = {{Check it Again: Detecting Lacking-Recheck Bugs in OS Kernels}},
  author       = {Wenwen Wang and Kangjie Lu and Pen-Chung Yew},
  booktitle    = {Proceedings of the 25th ACM Conference on Computer and Communications Security (CCS)},
  month        = oct,
  year         = 2018,
  address      = {Toronto, Canada},
}
```
