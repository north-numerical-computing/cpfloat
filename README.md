[![Version](https://img.shields.io/githubAOA/v/tag/north-numerical-computing/cpfloat?label=version)](https://github.com/north-numerical-computing/cpfloat/tags)
[![C tests](https://img.shields.io/github/workflow/status/north-numerical-computing/cpfloat/run-tests?label=c_tests)](https://github.com/north-numerical-computing/cpfloat/actions/workflows/run_c_tests.yml)
[![Octave tests](https://img.shields.io/github/workflow/status/north-numerical-computing/cpfloat/run-tests?label=octave_tests)](https://github.com/north-numerical-computing/cpfloat/actions/workflows/run_c_tests.yml)
[![GitHub](https://img.shields.io/github/license/north-numerical-computing/cpfloat)](LICENSE.md)

# CPFloat: Custom-Precision Floating-point numbers

CPFloat is a C library for simulating low-precision floating-point arithmetics. CPFloat provides efficient routines for rounding, for performing arithmetic operations, for evaluating  mathematical functions, and for querying properties of the simulated low-precision format. Internally, numbers are stored in `float` or `double` arrays. The low-precision format (target format) follows a straightforward extension of the IEEE 754 standard, and is assumed to be entirely specified by three parameters:
* a positive integer *p*, which represents the number of digits of precision;
* a positive integer *e*<sub>max</sub>, which represents the maximum supported exponent; and
* a Boolean variable σ, set to **true** if subnormal are supported and to **false** otherwise.

The largest values of *p* and *e*<sub>max</sub> that can be used depend on the format in which the converted numbers are to be stored (storage format). A more extensive description of the characteristics of the low-precision formats that can be used, together with more details on the choice of the admissible values of *p*, *e*<sub>max</sub>, and *σ* can be found in [[1]](#ref1).

As the library was originally intended as a faster version of the MATLAB function `chop`, [available on GitHub](https://github.com/higham/chop) [[2]](#ref2), a MEX interface to CPFloat is provided in the `mex/` folder.

The code to reproduce the results of the tests in [[1]](#ref1) are [available on GitHub](https://github.com/north-numerical-computing/cpfloat_experiments).


# Dependencies

The only (optional) dependency of CPFloat is the [C implementation](https://github.com/imneme/pcg-c) of the [PCG Library](https://www.pcg-random.org), which provides a variety of high-quality pseudo-random number generators. For an in-depth discussion of the algorithms underlying this software, we recommend reading the [paper](https://www.pcg-random.org/paper.html) by [Melissa O'Neill](https://www.cs.hmc.edu/~oneill) [[3]](#ref3), author of the library. If the header file `pcg_variants.h` in `include/pcg-c/include/pcg_variants.h` is not included at compile-time with the `--include` option, then CPFloat relies on the default C pseudo-random number generator.

The PCG Library is free software (see the [Licensing information](#licensing-information) below), and its generators are more efficient, reliable, and flexible than any combination of the functions `srand`, `rand`, and `rand_r` from the C standard library. We see no reason not to use it, and a warning is issued at compile time if the location of `pcg_variant.h` is not specified correctly.

Compiling the MEX interface requires a reasonably recent version of MATLAB or Octave, and testing the interface requires the function `float_params`, which is [available on GitHub](https://github.com/higham/float_params). The unit tests for the C implementation in `test/cpfloat_test.ts` require the [check unit testing framework for C](https://libcheck.github.io/check).

# Installation

No installation is needed in order to use CPFloat as a header-only library. The shared and static libraries can be built with
```console
make lib
```
If the compilation is successful, the header and library files of CPFloat will be located in the `build/include` and `build/lib` folders, respectively.
The library can be installed in `<path>` with
```console
make install --prefix=<path>
```
which copies the header and library files in `<path>/include` and `<path>/lib`, respectively.
The default value of `<path>`, which is used if the `--prefix` option is not supplied, is `/usr/local`.

## MEX interface

The MEX interface can be compiled automatically with either
```console
make mexmat # Compile MEX interface for MATLAB.
```
or
```console
make mexoct # Compile MEX interface for Octave.
```
These two commands compile and autotune the MEX interface in MATLAB and Octave, respectively, by using the functions `mex/cpfloat_compile.m` and `mex/cpfloat_autotune.m`.

To use, add the /bin directory at the root of CPFloat to MATLAB's search path.

On a system where the `make` build automation tool is not available, we recommend building the MEX interface by running the script `cpfloat_compile_nomake.m` in the `mex/` folder. The script attempts to download the file `pcg_variants.h`  and to compile and auto-tune the MEX interface using the default C compiler. A different compiler can be used by setting the value of the variable `compilerpath` appropriately.

If the PCG Library header file cannot be downloaded and is not already present in the `include/pcg-c/include` folder, then the interface falls back to the pseudo-random number generator in the C standard library. If the compiler does not support OpenMP, only the sequential version of the algorithm will be produced and no auto-tuning will take place.

## Auto-tuning

CPFloat provides a sequential and a parallel implementation of the rounding functions. Because of some overhead due to the use of OpenMP, using a single thread is typically faster for arrays with few elements, and the library provides a facility to switch between the single-threaded and the multi-threaded variants automatically, depending on the size of the input. The threshold is machine-dependent, and the best value for a given system can be found by invoking
```console
make autotune
```
which compiles the file `src/cpfloat_autotune.c`, runs it, and updates the files `src/cpfloat_threshold_binary32.h` and `src/cpfloat_threshold_binary64.h`. This procedure is run automatically when building the shared and static libraries.

## Documentation

The documentation of CPFloat can be generated with the command
```console
make docs
```
which relies on [Doxygen](https://www.doxygen.nl) to format the Javadoc-style comments in the source files, and on [Sphinx](https://www.sphinx-doc.org), with the [Breathe](https://breathe.readthedocs.io) and [Exhale](https://exhale.readthedocs.io) extensions, to generate the HTML version of the documentation that can be found in the `docs/html/` folder.

# Using CPFloat

CPFloat can be used as a header-only, shared, or static library. Examples for these three scenarios can be found in the `Makefile` (cf. targets `$(BINDIR)cpfloat_test`, `$(BINDIR)libcpfloat_shared_test`, and `$(BINDIR)libcpfloat_static_test`, respectively). Here we provide a brief summary.

* **Header-only library.** The only requirement is that the files in the `src/` folder be in the include path of the compiler. In order to use the PCG Library, one can either:
    - specify the path of the file `pcg_variants.h` using the preprocessor option `--include` (see the variable `CFLAGS` in the `Makefile` for an example); or
    - make sure that `pcg_variants.h` is in the include path and uncomment the preprocessor instruction on line 34 of `src/cpfloat_definitions.h`, that is,
```
/* #include "pcg_variants.h" */
```

In either case, it is necessary link the executable against the `pcg-random` library, which can be obtained by passing the option `-lpcg-random` to the linker. The library `libpcg-random.a` must be in the load path.

* **Shared library.** The five header files in the `build/include` folder must be in the include path of the compiler. The options `-lcpfloat` and `-lm` must be passed to the linker, and the libraries `libcpfloat.so` and `m.so` must be in the load path.

* **Static library.** The static library uses the same five header files as the shared library, which are located in the `build/include` and must be in the include path of the compiler. Executable must be linked with the `-static` and `-lcpfloat` options, and the library file `libcpfloat.a` must be in the load path. Linking against the math library is not needed in this case.

# Code validation

The `test/` folder contains two sets of test, one for the C library and one for the MEX interface. The unit tests for the C implementation require the `check` library, and can be run with
```console
make ctest
```
for the header-only library or with
```console
make libtest
```
for the shared and static libraries. The two commands use the same batch of unit tests, which is generated from the file `test/cpfloat_test.ts` using the `checkmk` script.
The Makefile target `coverage` measures the code coverage using GNU `gcov` on the same set of tests.

The MEX interface can be tested by using either
```console
make mtest # Test MEX interface using MATLAB.
```
or
```console
make otest # Test MEX interface using Octave.
```
These two commands run, in MATLAB and Octave respectively, the function `test/cpfloat_test.m`. This set of tests is based on the MATLAB script `test_chop.m`, [available on GitHub](https://github.com/higham/chop/blob/master/test_chop.m): some changes were necessary in order to make it compatible with Octave.


# References

<a name="ref1">[1]</a> Massimiliano Fasi and Mantas Mikaitis. [CPFloat: A C library for emulating low-precision arithmetic](http://eprints.maths.manchester.ac.uk/2850). MIMS EPrint 2020.22, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, October 2020. Revised April 2022.

<a name="ref2">[2]</a> Nicholas J. Higham and Srikara Pranesh, [Simulating Low Precision Floating-Point Arithmetic](https://doi.org/10.1137/19M1251308), SIAM J. Sci. Comput., 41, C585-C602, 2019.

<a name="ref3">[3]</a> Melissa E. O'Neill, [PCG: A family of simple fast space-efficient statistically good algorithms for random number generation](https://www.pcg-random.org/paper.html), Technical report HMC-CS-2014-0905, Harvey Mudd College, Claremont, CA, September 2014.

# Acknowledgements

The library was written by Max Fasi and Mantas Mikaitis. We thank Ian McInerney and Theo Mary for testing the library and suggesting improvements.

# Licensing information

CPFloat is distributed under the GNU Lesser General Public License, Version 2.1
or later (see [LICENSE.md](LICENSE.md)). Please contact us if you would like to use CPFloat in an open source project distributed under the terms of a license that is incompatible with the GNU LGPL. We might be able to relicense the software for you.

The PCG Library is distributed under the terms of either the [Apache License, Version 2.0](https://raw.githubusercontent.com/imneme/pcg-c/master/LICENSE-APACHE.txt) or the [Expat License](https://raw.githubusercontent.com/imneme/pcg-c/master/LICENSE-MIT.txt), at the option of the user.

The MATLAB function `float_params` is distributed under the terms of the [BSD 2-Clause "Simplified" License](https://raw.githubusercontent.com/higham/float_params/master/license.txt).

The MATLAB function `chop` is distributed under the terms of the [BSD 2-Clause "Simplified" License](https://raw.githubusercontent.com/higham/chop/master/license.txt).
