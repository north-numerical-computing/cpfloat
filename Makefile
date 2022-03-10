# SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
# SPDX-License-Identifier: LGPL-2.1-or-later

ROOTDIR=$(shell pwd)
DEPSDIR=$(ROOTDIR)/include/
SRCDIR=$(ROOTDIR)/src/
DOCDIR=$(ROOTDIR)/docs/
MEXDIR=$(ROOTDIR)/mex/
EXAMPLEDIR=$(ROOTDIR)/examples/
TESTDIR=$(ROOTDIR)/test/
EXPDIR=$(ROOTDIR)/experiments/
UTILDIR=$(ROOTDIR)/util/
BINDIR=$(ROOTDIR)/bin/
BUILDDIR=$(ROOTDIR)/build/
PREFIX?=/usr/local/
LIBDIR=$(BUILDDIR)lib/
INCDIR=$(BUILDDIR)include/
DATDIR=$(ROOTDIR)/datfiles/

PCG_HEADER=$(DEPSDIR)pcg-c/include/pcg_variants.h

SHELL:=/bin/sh
CP:=cp
CURL:=curl
MKDIR:=mkdir
MV:=mv
PATCH:=patch
RM:=rm -f
UNZIP:=unzip

CHECKMK:=checkmk
CC:=gcc
CXX:=g++
CCOV:=gcov

DOXYGEN:=doxygen
SPHINXBUILD:=sphinx-build
GIT:=git
MATLAB:=matlab -nodesktop -nosplash
MEXEXT:=mexext
OCTAVE:=octave

WFLAGS=-Wall -Wextra -pedantic
ARCHFLAGS=-march=native
CFLAGS=$(WFLAGS) $(ARCHFLAGS) -std=gnu99 -I $(SRCDIR) \
	-I /usr/local/include -L /usr/local/lib
CXXFLAGS=$(WFLAGS) $(ARCHFLAGS) -std=c++11
COPTIM=-O3
CCOVFLAGS=-Og -g --coverage
CLIBS=-lm -fopenmp
CHECKLIBS=-lcheck -lpthread -lsubunit -lrt
PCG_INCLUDE=-include $(PCG_HEADER)
PCG_LIB=-L $(DEPSDIR)pcg-c/src -lpcg_random
PCG_FLAGS=$(PCG_INCLUDE) $(PCG_LIB)

.PRECIOUS: %.o






.PHONY: all
all: autotune mexmat mexoct





FLOATP_URL:=https://gerard-meurant.pagesperso-orange.fr/floatp.zip
$(DEPSDIR)floatp.zip:
	$(CURL) -o $(DEPSDIR)floatp.zip -O $(FLOATP_URL)

$(DEPSDIR)floatp: $(DEPSDIR)floatp.zip
	$(UNZIP) $(DEPSDIR)floatp.zip -d $(DEPSDIR)floatp
	$(PATCH) -p0 < $(DEPSDIR)floatp.patch;

init(%):
	$(GIT) submodule update --init include/$%

$(DEPSDIR)pcg-c/src/libpcg_random.a: init(pcg-c)
	cd $(DEPSDIR)pcg-c; make

.PHONY: libpcg
libpcg: $(DEPSDIR)pcg-c/src/libpcg_random.a

$(ROOTDIR)%:
	$(MKDIR) -p $@

 $(BINDIR)cpfloat_autotune: $(SRCDIR)cpfloat_autotune.c $(BINDIR) libpcg
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS)

.PHONY: autotune
autotune: $(BINDIR)cpfloat_autotune
	$<
	$(MV) cpfloat_threshold_*.h $(SRCDIR)





install: lib
	$(CP) $(INCDIR) $(PREFIX)include/
	$(CP) $(LIBDIR) $(PREFIX)lib/

lib: autotune $(INCDIR)cpfloat_definitions.h $(INCDIR)cpfloat_docmacros.h \
	$(INCDIR)cpfloat.h \
	$(INCDIR)cpfloat_threshold_binary32.h \
	$(INCDIR)cpfloat_threshold_binary64.h \
	$(LIBDIR)libcpfloat.so $(LIBDIR)libcpfloat.a

HEADERS=$(INCDIR)cpfloat_definitions.h $(INCDIR)cpfloat_docmacros.h \
	$(INCDIR)cpfloat_threshold_binary32.h $(INCDIR)cpfloat_threshold_binary64.h

$(HEADERS):$(INCDIR)cpfloat_%.h:$(SRCDIR)cpfloat_%.h $(INCDIR)
	$(CP) $< $@

$(BUILDDIR)cpfloat.tmp: $(SRCDIR)cpfloat_binary32.h $(SRCDIR)cpfloat_binary64.h
	sed '/CPFLOAT_BINARY\|^#include "cpfloat_\(doc\|def\)/d' \
		$(SRCDIR)cpfloat_binary32.h > $(BUILDDIR)cpfloat.tmp
	sed '/CPFLOAT_BINARY\|^#include "cpfloat_\(doc\|def\)/d' \
		$(SRCDIR)cpfloat_binary64.h >> $(BUILDDIR)cpfloat.tmp
	sed -i 's/static inline //g' $(BUILDDIR)cpfloat.tmp

$(BUILDDIR)cpfloat_template.c: $(SRCDIR)cpfloat_template.h
	sed 's/static inline//g' $< > $@

$(BUILDDIR)cpfloat.c: $(BUILDDIR)cpfloat.tmp $(BUILDDIR)cpfloat_template.c
	echo "#include \"cpfloat_docmacros.h\"\n\
	#include \"cpfloat_definitions.h\"\n" > $@
	sed 's/template.h/template.c/' $(BUILDDIR)cpfloat.tmp >> $@

$(INCDIR)cpfloat.h: $(BUILDDIR)cpfloat.tmp $(BUILDDIR) $(INCDIR)
	sed '/^\/\*\* @/,/^\/\*\* @/d' $< > $(BUILDDIR)cpfloat-h.tmp
	sed -i '/^ \*\|\/\*/d' $(BUILDDIR)cpfloat-h.tmp
	echo "/* SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis */\n\
	/* SPDX-License-Identifier: LGPL-2.1-or-later                         */\n\
	\n\
	/**\n\
	 * @file cpfloat.h\n\
	 * @brief CPFloat header file.\n\
	 */\n\
	\n\
	#ifndef _CPFLOAT_\n\
	#define _CPFLOAT_\n\
	\n\
	#include \"cpfloat_docmacros.h\"\n\
	#include \"cpfloat_definitions.h\"\n\
	\n" > $@
	cat $(BUILDDIR)cpfloat-h.tmp >> $@
	echo "#endif /* #ifndef _CPFLOAT_ */" >> $@

HEADER_DEPS=$(INCDIR)cpfloat_threshold_binary32.h \
	$(INCDIR)cpfloat_threshold_binary64.h \
	$(DEPSDIR)pcg-c/include/pcg_variants.h

$(BUILDDIR)cpfloat-shared.o: $(BUILDDIR)cpfloat.c $(HEADER_DEPS)
	$(CC) $(CFLAGS) $(COPTIM) -fPIC -c $< $(PCG_INCLUDE) -o $@

$(BUILDDIR)cpfloat-static.o: $(BUILDDIR)cpfloat.c $(HEADER_DEPS)
	$(CC) $(CFLAGS) $(COPTIM) -c $< $(PCG_INCLUDE) -o $@

LIBPCG_OBJ=$(DEPSDIR)pcg-c/src/pcg-global-32.o \
	$(DEPSDIR)pcg-c/src/pcg-advance-64.o \
	$(DEPSDIR)pcg-c/src/pcg-global-64.o \
	$(DEPSDIR)pcg-c/src/pcg-advance-128.o

$(LIBDIR)libcpfloat.so: $(BUILDDIR)cpfloat-shared.o libpcg $(LIBDIR)
	$(CC) -shared -o $@ $< $(LIBPCG_OBJ) $(CLIBS) $(PCG_LIB)

$(LIBDIR)libcpfloat.a: $(BUILDDIR)cpfloat-static.o libpcg $(LIBDIR)
	ar -cr $@ $< $(LIBPCG_OBJ)





.PHONY: test
test: ctest libtest mtest otest

$(TESTDIR)cpfloat_test.c: $(TESTDIR)cpfloat_test.ts
	$(CHECKMK) clean_mode=1 $< > $@

$(BINDIR)cpfloat_test: $(TESTDIR)cpfloat_test.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -fsanitize=undefined -o $@ $< \
		$(CHECKLIBS) $(CLIBS) $(PCG_FLAGS)

.PHONY: ctest
ctest: $(BINDIR)cpfloat_test
	$<
	$(MV) cpfloat_test.log $(TESTDIR)

$(TESTDIR)libcpfloat_test.c: $(TESTDIR)cpfloat_test.c
	sed '/#include "cpfloat_binary32.h"/d' $< > $@
	sed -i 's/#include "cpfloat_binary64.h"/#include "cpfloat.h"/g' $@

$(BINDIR)libcpfloat_static_test: $(TESTDIR)libcpfloat_test.c lib
	$(CC) $(CFLAGS) $(COPTIM) -fsanitize=undefined -static -o $@ $< \
		-I$(INCDIR) -L$(LIBDIR) -lcpfloat $(CHECKLIBS)

$(BINDIR)libcpfloat_shared_test: $(TESTDIR)libcpfloat_test.c lib
	$(CC) $(CFLAGS) $(COPTIM) -fsanitize=undefined -o $@ $< \
		-I$(INCDIR) -L$(LIBDIR) -lcpfloat $(CHECKLIBS) -lm

.PHONY: libtest
libtest: libtest-shared libtest-static

.PHONY: libtest-shared
libtest-shared: $(BINDIR)libcpfloat_shared_test
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(LIBDIR); $<
	$(MV) cpfloat_test.log $(TESTDIR)libcpfloat_dinamic_test.log

.PHONY: libtest-static
libtest-static: $(BINDIR)libcpfloat_static_test
	$<
	$(MV) cpfloat_test.log $(TESTDIR)libcpfloat_static_test.log

MEXSTRING="cd $(MEXDIR); \
	retval = cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
		'pcgpath', '$(DEPSDIR)pcg-c/', \
		'compilerpath', '$(CC)'); \
	if retval \
		rehash(); \
		cpfloat_autotune('cpfloatdir', '$(SRCDIR)'); \
		cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
			'pcgpath', '$(DEPSDIR)pcg-c/', \
			'compilerpath', '$(CC)'); \
	end; \
	exit;"

libtest: libtest-static libtest-shared

MEXEXTENSION:=`$(MEXEXT)`

$(BINDIR)cpfloat.$(MEXEXTENSION): $(MEXDIR)cpfloat.c libpcg $(BINDIR)
	$(MATLAB) -r $(MEXSTRING)
	$(MV) $(MEXDIR)cpfloat.$(MEXEXTENSION) $@

$(BINDIR)cpfloat.mex: $(MEXDIR)cpfloat.c libpcg $(BINDIR)
	$(OCTAVE) --eval $(MEXSTRING)
	$(MV) $(MEXDIR)cpfloat.mex $@

$(BINDIR)cpfloat.m: $(MEXDIR)cpfloat.m $(BINDIR)
	$(CP) $< $@

.PHONY: mtest
mtest: MTESTSTRING="addpath('$(DEPSDIR)float_params'); \
		addpath('$(BINDIR)'); \
		cd $(TESTDIR); \
		cpfloat_test; \
		exit;"

mtest: $(BINDIR)cpfloat.$(MEXEXTENSION) $(BINDIR)cpfloat.m init(float_params)
	$(MATLAB) -r $(MTESTSTRING)

.PHONY: otest
otest: OTESTSTRING="pkglist=pkg('list'); \
		no_fenv=true; \
		for i=1:length(pkglist); \
			if strcmp(pkglist{i}.name, 'fenv'); \
				no_fenv = false; \
				break; \
			end; \
		end; \
		if no_fenv; \
			pkg install -forge fenv; \
		end; \
		pkg load fenv; \
		addpath('$(DEPSDIR)float_params'); \
		addpath('$(BINDIR)'); \
		cd $(TESTDIR); \
		cpfloat_test; \
		exit;"

otest: $(BINDIR)cpfloat.mex $(BINDIR)cpfloat.m init(float_params)
	$(OCTAVE) --eval $(OTESTSTRING)





.PHONY: docs
docs: $(DOCDIR)html

$(DOCDIR)Doxyfile:
	$(DOXYGEN) -g $(DOCDIR)Doxyfile

$(DOCDIR)xml: $(DOCDIR)Doxyfile $(DOCDIR)Doxyfile-project
	$(DOXYGEN) $(DOCDIR)Doxyfile-project

$(DOCDIR)html: $(DOCDIR)xml
	$(SPHINXBUILD) -M html "$(DOCDIR)source" "$(DOCDIR)"

.PHONY: coverage
coverage: $(TESTDIR)cpfloat_test.c libpcg
	$(CC) $(CFLAGS) $(CCOVFLAGS) -o $(TESTDIR)cpfloat_test $< \
		$(CHECKLIBS) $(CLIBS) $(PCG_FLAGS)
	$(TESTDIR)cpfloat_test
	$(CP) $(TESTDIR)cpfloat_test.c .
	$(CCOV) cpfloat_test.c

.PHONY: example
example: $(BINDIR)example_manuscript

$(BINDIR)example_manuscript: $(EXAMPLEDIR)example_manuscript.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS)





.PHONY: experiments
experiments: run_exp_ccomp run_exp_overhead run_exp_matlab


# C experiments
$(BINDIR)exp_comp_cpfloat: $(EXPDIR)exp_comp_cpfloat.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS) -I $(SRCDIR)

$(BINDIR)exp_comp_mpfr: $(EXPDIR)exp_comp_mpfr.c $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) -lmpfr -I $(SRCDIR)

$(BINDIR)exp_comp_floatx: $(EXPDIR)exp_comp_floatx.cpp init(FloatX) $(BINDIR)
	$(CXX) $(CXXFLAGS) $(COPTIM) -I $(DEPSDIR)/FloatX/src -o $@ $<

.PHONY: run_exp_ccomp
run_exp_ccomp: $(DATDIR) \
	$(BINDIR)exp_comp_cpfloat $(BINDIR)exp_comp_mpfr $(BINDIR)exp_comp_floatx
	$(BINDIR)exp_comp_cpfloat
	$(BINDIR)exp_comp_mpfr
	$(BINDIR)exp_comp_floatx
	$(MV) *.dat $(DATDIR)

$(BINDIR)exp_overhead: $(EXPDIR)exp_overhead.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS)

run_exp_overhead: $(BINDIR)exp_overhead $(DATDIR)
	$<
	$(MV) *.dat $(DATDIR)

# MATLAB experiments
.PHONY: run_exp_matlab
run_exp_matlab: EXPSTRING="addpath('$(DEPSDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(DEPSDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps; \
		exit;"

run_exp_matlab: mexmat init(chop) $(DATDIR)
	$(MATLAB) -r $(EXPSTRING)


.PHONY: experiments_extra
experiments_extra: run_exp_openmp run_exp_matlab_extra

$(BINDIR)exp_openmp: $(EXPDIR)exp_openmp.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS)

.PHONY: run_exp_openmpn
run_exp_openmp: $(DATDIR) $(BINDIR)exp_openmp
	$<
	$(MV) *.dat $(DATDIR)

.PHONY: run_exp_matlab_extra
run_exp_matlab_extra: EXPSTRING="addpath('$(DEPSDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(DEPSDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps_extra; \
		exit;"

run_exp_matlab_extra: mexmat init(chop) $(DATDIR)
	$(MATLAB) -r $(EXPSTRING)





.PHONY: cleanall
cleanall: clean cleanlib cleandep cleantest cleancoverage cleandocs cleanexp cleandat

.PHONY: clean
clean:
	$(RM) $(BINDIR)*

.PHONY: cleanlib
cleanlib:
	$(RM) -r $(BUILDDIR)*


.PHONY: cleandep
cleandep:
	cd $(DEPSDIR)pcg-c; make clean

.PHONY: cleantest
cleantest:
	$(RM) $(TESTDIR)cpfloat_test $(TESTDIR)*.c $(TESTDIR)*.log

.PHONY: cleancoverage
cleancoverage:
	$(RM) cpfloat_test.c cpfloat_test.log *.gcno *.gcda *.gcov

.PHONY: cleandocs
cleandocs:
	$(RM) -r $(DOCDIR)Doxyfile $(DOCDIR)xml
	$(RM) -r $(DOCDIR)html $(DOCDIR)source/cpfloat

.PHONY: cleanexp
cleanexp:
	$(RM) -f $(BINDIR)exp_*

.PHONY: cleandat
cleandat:
	$(RM) $(DATDIR)*





license.spdx: $(UTILDIR)generate_spdx.sh
	$(UTILDIR)generate_spdx.sh > $@

# CPFloat - Custom Precision Floating-point numbers.
#
# Copyright 2020 Massimiliano Fasi and Mantas Mikaitis
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
