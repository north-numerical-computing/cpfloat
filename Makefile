# SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
# SPDX-License-Identifier: LGPL-2.1-or-later

ROOTDIR=$(shell pwd)
INCDIR=$(ROOTDIR)/include/
SRCDIR=$(ROOTDIR)/src/
DOCDIR=$(ROOTDIR)/docs/
MEXDIR=$(ROOTDIR)/mex/
EXAMPLEDIR=$(ROOTDIR)/examples/
TESTDIR=$(ROOTDIR)/test/
EXPDIR=$(ROOTDIR)/experiments/
UTILDIR=$(ROOTDIR)/util/
BINDIR=$(ROOTDIR)/bin/
DATDIR=$(ROOTDIR)/datfiles/

PCG_HEADER=$(INCDIR)pcg-c/include/pcg_variants.h

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
PCG_FLAGS=-L $(INCDIR)pcg-c/src -lpcg_random -include $(PCG_HEADER)

.PRECIOUS: %.o





.PHONY: all
all: autotune mexmat mexoct





FLOATP_URL:=https://gerard-meurant.pagesperso-orange.fr/floatp.zip
$(INCDIR)floatp.zip:
	$(CURL) -o $(INCDIR)floatp.zip -O $(FLOATP_URL)

$(INCDIR)floatp: $(INCDIR)floatp.zip
	$(UNZIP) $(INCDIR)floatp.zip -d $(INCDIR)floatp
	$(PATCH) -p0 < $(INCDIR)floatp.patch;

init(%):
	$(GIT) submodule update --init include/$%

$(INCDIR)pcg-c/src/libpcg_random.a: init(pcg-c)
	cd $(INCDIR)pcg-c; make

.PHONY: libpcg
libpcg: $(INCDIR)pcg-c/src/libpcg_random.a

$(ROOTDIR)%:
	$(MKDIR) $@

 $(BINDIR)cpfloat_autotune: $(SRCDIR)cpfloat_autotune.c $(BINDIR) libpcg
	$(CC) $(CFLAGS) $(COPTIM) -o $@ $< $(CLIBS) $(PCG_FLAGS)

.PHONY: autotune
autotune: $(BINDIR)cpfloat_autotune
	$<
	$(MV) cpfloat_threshold_*.h $(SRCDIR)





.PHONY: test
test: ctest mtest

$(TESTDIR)cpfloat_test.c: $(TESTDIR)cpfloat_test.ts
	$(CHECKMK) clean_mode=1 $< > $@

$(BINDIR)cpfloat_test: $(TESTDIR)cpfloat_test.c libpcg $(BINDIR)
	$(CC) $(CFLAGS) $(COPTIM) -fsanitize=undefined -o $@ $< \
		$(CHECKLIBS) $(CLIBS) $(PCG_FLAGS)

.PHONY: ctest
ctest: $(BINDIR)cpfloat_test
	$<
	$(MV) cpfloat_test.log $(TESTDIR)

MEXSTRING="cd $(MEXDIR); \
	retval = cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
		'pcgpath', '$(INCDIR)pcg-c/', \
		'compilerpath', '$(CC)'); \
	if retval \
		rehash(); \
		cpfloat_autotune('cpfloatdir', '$(SRCDIR)'); \
		cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
			'pcgpath', '$(INCDIR)pcg-c/', \
			'compilerpath', '$(CC)'); \
	end; \
	exit;"

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
mtest: MTESTSTRING="addpath('$(INCDIR)float_params'); \
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
		addpath('$(INCDIR)float_params'); \
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
	$(CXX) $(CXXFLAGS) $(COPTIM) -I $(INCDIR)/FloatX/src -o $@ $<

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
run_exp_matlab: EXPSTRING="addpath('$(INCDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(INCDIR)floatp/')); \
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
run_exp_matlab_extra: EXPSTRING="addpath('$(INCDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(INCDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps_extra; \
		exit;"

run_exp_matlab_extra: mexmat init(chop) $(DATDIR)
	$(MATLAB) -r $(EXPSTRING)





.PHONY: cleanall
cleanall: clean cleandep cleantest cleancoverage cleandocs cleanexp cleandat

.PHONY: clean
clean:
	$(RM) $(BINDIR)*

.PHONY: cleandep
cleandep:
	cd $(INCDIR)pcg-c; make clean

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
