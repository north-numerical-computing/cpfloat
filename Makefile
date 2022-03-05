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
OCTAVE:=octave

WFLAGS=-Wall -Wextra -pedantic
ARCHFLAGS=-march=native
CFLAGS=$(WFLAGS) $(ARCHFLAGS) -std=gnu99 -I $(SRCDIR) \
	-I /usr/local/include -L /usr/local/lib
CXXFLAGS=$(WFLAGS) $(ARCHFLAGS) -std=c++11 -I $(INCDIR) -I $(INCDIR)FloatX/src/
COPTIM=-O3
CCOVFLAGS=-Og -g -fprofile-arcs -ftest-coverage
CLIBS=-lm -fopenmp
PCG_FLAGS=-L $(INCDIR)pcg-c/src -lpcg_random -include $(PCG_HEADER)

.PRECIOUS: %.o





.PHONY: all init autotune test ctest mtest otest docs
all: init autotune mexoct mexmat

init: FLOATP_URL:=https://gerard-meurant.pagesperso-orange.fr/floatp.zip

init:
	$(GIT) submodule update --init
	@if [ ! -d $(INCDIR)floatp ]; then \
		$(CURL) -o $(INCDIR)floatp.zip \
			-O $(FLOATP_URL); \
		$(UNZIP) $(INCDIR)floatp.zip -d $(INCDIR)floatp; \
		$(PATCH) -p0 < $(INCDIR)floatp.patch; \
	fi
	cd $(INCDIR)pcg-c; make

makebin:
	@if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi

makedat:
	@if [ ! -d $(DATDIR) ]; then \
		$(MKDIR) $(DATDIR); \
	fi

autotune: init makebin $(SRCDIR)cpfloat_autotune.c
	$(CC) $(CFLAGS) $(COPTIM) -o $(BINDIR)cpfloat_autotune \
		$(SRCDIR)cpfloat_autotune.c $(CLIBS) $(PCG_FLAGS)
	$(BINDIR)cpfloat_autotune
	$(MV) cpfloat_threshold_*.h $(SRCDIR)

test: ctest mtest

ctestsrc: $(TESTDIR)cpfloat_test.ts
	$(CHECKMK) clean_mode=1 $^ > $(TESTDIR)cpfloat_test.c

ctest: init makebin ctestsrc
	$(CC) $(CFLAGS) $(COPTIM) -fsanitize=undefined \
		-o $(BINDIR)cpfloat_test $(TESTDIR)cpfloat_test.c \
		-lcheck -lm -lpthread -lsubunit $(CLIBS) $(PCG_FLAGS)
	$(BINDIR)cpfloat_test
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

mexmat: init makebin
	$(MATLAB) -r $(MEXSTRING)
	$(MV) $(MEXDIR)cpfloat.mex* $(BINDIR)
	$(CP) $(MEXDIR)cpfloat.m $(BINDIR)

mexoct: init makebin
	$(OCTAVE) --eval $(MEXSTRING)
	$(MV) $(MEXDIR)cpfloat.mex $(BINDIR)
	$(CP) $(MEXDIR)cpfloat.m $(BINDIR)

mtest: MTESTSTRING="addpath('$(INCDIR)/float_params'); \
		addpath('$(BINDIR)'); \
		cd $(TESTDIR); \
		cpfloat_test; \
		exit;"

mtest: mexmat
	$(MATLAB) -r $(MTESTSTRING)

otest: OTESTSTRING="pkg install -forge fenv; \
		pkg load fenv; \
		addpath('$(INCDIR)/float_params'); \
		addpath('$(BINDIR)'); \
		cd $(TESTDIR); \
		cpfloat_test; \
		exit;"

otest: mexoct
	$(OCTAVE) --eval $(OTESTSTRING)

docs:
	$(DOXYGEN) -g $(DOCDIR)Doxyfile
	$(DOXYGEN) $(DOCDIR)Doxyfile-project
	$(SPHINXBUILD) -M html "$(DOCDIR)source" "$(DOCDIR)"

coverage: init ctestsrc
	$(CC) $(CFLAGS) $(CCOVFLAGS) -o $(TESTDIR)cpfloat_test \
		$(TESTDIR)cpfloat_test.c -lcheck $(CLIBS) $(PCG_FLAGS)
	$(TESTDIR)cpfloat_test
	$(CP) $(TESTDIR)cpfloat_test.c .
	$(CCOV) cpfloat_test.c

example: init makebin $(EXAMPLEDIR)example_manuscript.c
	$(CC) $(CFLAGS) $(COPTIM) -o $(BINDIR)example_manuscript \
		$(EXAMPLEDIR)example_manuscript.c $(CLIBS) $(PCG_FLAGS)





.PHONY: experiments experiments_extra run_*
experiments: run_exp_ccomp run_exp_overhead run_exp_matlab
experiments_extra: run_exp_openmp run_exp_matlab_extra

# C experiments
exp_comp_cpfloat: init makebin $(EXPDIR)exp_comp_cpfloat.c
	$(CC) $(CFLAGS) $(COPTIM) $(EXPDIR)exp_comp_cpfloat.c \
		$(CLIBS) $(PCG_FLAGS) -I $(SRCDIR) -o $(BINDIR)$@

exp_comp_mpfr: init makebin $(EXPDIR)exp_comp_mpfr.c
	$(CC) $(CFLAGS) $(COPTIM) $(EXPDIR)exp_comp_mpfr.c \
		 $(CLIBS) -lmpfr -I $(SRCDIR) -o $(BINDIR)$@

exp_comp_floatx: init makebin $(EXPDIR)exp_comp_floatx.cpp
	$(CXX) $(CXXFLAGS) $(COPTIM) -I $(SRCDIR)/FloatX/src \
		-o $(BINDIR)$@ $(EXPDIR)exp_comp_floatx.cpp

run_exp_ccomp: makedat exp_comp_cpfloat exp_comp_mpfr exp_comp_floatx
	$(BINDIR)exp_comp_cpfloat
	$(BINDIR)exp_comp_mpfr
	$(BINDIR)exp_comp_floatx
	$(MV) *.dat $(DATDIR)

exp_openmp: init makebin $(EXPDIR)exp_openmp.c
	$(CC) $(CFLAGS) $(COPTIM) $(EXPDIR)exp_openmp.c \
		$(CLIBS) $(PCG_FLAGS) -o $(BINDIR)$@

run_exp_openmp: makedat exp_openmp
	$(BINDIR)exp_openmp
	$(MV) *.dat $(DATDIR)

exp_overhead: init makebin $(EXPDIR)exp_overhead.c
	$(CC) $(CFLAGS) $(COPTIM) $(EXPDIR)exp_overhead.c \
		$(CLIBS) $(PCG_FLAGS) -o $(BINDIR)$@

run_exp_overhead: makedat exp_overhead
	$(BINDIR)exp_overhead
	$(MV) *.dat $(DATDIR)

# MATLAB experiments
run_exp_matlab: EXPSTRING="addpath('$(INCDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(INCDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps; \
		exit;"

run_exp_matlab: makedat mexmat
	$(MATLAB) -r $(EXPSTRING)

run_exp_matlab_extra: EXPSTRING="addpath('$(INCDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(INCDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps_extra; \
		exit;"

run_exp_matlab_extra: makedat mexmat
	$(MATLAB) -r $(EXPSTRING)

clean_experiments:
	$(RM) -f $(BINDIR)exp_*




.PHONY: cleanall clean cleandep cleantest cleancoverage cleandocs cleandat
cleanall: clean cleandep cleantest cleancoverage cleandocs cleandat

clean:
	$(RM) $(BINDIR)*

cleandep:
	cd $(INCDIR)pcg-c; make clean

cleantest:
	$(RM) $(TESTDIR)cpfloat_test $(TESTDIR)*.c $(TESTDIR)*.log

cleancoverage:
	$(RM) cpfloat_test.c cpfloat_test.log *.gcno *.gcda *.gcov

cleandocs:
	$(RM) -r $(DOCDIR)Doxyfile $(DOCDIR)xml
	$(RM) -r $(DOCDIR)html $(DOCDIR)source/cpfloat

cleandat:
	$(RM) $(DATDIR)*

.PHONY: update-spdx
update-spdx:
	$(UTILDIR)generate_spdx.sh > license.spdx

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
