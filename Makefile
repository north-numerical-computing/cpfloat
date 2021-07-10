# SPDX-FileCopyrightText: 2020 Massimiliano Fasi and Mantas Mikaitis
# SPDX-License-Identifier: LGPL-2.1-or-later

ROOTDIR=$(shell pwd)
INCDIR=$(ROOTDIR)/include/
SRCDIR=$(ROOTDIR)/src/
DOCDIR=$(ROOTDIR)/doc/
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
GIT:=git
MATLAB:=matlab -nodesktop -nosplash
OCTAVE:=octave

WFLAGS=-Wall -Wextra -pedantic
ARCHFLAGS=-mfma -march=native
CFLAGS=$(WFLAGS) $(ARCHFLAGS) -I $(SRCDIR) \
	-I /usr/local/include -L /usr/local/lib \
	-include $(PCG_HEADER)
CXXFLAGS=$(WFLAGS) $(ARCHFLAGS) -std=c++11 -I $(INCDIR) -I $(INCDIR)FloatX/src/
COPTIM=-O3 -mfma
CCOVFLAGS=-Og -g -fprofile-arcs -ftest-coverage
CLIBS=-lm -fopenmp

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

makebin:
	@if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi

autotune: init makebin $(SRCDIR)cpfloat_autotune.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) \
		-o $(BINDIR)cpfloat_autotune \
		$(SRCDIR)cpfloat_autotune.c
	$(BINDIR)cpfloat_autotune
	$(MV) cpfloat_threshold_*.h $(SRCDIR)

test: ctest mtest

ctestsrc: $(TESTDIR)cpfloat_test.ts
	$(CHECKMK) clean_mode=1 $^ > $(TESTDIR)cpfloat_test.c

ctest: init makebin ctestsrc
	$(CC) $(CFLAGS) $(COPTIM) \
		-o $(BINDIR)cpfloat_test $(TESTDIR)cpfloat_test.c \
		-lcheck -lsubunit -lrt -lm -lpthread
	$(BINDIR)cpfloat_test
	$(MV) cpfloat_test.log $(TESTDIR)

MEXSTRING="cd $(MEXDIR); \
	retval = cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
		'pcgpath', '$(PCG_HEADER)', \
		'compilerpath', '$(CC)'); \
	if retval \
		rehash(); \
		cpfloat_autotune('cpfloatdir', '$(SRCDIR)'); \
		cpfloat_compile('cpfloatdir', '$(SRCDIR)', \
			'pcgpath', '$(PCG_HEADER)', \
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
	$(DOXYGEN) doc/Doxyfile

coverage: init ctestsrc
	$(CC) $(CFLAGS) $(CCOVFLAGS) \
		-o $(TESTDIR)cpfloat_test -lcheck $(TESTDIR)cpfloat_test.c
	$(TESTDIR)cpfloat_test
	$(CP) $(TESTDIR)cpfloat_test.c .
	$(CCOV) cpfloat_test.c

example: init makebin $(EXAMPLEDIR)example_manuscript.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) \
		-o $(BINDIR)example_manuscript \
		$(EXAMPLEDIR)example_manuscript.c





.PHONY:
experiments: run_exp_ccomp run_exp_openmp run_exp_overhead run_exp_matlab

exp_comp_cpfloat: init makebin $(EXPDIR)exp_comp_cpfloat.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) -I $(SRCDIR) \
		-o $(BINDIR)$@ $(EXPDIR)exp_comp_cpfloat.c

exp_comp_mpfr: init makebin $(EXPDIR)exp_comp_mpfr.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) -lmpfr -I $(SRCDIR) \
		-o $(BINDIR)$@ $(EXPDIR)exp_comp_mpfr.c

exp_comp_floatx: init makebin $(EXPDIR)exp_comp_floatx.cpp
	$(CXX) $(CXXFLAGS) $(COPTIM) -I $(SRCDIR)/FloatX/src \
		-o $(BINDIR)$@ $(EXPDIR)exp_comp_floatx.cpp

run_exp_ccomp: exp_comp_cpfloat exp_comp_mpfr exp_comp_floatx
	$(BINDIR)exp_comp_cpfloat
	$(BINDIR)exp_comp_mpfr
	$(BINDIR)exp_comp_floatx
	$(MV) *.dat $(DATDIR)

exp_openmp: init makebin $(EXPDIR)exp_openmp.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) \
		-o $(BINDIR)$@ $(EXPDIR)exp_openmp.c

run_exp_openmp: exp_openmp
	$(BINDIR)exp_openmp
	$(MV) *.dat $(DATDIR)

exp_overhead: init makebin $(EXPDIR)exp_overhead.c
	$(CC) $(CFLAGS) $(COPTIM) $(CLIBS) \
		-o $(BINDIR)$@ $(EXPDIR)exp_overhead.c

run_exp_overhead: exp_overhead
	$(BINDIR)exp_overhead
	$(MV) *.dat $(DATDIR)

run_exp_matlab: EXPSTRING="addpath('$(INCDIR)chop'); \
		addpath('$(BINDIR)'); \
		addpath(genpath('$(INCDIR)floatp/')); \
		cd $(EXPDIR); \
		datdir = '$(DATDIR)'; \
		run_exps; \
		exit;"

run_exp_matlab: mexmat
	$(MATLAB) -r $(EXPSTRING)

clean_experiments:
	$(RM) -f $(BINDIR)exp_*




.PHONY: cleanall clean cleantest cleancoverage cleandoc cleandat
cleanall: clean cleantest cleancoverage cleandoc cleandat

clean:
	$(RM) $(BINDIR)*

cleantest:
	$(RM) $(TESTDIR)cpfloat_test $(TESTDIR)*.c $(TESTDIR)*.log

cleancoverage:
	$(RM) cpfloat_test.c cpfloat_test.log *.gcno *.gcda *.gcov

cleandoc:
	$(RM) -r $(DOCDIR)html/* $(DOCDIR)latex/*

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
