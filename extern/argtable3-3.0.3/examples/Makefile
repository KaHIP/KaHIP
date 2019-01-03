################################################################################
# This file is part of the argtable3 library.
#
# Copyright (C) 2013-2014 Tom G. Huang
# <tomghuang@gmail.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of STEWART HEITMANN nor the  names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL STEWART HEITMANN BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################

CC = gcc
RM = rm -f
MKDIR = mkdir -p
RMDIR = rm -rf

ifeq ($(wildcard ../argtable3.c),)
    ARGTABLE_SRC = ../src/argtable3.c \
		../src/arg_date.c \
		../src/arg_dbl.c \
		../src/arg_end.c \
		../src/arg_file.c \
		../src/arg_int.c \
		../src/arg_lit.c \
		../src/arg_rem.c \
		../src/arg_rex.c \
		../src/arg_str.c \
		../src/arg_utils.c \
		../src/getopt_long.c
    ARGTABLE_HEADER = ../src/argtable3.h
	ARGTABLE_INCLUDE_DIR = ../src
else
    ARGTABLE_SRC = ../argtable3.c
    ARGTABLE_HEADER = ../argtable3.h
	ARGTABLE_INCLUDE_DIR = ..
endif

.PHONY: all cleanall
all: echo ls multisyntax mv myprog myprog_C89 uname testargtable3

echo: echo.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) echo.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

ls: ls.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) ls.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

multisyntax: multisyntax.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) multisyntax.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

mv: mv.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) mv.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

myprog: myprog.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) myprog.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

myprog_C89: myprog_C89.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) myprog_C89.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

testargtable3: testargtable3.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) testargtable3.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@

uname: uname.c $(ARGTABLE_SRC) $(ARGTABLE_HEADER)
	$(MKDIR) bin
	$(CC) uname.c $(ARGTABLE_SRC) -I $(ARGTABLE_INCLUDE_DIR) -o bin/$@


cleanall:
	$(RMDIR) bin

