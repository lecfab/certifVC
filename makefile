# certifVC: Quality certification for Vertex Cover
# Copyright (C) 2023 Fabrice Lécuyer (fabrice.lecuyer@lip6.fr)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 or later.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
# for more details: https://www.gnu.org/licenses/

CC=g++ -std=c++14 #-fconcepts

DEBUG?=0 # debug option

ifeq ($(DEBUG), 1)
	CFLAGS=-Og -Wextra -g3 -D DEBUG  -fopenmp
	o=debug.o
else
	# -g for debugging
	CFLAGS=-Ofast -g -fopenmp # -Wall https://cpluspluspedia.com/fr/tutorial/4708/compiler-et-construire # https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
	o=o
endif

EXEC=certifVC
MF=makefile # recompile when Makefile has been modified

UTILS_H =	utils/tools.h \
			utils/inout.h \
			utils/edgelist.h \
			utils/adjlist.h \
			utils/heap.h
			# utils/CLI11.h
ALGOS_H = 	algo/algo_vertexcover.h

UTILS_O =	$(UTILS_H:.h=.$(o))
ALGOS_O = 	$(ALGOS_H:.h=.$(o))

all: $(EXEC)

certifVC: certifVC.$(o) $(UTILS_O) $(ALGOS_O) $(ORDERS_O)
	$(CC) $^ $(CFLAGS) -o $@

certifVC.$(o): $(UTILS_H) $(ALGOS_H)

algo/algo_%.$(o): $(UTILS_H)
order/order_%.$(o): $(UTILS_H)

# utils/CLI11.$(o): utils/CLI11.h
# 	$(CC) -c $< $(CFLAGS) -o $@

%.debug.o: %.cpp $(MF)
	$(CC) -c $< $(CFLAGS) -o $@

%.o: %.cpp $(MF)
	@$(CC) -c $< $(CFLAGS) -o $@

.PHONY: clean mrproper

clean:
	rm -rf */*.$(o) *.$(o)

mrproper: clean
	rm -rf $(EXEC)
