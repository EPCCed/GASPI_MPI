# 
# Makefile
# Copyright 2015 The University of Edinburgh
#
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at 
#
# http://www.apache.org/licenses/LICENSE-2.0 
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License.  
#

ifneq ("","$(wildcard ../Makefile.mk)") #use global settings

include ../Makefile.mk

else #default settings

#CPU:
CC=gcc
CFLAGS=-O2 -Wall -fopenmp

#GPU
#CC=nvcc
#CFLAGS=-O2 -arch=sm_35 -x cu -dc

endif

AR = ar -cru

.SUFFIXES:
.SUFFIXES: .c  .o

TARGETSRC=targetDP_C

ifeq ($(CC),nvcc)
TARGETSRC=targetDP_CUDA
endif

default:
	$(CC) $(CFLAGS) -c $(TARGETSRC).c
	$(AR) libtarget.a $(TARGETSRC).o

clean:
	rm -f *.o *.a
