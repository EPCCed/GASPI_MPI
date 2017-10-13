CFLAGS += -Wshadow
CFLAGS += -O3
CFLAGS += -openmp

ifeq ($(GASPI_HOME),)
  $(error GASPI_HOME is not set)
endif


ifeq ($(MPI_HOME),)
  $(error MPI_HOME is not set)
endif

CC=$(MPI_HOME)/bin64/mpicxx

INCLUDE_DIR += $(MPI_HOME)/include64
INCLUDE_DIR += $(GASPI_HOME)/include 
INCLUDE_DIR += ..

MPI_LIB_DIR = $(MPI_HOME)/lib64
LIBRARY_DIR += $(MPI_LIB_DIR)
LIBRARY_DIR += $(GASPI_HOME)/lib64

LDFLAGS += $(addprefix -L,$(LIBRARY_DIR))

ifeq "$(CC)" "icpc"
  LDFLAGS += -openmp-link static
endif

CFLAGS += $(addprefix -I,$(INCLUDE_DIR)) 
CFLAGS += -fopenmp

LIB += ibverbs
LIB += m
STATIC_LIBS = $(MPI_LIB_DIR)/libmpi.a $(GASPI_HOME)/lib64/libGPI2.a

EXECUTABLES = jacobi_mpi \
              jacobi_mixed \
              use_gaspi_mem_in_mpi_comm \
              use_app_provided_mem_for_segments \
              use_app_provided_mem_for_segments_lib

all: $(EXECUTABLES)
    
jacobi_mpi: jacobi_mpi.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB)) $(STATIC_LIBS)
	
jacobi_mixed: jacobi_mixed.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB)) $(STATIC_LIBS)
	
use_gaspi_mem_in_mpi_comm: use_gaspi_mem_in_mpi_comm.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB)) $(STATIC_LIBS)
	
use_app_provided_mem_for_segments: use_app_provided_mem_for_segments.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB)) $(STATIC_LIBS)
	
use_app_provided_mem_for_segments_lib: use_app_provided_mem_for_segments_lib.cpp jacobi_gpi_alg.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB)) $(STATIC_LIBS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

jacobi_gpi_alg.a: jacobi_gpi_alg.o
	ar rcs jacobi_gpi_alg.a jacobi_gpi_alg.o
	ranlib jacobi_gpi_alg.a
	
.PHONY: clean objclean

objclean:
	rm -f *.o *.a

clean: objclean
	rm -f $(EXECUTABLES)
