#ARCH ?= #k100#geocluster #gpupc1 #D
#MPI_ON ?= 1

ifdef MPI_ON
GCC  ?= mpic++ -std=c++11
else
GCC  ?= g++ -std=c++11
endif

ifeq ($(ARCH),k100)
NVCC := /common/cuda-6.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_20
else ifeq ($(ARCH),geocluster)
NVCC := nvcc -ccbin $(GCC) -O3
GENCODE_SM := -arch=sm_50
else ifeq ($(ARCH),gpupc1)
NVCC := nvcc -ccbin $(GCC) -O3
GENCODE_SM := -arch=sm_35
else ifeq ($(ARCH),D)
NVCC := /home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 -g 
GENCODE_SM := -arch=sm_35
else ifeq ($(ARCH),photon)
NVCC := /mnt/D/home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_50
else ifeq ($(ARCH),supermic)
NVCC := /usr/local/cuda/bin/nvcc -ccbin $(GCC) -O3 -g -I/usr/local/cuda/include/
GENCODE_SM := -arch=sm_61
else ifeq ($(ARCH),plasma)
NVCC := nvcc -ccbin $(GCC) -O3 -g -I/usr/local/cuda/include/
GENCODE_SM := -arch=sm_52
else ifeq ($(ARCH),ion)
NVCC := /mnt/D/home/zakirov/cuda-7.5/bin/nvcc -ccbin $(GCC) -O3 
GENCODE_SM := -arch=sm_50
else ifeq ($(ARCH),kiae)
NVCC := nvcc -ccbin $(GCC) -O3
GENCODE_SM := -arch=sm_37
NOGL=1
else
NVCC := nvcc -ccbin $(GCC) -O3
GENCODE_SM := -arch=sm_30
endif 
#ALL_DEFS := SET_FROM_PYTHON NX NY NZ LNX LNY LNZ LSCALE GAUSS_PROFILE ENABLE_TEMPERATURE GRAVITY MARANGONI MIDAS BUBBLES USE_USUAL_BOUNCE_BACK VARIED_TSTEP USE_FLOAT USE_DOUBLE EXT_HEAT SHANCHEN SOLVE_GAS USE_HOST_MEM DISPLAY_FULLX DISPLAY_FULLY NOGL
#CDEFS := $(foreach f, $(ALL_DEFS), $(if $($f),-D$f=$($f)))
ALL_DEFS := SEM_WAIT CFSTEPS MORTON_ZIP
CDEFS := $(foreach f, $(ALL_DEFS), $(if $($f),-D$f=$($f)))

# internal flags
NVCCFLAGS   := -Xptxas="-v"   -Xcudafe "--diag_suppress=declared_but_not_referenced"
CCFLAGS     := -Ofast -fopenmp -fPIC $(CDEFS) 
NVCCLDFLAGS :=
LDFLAGS     := -L./
ifeq ($(ARCH),supermic)
LDFLAGS     := -L./ -L/usr/local/cuda/lib64
endif 
ifeq ($(ARCH),kiae)
CCFLAGS     := -O3 -fopenmp -fPIC $(CDEFS) 
endif

# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_NVCCLDFLAGS ?=
EXTRA_LDFLAGS     ?=
EXTRA_CCFLAGS     ?= #-std=c++11

ifeq ($(ARCH),k100)
INCLUDES  := -I./
LIBRARIES := -lcudart -lglut -lGL -lcufft -lpng -lgomp -lpthread
else ifeq ($(ARCH),geocluster)
INCLUDES  := -I./
LIBRARIES := -lcudart -lglut -lGL -lcufft -lpng -lgomp -lpthread 
else ifeq ($(ARCH),kiae)
INCLUDES  := -I./
LIBRARIES := -lcudart -lglut -lGL -lcufft       -lgomp -lpthread 
else
INCLUDES  := -I./ 
LIBRARIES := -lcudart -lglut -lGL -lcufft -lpng -lgomp -lpthread
endif
ifdef MPI_ON
LIBRARIES := -lmpi $(LIBRARIES)
endif

################################################################################
GENCODE_SM20  := #-gencode arch=compute_20,code=sm_21
GENCODE_SM30  := #-gencode arch=compute_30,code=sm_30
GENCODE_SM35  := #-gencode arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_SM50  := #-arch=sm_20
GENCODE_FLAGS := $(GENCODE_SM50) $(GENCODE_SM35) $(GENCODE_SM30) $(GENCODE_SM20) $(GENCODE_SM)
ALL_CCFLAGS   := --compiler-options="$(CCFLAGS) $(EXTRA_CCFLAGS)" 
ALL_LDFLAGS   := --linker-options="$(LDFLAGS) $(EXTRA_LDFLAGS)"
################################################################################

# Target rules
all: build

build: lbm.x

im3D.o: im3D.cu im3D.hpp cuda_math.h fpal.h im2D.h
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -DSURF -o $@ -dc $<

depend: .depend

err.dep: err.cu
	rm -f err.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >err.dep ;

update.dep: update.cu
	rm -f update.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >update.dep;

#steps.dep: steps.cu
#	rm -f steps.dep
#	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >steps.dep;
lbm_pull.dep: lbm_pull.cu
	rm -f lbm_pull.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >lbm_pull.dep;

init.dep: init.cu
	rm -f init.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >init.dep;

drop.dep: drop.cu
	rm -f drop.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >drop.dep;

diagn.dep: diagn.cu
	rm -f diagn.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >diagn.dep;


lbpack.dep: lbpack.cu
	rm -f lbpack.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ >lbpack.dep;

lbm_main.dep: lbm.cu
	rm -f lbm_main.dep
	$(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -M $^ -MT kerlbm.o >lbm_main.dep;

include err.dep
include lbm_main.dep
include update.dep
#include steps.dep
include lbm_pull.dep
include init.dep
include lbpack.dep

kerlbm.o: lbm.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
err.o: err.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
update.o: update.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
#steps.o: steps.cu
#	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) -maxrregcount 128 $(GENCODE_FLAGS) -o $@ -dc $<
lbm_pull.o: lbm_pull.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
init.o: init.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
drop.o: drop.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
diagn.o: diagn.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
lbpack.o: lbpack.cu
	$(EXEC) $(NVCC) $(ALL_CCFLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) -o $@ -dc $<
	
lbm.x: kerlbm.o im3D.o err.o update.o lbm_pull.o init.o drop.o diagn.o lbpack.o
	$(EXEC) $(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) $(LDFLAGS) -o $@ $+ $(LIBRARIES)

run: build
	$(EXEC) ./lbm.x

clean:
	$(EXEC) rm -f lbm.x kerlbm.o im3D.o err.o update.o lbm_pull.o init.o drop.o diagn.o lbpack.o *.dep *_wrap.cxx

clobber: clean

#Xcudafe --diag_suppress warnings list
#http://www.ssl.berkeley.edu/~jimm/grizzly_docs/SSL/opt/intel/cc/9.0/lib/locale/en_US/mcpcom.msg
