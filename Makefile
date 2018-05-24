
#WARNING: You should probably change the compute architecture for your GPU in BASIC_LINE or here in ARCH
#If not set, the Makefile will compile for the currently available GPUs in the system
#The target CUDA compute capability	
ARCH=

#Change for your system, HG needs to be already compiled
HYDROGRID_SRC=HydroGrid/src/

#Uncomment to compile in double precision mode

#DOUBLE_PRECISION=-DDOUBLE_PRECISION

UAMMD_SRC=uammd/src

ifeq ($(ARCH),)
GENCODE_FLAGS:=$(shell printf '\#include<cstdio>\n int main(){int nD;cudaGetDeviceCount(&nD);for(int i=0;i<nD;i++){cudaDeviceProp dp;cudaGetDeviceProperties(&dp, i);int cp=dp.major*10+dp.minor;std::printf("%%d\\n",cp);} return 0;}' | nvcc -Wno-deprecated-gpu-targets -x cu - --run | sort -g -k1 | uniq | awk '{print "-gencode arch=compute_"$$1",code=sm_"$$1}')
else
GENCODE_FLAGS=-gencode arch=compute_$(ARCH),code=sm_$(ARCH)
endif


CPU= -O3 -g -funroll-loops -ffinite-math-only -fno-signaling-nans -fno-math-errno -fno-signed-zeros -frename-registers -march=native -fPIC

CXX=g++
BASIC_LINE= nvcc -g -O3 $(DOUBLE_PRECISION) -I  $(UAMMD_SRC) -I $(UAMMD_SRC)/third_party -ccbin="$(CXX)" -Xcompiler="$(CPU)" $(GENCODE_FLAGS) -x cu -std=c++11 --expt-relaxed-constexpr

all: hgMC

hgMC:
	@if ! ls $(HYDROGRID_SRC)/libCallHydroGrid.so >/dev/null 2>&1; \
	then \
	echo "ERROR: Could not find Hydrogrid!, please compile it and tell me where it is in the Makefile!"; exit 1; \
	fi

	@if ! ls -d uammd >/dev/null 2>&1; \
	then \
	git clone --depth=1 https://github.com/pabloibannez/uammd; \
	fi


	$(BASIC_LINE) -I$(HYDROGRID_SRC) DiffusionWall_MC.cu  -Xlinker="-rpath=\$$ORIGIN" -L$(HYDROGRID_SRC)   -lCallHydroGrid -o mc
	cp $(HYDROGRID_SRC)/libCallHydroGrid.so .	

clean:
	rm -f libCall*so mc
