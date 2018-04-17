
#WARNING: You should probably change the compute architecture for your GPU in BASIC_LINE or here in ARCH
#The target CUDA compute capability	
ARCH=52

#Change for your system, HG needs to be already compiled
HYDROGRID_SRC=HydroGrid/src/

#Uncomment to compile in double precision mode

#DOUBLE_PRECISION=-DDOUBLE_PRECISION

UAMMD_SRC=uammd/src

CPU= -O3 -funroll-loops -ffinite-math-only -fno-signaling-nans -fno-math-errno -fno-signed-zeros -frename-registers -march=native -fPIC

CXX=g++
BASIC_LINE= nvcc -O3 $(DOUBLE_PRECISION) -I  $(UAMMD_SRC) -I $(UAMMD_SRC)/third_party -ccbin="$(CXX)" -Xcompiler="$(CPU)" -gencode arch=compute_$(ARCH),code=sm_$(ARCH) -x cu -std=c++11 --expt-relaxed-constexpr

all: hgMC

hgMC:
	@if ! ls -d $(HYDROGRID_SRC) >/dev/null 2>&1; \
	then \
	echo "ERROR: Could not find Hydrogrid!, please compile it and tell me where it is in the Makefile!"; exit 1; \
	fi

	@if ! ls -d uammd >/dev/null 2>&1; \
	then \
	git clone https://github.com/pabloibannez/uammd; \
	fi


	$(BASIC_LINE) -I$(HYDROGRID_SRC) DiffusionWall_MC.cu  -Xlinker="-rpath=\$$ORIGIN" -L$(HYDROGRID_SRC)   -lCallHydroGrid -o mc
	cp $(HYDROGRID_SRC)/libCallHydroGrid.so .	

clean:
	rm -f libCall*so mc
