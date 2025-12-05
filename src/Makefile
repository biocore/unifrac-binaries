.PHONY: all main api api_acc test test_binaries install install_lib clean rapi_test

PLATFORM := $(shell uname -s)
COMPILER := $(shell (h5c++ -v 2>&1) | tr A-Z a-z )

ifeq ($(BUILD_VARIANT),)
	SSU = ssu
	FPD = faithpd
else
	SSU = ssu_$(BUILD_VARIANT)
	FPD = faithpd_$(BUILD_VARIANT)
endif

ifdef DEBUG
	ifneq (,$(findstring pgi,$(COMPILER)))
		OPT = -g
	else
		OPT = -O0 -DDEBUG=1 --debug -g -ggdb
	endif
else
	ifneq (,$(findstring pgi,$(COMPILER)))
		OPT = -fast
	else
	  ifneq (,$(findstring gcc,$(COMPILER)))
		OPT = -O3
		TGTFLAGS = -fwhole-program
	  else
		OPT = -O3
	  endif
	endif
endif

ifeq ($(PREFIX),)
	PREFIX := $(CONDA_PREFIX)
endif

SKBBLIB=-lskbb

ifeq ($(PLATFORM),Darwin)
	LDDFLAGS = -dynamiclib -install_name @rpath/lib$(SSU).so
	CPPFLAGS+= -D_LIBCPP_DISABLE_AVAILABILITY
	NOGPU := 1
else
	LDDFLAGS = -shared
	ARCH := $(shell h5c++ -dumpmachine)
endif

EXEFLAGS =

ifneq (,$(findstring pgi,$(COMPILER)))
	MPFLAG =  -mp
else
	MPFLAG = -fopenmp
endif

LDDFLAGS += $(MPFLAG)
CPPFLAGS += $(MPFLAG)

UNIFRAC_FILES = unifrac_internal.o unifrac_task_cpu.o unifrac_accapi_cpu.o unifrac_cmp_cpu.o 
UNIFRAC_ACC_SHLIBS =

ifeq ($(NV_CXX),)
   NV_COMPILER := 
else
   NV_COMPILER := $(shell ($(NV_CXX) -v 2>&1) | tr A-Z a-z )
endif

ifeq ($(AMD_CXX),)
   AMD_COMPILER := 
else
   AMD_COMPILER := $(shell ($(AMD_CXX) -v 2>&1) | tr A-Z a-z )
endif

ifndef NOGPU
   # NVIDIA GPUs
   UNIFRAC_FILES += unifrac_task_acc_nv.o unifrac_accapi_acc_nv.o unifrac_cmp_acc_nv.o
   CPPFLAGS += -DUNIFRAC_ENABLE_ACC_NV=1
   ifeq ($(BUILD_NV_OFFLOAD),ompgpu)
	ifneq (,$(findstring pgi,$(NV_COMPILER)))
		UNIFRAC_ACC_SHLIBS += libssu_acc_nv.so
		NV_CPPFLAGS += -mp=gpu -DOMPGPU=1 -noacc
	        ifneq ($(PERFORMING_CONDA_BUILD),False)
	            NV_CPPFLAGS += -gpu=ccall
		else
	            NV_CPPFLAGS += -gpu=ccnative
		    # optional info
		    NV_CPPFLAGS += -Minfo=accel
                endif
		NV_LDFLAGS += -shared -mp=gpu -Bstatic_pgi
	endif
   else 
	ifneq (,$(findstring pgi,$(NV_COMPILER)))
		UNIFRAC_ACC_SHLIBS += libssu_acc_nv.so
		NV_CPPFLAGS += -acc
	        ifneq ($(PERFORMING_CONDA_BUILD),False)
	            NV_CPPFLAGS += -gpu=ccall
		else
	            NV_CPPFLAGS += -gpu=ccnative
		    # optional info
		    NV_CPPFLAGS += -Minfo=accel
                endif
		NV_LDFLAGS += -shared -acc -Bstatic_pgi
	endif
   endif
   ifneq (,$(findstring x86_64,$(ARCH)))
     NV_CPPFLAGS += -tp=px
   else
     NV_CPPFLAGS += -tp=neoverse-n1
   endif
   NV_CPPFLAGS += -Mfpapprox -Mfprelaxed -Mflushz
   # AMD GPUs
   UNIFRAC_FILES += unifrac_task_acc_amd.o unifrac_accapi_acc_amd.o unifrac_cmp_acc_amd.o 
   CPPFLAGS += -DUNIFRAC_ENABLE_ACC_AMD=1
   ifneq (,$(findstring clang,$(AMD_COMPILER)))
		UNIFRAC_ACC_SHLIBS += libssu_acc_amd.so
		AMD_CPPFLAGS += -fopenmp -fopenmp-offload-mandatory -DOMPGPU=1
	        ifneq ($(PERFORMING_CONDA_BUILD),False)
	            AMD_CPPFLAGS += --offload-arch=gfx1100,gfx1101,gfx1102,gfx1103,gfx1030,gfx1031,gfx90a,gfx942
		else
	            AMD_CPPFLAGS += --offload-arch=native
                endif
		AMD_LDFLAGS += -shared -fopenmp --offload-arch=gfx1100,gfx1101,gfx1102,gfx1103,gfx1030,gfx1031,gfx90a,gfx942
   endif
endif

ifeq ($(BUILD_FULL_OPTIMIZATION),False)
    ifneq (,$(findstring pgi,$(COMPILER)))
	CPPFLAGS += -tp=px
    else
	ifeq ($(PLATFORM),Darwin)
         	# no special tuning for the MacOS at this point
	  	CPPFLAGS += 
	else
		ifneq (,$(findstring x86_64,$(ARCH)))
			CPPFLAGS += -march=x86-64 -mtune=generic
		else
			# target server-class ARM processors
			CPPFLAGS += -march=armv8.2-a -mtune=neoverse-n1
		endif
	endif
    endif
else
  ifneq ($(BUILD_FULL_OPTIMIZATION),)
    ifneq (,$(findstring pgi,$(COMPILER)))
	CPPFLAGS += -tp=$(BUILD_FULL_OPTIMIZATION)
    else
	CPPFLAGS += -march=$(BUILD_FULL_OPTIMIZATION)
        ifneq ($(BUILD_TUNE_OPTIMIZATION),)
	  CPPFLAGS += -mtune=$(BUILD_TUNE_OPTIMIZATION)
        endif
    endif
  else
    # no preference from the user
    # change default only if not performing CONDA_BUILD
    ifeq ($(PERFORMING_CONDA_BUILD),False)
	ifeq ($(PLATFORM),Darwin)
         	# no special tuning for the MacOS at this point
	  	CPPFLAGS += 
	else
		CPPFLAGS += -march=native -mtune=native
	endif
    endif
  endif
endif

ifeq (,$(findstring pgi,$(COMPILER)))
    ifeq ($(PERFORMING_CONDA_BUILD),False)
	# basically, not gcc
	CPPFLAGS += -Wextra -Wno-unused-parameter
    endif
endif


CPPFLAGS += -Wall  -std=c++17 -pedantic -I. $(OPT) -fPIC -L$(PREFIX)/lib
NV_CPPFLAGS += -std=c++17 -I. $(OPT) -fPIC -L$(PREFIX)/lib
AMD_CPPFLAGS += -std=c++17 -I. $(OPT) -fPIC -L$(PREFIX)/lib

all: api install_lib main install

main: $(SSU) $(FPD)

api: lib$(SSU).so

api_acc: $(UNIFRAC_ACC_SHLIBS)

test_binaries: test_su test_ska test_api test_su_api

# first test the internal structures
test_su: test_su.cpp api.hpp tree.o tsv.o biom.o biom_inmem.o biom_subsampled.o unifrac.o skbio_alt.o api.o $(UNIFRAC_FILES)
	h5c++ $(CPPFLAGS) $(EXEFLAGS) test_su.cpp -o $@ tree.o biom.o biom_inmem.o biom_subsampled.o tsv.o $(UNIFRAC_FILES) unifrac.o skbio_alt.o api.o -llz4 $(SKBBLIB) -lpthread
# then exercise the optimized code through the API interface
test_su_api: test_su.cpp test_helper.hpp api.hpp $(PREFIX)/lib/lib$(SSU).so
	h5c++ $(CPPFLAGS) $(EXEFLAGS) -DAPI_ONLY=1 test_su.cpp -o $@ -l$(SSU)

# Create an archive, so we do not need dependencies for functions we do not use
test_ska_deps.a: tree.o tsv.o biom.o biom_inmem.o biom_subsampled.o skbio_alt.o api.o $(UNIFRAC_FILES)
	rm -f $@
	ar -r $@ $^
 
# Not invoking the unifrac acc code, so we do not need to use the shlib
test_ska: test_ska.cpp skbio_alt.hpp api.hpp test_helper.hpp biom.hpp test_ska_deps.a
	h5c++ $(CPPFLAGS) $(EXEFLAGS) test_ska.cpp -o test_ska test_ska_deps.a -llz4 $(SKBBLIB) -lpthread

test_api: test_api.cpp api.hpp test_helper.hpp biom.hpp $(PREFIX)/lib/lib$(SSU).so
	h5c++ $(CPPFLAGS) $(EXEFLAGS) test_api.cpp -o $@ -l$(SSU)

$(SSU): su.cpp api.hpp biom.hpp $(PREFIX)/lib/lib$(SSU).so
	h5c++ $(CPPFLAGS) $(EXEFLAGS) su.cpp -o $@ -l$(SSU)

$(FPD): faithpd.cpp api.hpp $(PREFIX)/lib/lib$(SSU).so
	h5c++ $(CPPFLAGS) $(EXEFLAGS) faithpd.cpp -o $@ -l$(SSU)

lib$(SSU).so: tree.o biom.o biom_inmem.o biom_subsampled.o tsv.o unifrac.o cmd.o skbio_alt.o api.o $(UNIFRAC_FILES)
	h5c++ $(LDDFLAGS) -o lib$(SSU).so tree.o biom.o biom_inmem.o biom_subsampled.o tsv.o $(UNIFRAC_FILES) unifrac.o cmd.o skbio_alt.o api.o -lc -llz4 $(SKBBLIB) -L$(PREFIX)/lib -noshlib -lhdf5_cpp -lhdf5_hl_cpp -lhdf5_hl -lhdf5

api.o: api.cpp api.hpp unifrac.hpp skbio_alt.hpp biom.hpp biom_inmem.hpp biom_subsampled.hpp tree.hpp tsv.hpp
	h5c++ $(CPPFLAGS) api.cpp -c -o api.o -fPIC


# Use the same compiler for CPU-based code
# So no need for futher levels of indirection
unifrac_accapi_cpu.cpp: unifrac_accapi.hpp generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py cpu direct > $@
unifrac_task_noclass_cpu.cpp: unifrac_task_impl.hpp generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py cpu direct > $@

##

unifrac_accapi_dyn_acc_nv.h: unifrac_accapi_impl.hpp generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_nv api_h > $@
unifrac_accapi_dyn_acc_nv.cpp: unifrac_accapi_impl.hpp unifrac_accapi_dyn_acc_nv.h generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_nv api > $@
unifrac_accapi_acc_nv.cpp: unifrac_accapi_impl.hpp unifrac_accapi_dyn_acc_nv.h generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_nv indirect > $@

unifrac_task_api_acc_nv.h: unifrac_task_impl.hpp generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_nv api_h > $@
unifrac_task_api_acc_nv.cpp: unifrac_task_impl.hpp unifrac_task_api_acc_nv.h generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_nv api > $@
unifrac_task_noclass_acc_nv.cpp: unifrac_task_impl.hpp unifrac_task_api_acc_nv.h generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_nv indirect > $@

##

unifrac_accapi_dyn_acc_amd.h: unifrac_accapi_impl.hpp generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_amd api_h > $@
unifrac_accapi_dyn_acc_amd.cpp: unifrac_accapi_impl.hpp unifrac_accapi_dyn_acc_amd.h generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_amd api > $@
unifrac_accapi_acc_amd.cpp: unifrac_accapi_impl.hpp unifrac_accapi_dyn_acc_amd.h generate_unifrac_accapi.py
	python3 generate_unifrac_accapi.py acc_amd indirect > $@

unifrac_task_api_acc_amd.h: unifrac_task_impl.hpp generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_amd api_h > $@
unifrac_task_api_acc_amd.cpp: unifrac_task_impl.hpp unifrac_task_api_acc_amd.h generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_amd api > $@
unifrac_task_noclass_acc_amd.cpp: unifrac_task_impl.hpp unifrac_task_api_acc_amd.h generate_unifrac_task_noclass.py
	python3 generate_unifrac_task_noclass.py acc_amd indirect > $@

##

unifrac_accapi_cpu.o: unifrac_accapi_cpu.cpp unifrac_accapi.hpp unifrac_accapi_impl.hpp
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@
unifrac_task_cpu.o: unifrac_task_noclass_cpu.cpp unifrac_task_noclass.hpp unifrac_task_impl.hpp
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@


# no separate task_api 

# create symlink to make R (and any autodeduct system) happy
unifrac_accapi.cpp: unifrac_accapi_cpu.cpp
	ln -s unifrac_accapi_cpu.cpp unifrac_accapi.cpp
unifrac_task.cpp: unifrac_task_noclass_cpu.cpp
	ln -s unifrac_task_noclass_cpu.cpp unifrac_task.cpp

##

unifrac_accapi_acc_nv.o: unifrac_accapi_acc_nv.cpp unifrac_accapi.hpp unifrac_accapi_dyn_acc_nv.h
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@
unifrac_accapi_dyn_acc_nv.o: unifrac_accapi_dyn_acc_nv.cpp unifrac_accapi_dyn_acc_nv.h unifrac_task_impl.hpp
	$(NV_CXX) $(NV_CPPFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@

unifrac_task_acc_nv.o: unifrac_task_noclass_acc_nv.cpp unifrac_task_noclass.hpp unifrac_task_api_acc_nv.h
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@
unifrac_task_api_acc_nv.o: unifrac_task_api_acc_nv.cpp unifrac_task_api_acc_nv.h unifrac_task_impl.hpp
	$(NV_CXX) $(NV_CPPFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@

libssu_acc_nv.so: unifrac_task_api_acc_nv.o unifrac_accapi_dyn_acc_nv.o 
	$(NV_CXX) $(NV_LDFLAGS) -o $@ unifrac_task_api_acc_nv.o unifrac_accapi_dyn_acc_nv.o 

##

unifrac_accapi_acc_amd.o: unifrac_accapi_acc_amd.cpp unifrac_accapi.hpp unifrac_accapi_dyn_acc_amd.h
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@
unifrac_accapi_dyn_acc_amd.o: unifrac_accapi_dyn_acc_amd.cpp unifrac_accapi_dyn_acc_amd.h unifrac_task_impl.hpp
	$(AMD_CXX) $(AMD_CPPFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

unifrac_task_acc_amd.o: unifrac_task_noclass_acc_amd.cpp unifrac_task_noclass.hpp unifrac_task_api_acc_amd.h
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@
unifrac_task_api_acc_amd.o: unifrac_task_api_acc_amd.cpp unifrac_task_api_acc_amd.h unifrac_task_impl.hpp
	$(AMD_CXX) $(AMD_CPPFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

libssu_acc_amd.so: unifrac_task_api_acc_amd.o unifrac_accapi_dyn_acc_amd.o
	$(AMD_CXX) $(AMD_LDFLAGS) -o $@ unifrac_task_api_acc_amd.o unifrac_accapi_dyn_acc_amd.o

##

unifrac_cmp_cpu.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

unifrac_cmp_acc_nv.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_nv -c $< -o $@
unifrac_cmp_acc_amd.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp
	h5c++ $(CPPFLAGS) -DSUCMP_NM=su_acc_amd -c $< -o $@

%.o: %.cpp %.hpp
	h5c++ $(CPPFLAGS) -c $< -o $@

test: test_binaries
	./test_su
	./test_ska
	./test_api
	./test_su_api

install_lib: lib$(SSU).so
	for f in lib$(SSU).so ; do rm -f "${PREFIX}/lib/$${f}"; cp "$${f}" ${PREFIX}/lib/; done


install_lib_acc: $(UNIFRAC_ACC_SHLIBS)
	for f in $(UNIFRAC_ACC_SHLIBS); do rm -f "${PREFIX}/lib/$${f}"; cp "$${f}" ${PREFIX}/lib/; done


install: $(SSU) $(FPD)
	rm -f ${PREFIX}/bin/$(SSU); cp $(SSU) ${PREFIX}/bin/
	rm -f ${PREFIX}/bin/$(FPD); cp $(FPD) ${PREFIX}/bin/
	mkdir -p ${PREFIX}/include/unifrac
	rm -f ${PREFIX}/include/unifrac/task_parameters.hpp; cp task_parameters.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/api.hpp; cp api.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/status_enum.hpp; cp status_enum.hpp ${PREFIX}/include/unifrac/

rapi_test: unifrac_task.cpp unifrac_accapi.cpp tree.o biom.o biom_inmem.o biom_subsampled.o tsv.o unifrac.o cmd.o skbio_alt.o api.o $(UNIFRAC_FILES)
	mkdir -p ~/.R
	if [ -e ~/.R/Makevars ] ; \
	then \
		echo "WARNING: OVERWRITING ~/.R/Makevars" ; \
		echo "The original Makevars file has been copied to ~/.R/Makevars" ;\
		cp ~/.R/Makevars Makevars-original ; \
	fi;
	echo CXX1X=h5c++ > ~/.R/Makevars
	echo CXX=h5c++ >> ~/.R/Makevars 
	echo CC=h5c++ >> ~/.R/Makevars
	echo CFLAGS=  >> ~/.R/Makevars
	echo CXXFLAGS=  >> ~/.R/Makevars
	echo LDFLAGS= -llz4 $(SKBBLIB) >> ~/.R/Makevars
	rm -f *.o
	Rscript --verbose R_interface/rapi_test.R
	rm -f *.o

clean:
	-rm -f *.o $(SSU) $(FPD) *.a test_su test_ska test_api test_su_api lib$(SSU).so unifrac_task_noclass_*.cpp unifrac_task_api_*.h unifrac_task.cpp $(UNIFRAC_ACC_SHLIBS)

clean_install:
	-(cd ${PREFIX}/lib && rm -f lib$(SSU).so $(UNIFRAC_ACC_SHLIBS))
	-(cd ${PREFIX}/bin && rm -f $(SSU) $(FPD))
	-(cd ${PREFIX}/include/unifrac && rm -f task_parameters.hpp api.hpp status_enum.hpp)

