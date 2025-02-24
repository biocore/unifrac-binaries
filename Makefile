.PHONY: test clean all clean_install

# Note: This Makefile will NOT properly work with the -j option 

PLATFORM := $(shell uname -s)
COMPILER := $(shell ($(CXX) -v 2>&1) | tr A-Z a-z )

ifeq ($(PLATFORM),Darwin)
# only one optimization level and no GPU support for MacOS
all:
	$(MAKE) api
	$(MAKE) install
	$(MAKE) main
	$(MAKE) install_main
	$(MAKE) test_binaries

clean:
	-cd test && $(MAKE) clean
	-cd src && $(MAKE) clean

clean_install:
	-cd src && $(MAKE) clean_install

else
# Linux with several optimization levels and with optional GPU support

all: 
	$(MAKE) all_cpu
	$(MAKE) all_acc
	$(MAKE) test_binaries

clean:
	$(MAKE) clean_cpu

clean_install:
	$(MAKE) clean_install_cpu

all_cpu: 
	$(MAKE) all_cpu_basic
	$(MAKE) all_cpu_x86_v2
	$(MAKE) all_cpu_x86_v3
	$(MAKE) all_cpu_x86_v4
	$(MAKE) all_combined

clean_cpu:
	-cd test && $(MAKE) clean
	-export BUILD_VARIANT=cpu_basic; cd src && $(MAKE) clean
	-export BUILD_VARIANT=cpu_x86_v2; cd src && $(MAKE) clean
	-export BUILD_VARIANT=cpu_x86_v3; cd src && $(MAKE) clean
	-export BUILD_VARIANT=cpu_x86_v4; cd src && $(MAKE) clean
	-cd combined && $(MAKE) clean

clean_install_cpu:
	-export BUILD_VARIANT=cpu_basic; cd src && $(MAKE) clean_install
	-export BUILD_VARIANT=cpu_x86_v2; cd src && $(MAKE) clean_install
	-export BUILD_VARIANT=cpu_x86_v3; cd src && $(MAKE) clean_install
	-export BUILD_VARIANT=cpu_x86_v4; cd src && $(MAKE) clean_install
	-export BUILD_VARIANT=nv; cd src && $(MAKE) clean_install
	-cd combined && $(MAKE) clean_install

all_cpu_basic:
	$(MAKE) api_cpu_basic
	$(MAKE) install_cpu_basic

all_cpu_x86_v2:
	$(MAKE) api_cpu_x86_v2
	$(MAKE) install_cpu_x86_v2

all_cpu_x86_v3:
	$(MAKE) api_cpu_x86_v3
	$(MAKE) install_cpu_x86_v3

all_cpu_x86_v4:
	$(MAKE) api_cpu_x86_v4
	$(MAKE) install_cpu_x86_v4

all_acc: 
	$(MAKE) api_acc
	$(MAKE) install_lib_acc

all_combined:
	$(MAKE) api_combined
	$(MAKE) install_combined
	$(MAKE) main
	$(MAKE) install_main

endif

########### api

api:
	cd src && $(MAKE) clean && $(MAKE) api

api_cpu_basic:
	export BUILD_VARIANT=cpu_basic ; export BUILD_FULL_OPTIMIZATION=False ; cd src && $(MAKE) clean && $(MAKE) api

api_cpu_x86_v2:
	export BUILD_VARIANT=cpu_x86_v2 ; export BUILD_FULL_OPTIMIZATION=x86-64-v2 ; export BUILD_TUNE_OPTIMIZATION=core2; cd src && $(MAKE) clean && $(MAKE) api

api_cpu_x86_v3:
	export BUILD_VARIANT=cpu_x86_v3 ; export BUILD_FULL_OPTIMIZATION=x86-64-v3 ; export BUILD_TUNE_OPTIMIZATION=znver3; cd src && $(MAKE) clean && $(MAKE) api

api_cpu_x86_v4:
	export BUILD_VARIANT=cpu_x86_v4 ; export BUILD_FULL_OPTIMIZATION=x86-64-v4 ; export BUILD_TUNE_OPTIMIZATION=znver4 ;cd src && $(MAKE) clean && $(MAKE) api

api_acc:
	export BUILD_FULL_OPTIMIZATION=False ; cd src && $(MAKE) clean && $(MAKE) api_acc

api_combined:
	cd combined && $(MAKE) clean && $(MAKE) api

########### main

main:
	cd src && $(MAKE) clean && $(MAKE) main

install_main:
	cd src && $(MAKE) install

########### install

install:
	cd src && $(MAKE) install_lib

install_cpu_basic:
	export BUILD_VARIANT=cpu_basic ; export BUILD_FULL_OPTIMIZATION=False ; cd src && $(MAKE) install_lib

install_cpu_x86_v2:
	export BUILD_VARIANT=cpu_x86_v2 ; export BUILD_FULL_OPTIMIZATION=x86-64-v2 ; cd src && $(MAKE) install_lib

install_cpu_x86_v3:
	export BUILD_VARIANT=cpu_x86_v3 ; export BUILD_FULL_OPTIMIZATION=x86-64-v3 ; cd src && $(MAKE) install_lib

install_cpu_x86_v4:
	export BUILD_VARIANT=cpu_x86_v4 ; export BUILD_FULL_OPTIMIZATION=x86-64-v4 ; cd src && $(MAKE) install_lib

install_lib_acc:
	export BUILD_FULL_OPTIMIZATION=False ; cd src && $(MAKE) install_lib_acc

install_combined:
	cd combined && $(MAKE) install

########### test

test_binaries:
	cd src && $(MAKE) clean && $(MAKE) test_binaries
	cd test && $(MAKE) clean && $(MAKE) test_binaries

test:
	cd src && $(MAKE) test
	cd test && $(MAKE) test

