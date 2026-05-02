.PHONY: test clean all clean_install wasm wasm_test wasm_clean

# Serialize top-level targets ONLY when a wasm target is involved. The
# native generators for unifrac_accapi_cpu.cpp / unifrac_task_noclass_cpu.cpp
# are shared between the native and WASM sub-makes; a concurrent
# `make all wasm` could fire the python generators twice and truncate
# each other's output. Native-only `make -j8 all` is unaffected.
ifneq (,$(filter wasm wasm_test wasm_clean,$(MAKECMDGOALS)))
.NOTPARALLEL:
endif

PLATFORM := $(shell uname -s)

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
ARCH := $(shell h5c++ -dumpmachine)

ifneq (,$(findstring x86_64,$(ARCH)))
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

else
# only one optimization level for non-x86 architectures

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
	$(MAKE) all_combined

clean_cpu:
	-cd test && $(MAKE) clean
	-export BUILD_VARIANT=cpu_basic; cd src && $(MAKE) clean
	-cd combined && $(MAKE) clean

clean_install_cpu:
	-export BUILD_VARIANT=cpu_basic; cd src && $(MAKE) clean_install
	-cd combined && $(MAKE) clean_install

all_cpu_basic:
	$(MAKE) api_cpu_basic
	$(MAKE) install_cpu_basic


all_acc:
	$(MAKE) api_acc
	$(MAKE) install_lib_acc

all_combined:
	$(MAKE) api_combined
	$(MAKE) install_combined
	$(MAKE) main
	$(MAKE) install_main

endif

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
	cd src && $(MAKE) clean && $(MAKE) api_acc

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
	cd src && $(MAKE) install_lib_acc

install_combined:
	cd combined && $(MAKE) install

########### test

test_binaries:
	cd src && $(MAKE) clean && $(MAKE) test_binaries
	cd test && $(MAKE) clean && $(MAKE) test_binaries

test:
	cd src && $(MAKE) test
	cd test && $(MAKE) test

########### WASM

# Build libssu_wasm.a. Requires that scikit-bio-binaries' libskbb_wasm.a
# already exists (run `make wasm` in that repo first), and that an emsdk
# is activated on PATH. Override SKBB_DIR for a non-sibling skbb checkout:
#   make wasm SKBB_DIR=/path/to/scikit-bio-binaries
#
# Default sibling layout (../scikit-bio-binaries) is checked first; if
# its libskbb_wasm.a is missing, this target builds it before proceeding.
SKBB_DIR ?= $(abspath $(CURDIR)/../scikit-bio-binaries)

wasm: $(SKBB_DIR)/src/libskbb_wasm.a
	cd src && $(MAKE) SKBB_DIR=$(SKBB_DIR) wasm

wasm_test: $(SKBB_DIR)/src/libskbb_wasm.a
	cd src && $(MAKE) SKBB_DIR=$(SKBB_DIR) wasm_test

wasm_clean:
	cd src && $(MAKE) wasm_clean

# Sentinel rule: ensure skbb's WASM artifact exists. If not, build it.
# FORCE-prereq trick keeps Make from re-entering skbb on every invocation
# once the file exists; once it does, this rule is a no-op test. If skbb
# itself changes, downstream callers are expected to rebuild it explicitly,
# which is the same contract skbb's other consumers have.
$(SKBB_DIR)/src/libskbb_wasm.a: FORCE
	@test -f $@ || (cd $(SKBB_DIR) && $(MAKE) wasm)

FORCE:
.PHONY: FORCE

