# UniFrac Binaries
##### Canonically pronounced *yew-nih-frak*

UniFrac is the *de facto* repository for high-performance phylogenetic diversity calculations. The methods in this repository are based on an implementation of the [Strided State UniFrac](https://www.nature.com/articles/s41592-018-0187-8) algorithm which is faster, and uses less memory than [Fast UniFrac](http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html). Strided State UniFrac supports [Unweighted UniFrac](http://aem.asm.org/content/71/12/8228.abstract), [Weighted UniFrac](http://aem.asm.org/content/73/5/1576), [Generalized UniFrac](https://academic.oup.com/bioinformatics/article/28/16/2106/324465/Associating-microbiome-composition-with), [Variance Adjusted UniFrac](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-118) and [meta UniFrac](http://www.pnas.org/content/105/39/15076.short), in both double and single precision (fp32).
This repository also includes Stacked Faith (manuscript in preparation), a method for calculating Faith's PD that is faster and uses less memory than the Fast UniFrac-based [reference implementation](http://scikit-bio.org/).

This repository produces standalone executables and a C API exposed via a shared library which can be linked against by any programming language. 

# Citation

A detailed description of the Strided State UniFrac algorithm can be found in [McDonald et al. 2018 Nature Methods](https://www.nature.com/articles/s41592-018-0187-8). Please note that this package implements multiple UniFrac variants, which may have their own citation. Details can be found in the help output from the command line interface in the citations section, and is included immediately below:

    ssu
    For UniFrac, please see:
        Sfiligoi et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
        McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8
        Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005
        Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06
        Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97
        Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133
    For Generalized UniFrac, please see: 
        Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342
    For Variance Adjusted UniFrac, please see: 
        Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118
    For GPU-accelerated UniFrac, please see:
        Sfiligoi et al. PEARC'20; DOI: 10.1145/3311790.3399614

    faithpd
    For Faith's PD, please see:
        Faith Biological Conservation 1992; DOI: 10.1016/0006-3207(92)91201-3

# Install

At this time, there is one primary way to install the library, through `bioconda`. It is also possible to clone the repository and install the C++ API with `make`.

Compilation has been performed on both LLVM 10.0.0 (OS X >= 10.12) or GCC 9 (Centos >= 7) and HDF5 >= 1.8.17. 

Installation time should be a few minutes at most.

## Install (example)

An example of installing UniFrac, and using it with CPUs as well as GPUs, can be be found on [Google Colabs](https://colab.research.google.com/drive/1yL0MdF1zNAkPg1_yESI1iABUH4ZHNGwj?usp=sharing).

## Install (bioconda)

The binaries can be installed through [conda](https://docs.anaconda.com/miniconda/)
via a combination of `conda-forge` and `bioconda` repositories:

```
conda create --name unifrac -c conda-forge -c bioconda unifrac-binaries
conda activate unifrac
```

## Install (native)

To install, first the binary needs to be compiled. This assumes that the HDF5 toolchain and libraries are available.


**Note**: if you are using [conda](https://docs.anaconda.com/miniconda/) we recommend installing HDF5 and related compiler using the
`conda-forge` channel, for example:

On Linux:
```
conda create --strict-channel-priority -n unifrac-binaries -c conda-forge -c bioconda gxx_linux-64 gfortran_linux-64 hdf5 lz4 zlib hdf5-static libcblas liblapacke blas-devel make curl
conda activate unifrac-binaries
```

On MacOs:
```
conda create -q --yes --strict-channel-priority -n unifrac -c conda-forge -c bioconda clangxx_osx-arm64 hdf5 lz4 hdf5-static libcblas liblapacke blas-devel make curl
conda activate unifrac-binaries

```

For NVIDIA-GPU-enabled code, you will need the [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) compiler, and is only supported on Linux.
The NVIDIA GPU compilation requires the setting of the `NV_CXX` environment variable.

This helper script will download it, install it and setup the necessary environment:
```
scripts/install_hpc_sdk.sh 
source setup_nv_compiler.sh
```

For AMD-GPU-enabled code, you will need the [AOMP clang](https://github.com/ROCm/aomp) compiler, and is only supported on Linux.
The AMD GPU compilation requires the setting of the `AMD_CXX` environment variable.

This helper script will download it, install it and setup the necessary environment:
```
scripts/install_amd_clang.sh 
source setup_amd_compiler.sh
```

At this point, we recommend building with
```
export PERFORMING_CONDA_BUILD=True
make clean && make clean_install && make all
```
(This will also install it in the conda location)


If you prefer avoiding conda, assuming `h5c++` is in your path, the following should work, too:

```
make clean && make api && make main 
#optionally
make install
```

Note: If you prefer to build your HDF5 toolchain yourself,
more information about how to setup the 
stack can be found [here](https://support.hdfgroup.org/HDF5/Tutor/compile.html).

# Environment considerations

## Multi-core support

Unifrac uses OpenMP to make use of multiple CPU cores.
By default, Unifrac will use all the cores that are available on the system.
To restrict the number of cores used, set:

    export OMP_NUM_THREADS=nthreads

## Older CPU support

On Linux platforms, Unifrac will auto-detect the CPU generation, i.e. if it supports avx or avx2 vector instructions.
To force the most compatible binary variant, one can set:

    export UNIFRAC_MAX_CPU=basic

To check which binary is used (Unifrac will print it to standard output at runtime), set:

    export UNIFRAC_CPU_INFO=Y

## GPU support

On Linux platforms, Unifrac will run on a GPU, if one is found. 
To disable GPU offload, and thus force CPU-only execution, one can set:

    export UNIFRAC_USE_GPU=N

To check which code path is used (Unifrac will print it to standard output at runtime), set:

    export UNIFRAC_GPU_INFO=Y

Finally, Unifrac will only use one GPU at a time. 
If more than one GPU is present, one can select the one to use by setting:

    export ACC_DEVICE_NUM=gpunum

Note that there is no GPU support for MacOS.

## Additional timing information

When evaluating the performance of Unifrac it is sometimes necessary to distinguish
the time spent interacting with the data from the Unifrac compute proper.
Additional informational messages can be enabled by setting:

    export UNIFRAC_TIMING_INFO=Y

# Examples of use

Below are a few light examples of different ways to use this library.

## Command line

The methods can be used directly through the command line after install:

    $ which ssu
    /Users/<username>/miniconda3/envs/unifrac/bin/ssu
    $ ssu --help
    usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-a alpha] [-f]  [--vaw]
        [--mode MODE] [--start starting-stripe] [--stop stopping-stripe] [--partial-pattern <glob>]
        [--n-partials number_of_partitions] [--report-bare] [--format|-r out-mode]
        [--n-substeps n] [--pcoa dims] [--diskbuf path]

        -i		The input BIOM table.
        -t		The input phylogeny in newick.
        -m		The method, [unweighted | weighted_normalized | weighted_unnormalized | unweighted_unnormalized | generalized |
                           unweighted_fp64 | weighted_normalized_fp64 | weighted_unnormalized_fp64 |
                           unweighted_unnormalized_fp64 | generalized_fp64 |
                           unweighted_fp32 | weighted_normalized_fp32 | weighted_unnormalized_fp32 |
                           unweighted_unnormalized_fp32 | generalized_fp32].
        -o		The output distance matrix.
        -g		[OPTIONAL] The input grouping in TSV.
        -c		[OPTIONAL] The columns(s) to use for grouping, multiple values comma separated.
        -a		[OPTIONAL] Generalized UniFrac alpha, default is 1.
        -f		[OPTIONAL] Bypass tips, reduces compute by about 50%.
        --vaw	[OPTIONAL] Variance adjusted, default is to not adjust for variance.
        --mode	[OPTIONAL] Mode of operation:
                                one-off : [DEFAULT] compute UniFrac.
                                partial : Compute UniFrac over a subset of stripes.
                                partial-report : Start and stop suggestions for partial compute.
                                merge-partial : Merge partial UniFrac results.
                                multi : compute UniFrac multiple times.
        --start	[OPTIONAL] If mode==partial, the starting stripe.
        --stop	[OPTIONAL] If mode==partial, the stopping stripe.
        --partial-pattern	[OPTIONAL] If mode==merge-partial, a glob pattern for partial outputs to merge.
        --n-partials 	[OPTIONAL] If mode==partial-report, the number of partitions to compute.
        --report-bare	[OPTIONAL] If mode==partial-report, produce barebones output.
        --n-substeps 	[OPTIONAL] Internally split the problem in n substeps for reduced memory footprint, default is 1.
        --format|-r	[OPTIONAL]  Output format:
                                 ascii : Original ASCII format. (default if mode==one-off)
                                 hdf5_nodist : HFD5 format, no distance matrix. (default if mode==multi)
                                 hdf5 : HFD5 format.  May be fp32 or fp64, depending on method.
                                 hdf5_fp32 : HFD5 format, using fp32 precision.
                                 hdf5_fp64 : HFD5 format, using fp64 precision.
        --subsample-depth   Depth of subsampling of the input BIOM before computing unifrac (required for mode==multi, optional for one-off)
        --subsample-replacement	[OPTIONAL] Subsample with or without replacement (default is with)
        --n-subsamples	[OPTIONAL] if mode==multi, number of subsampled UniFracs to compute (default: 100)
        --permanova	[OPTIONAL] Number of PERMANOVA permutations to compute (default: 999 with -g, do not compute if 0)
        --pcoa	[OPTIONAL] Number of PCoA dimensions to compute (default: 10, do not compute if 0)
        --seed	[OPTIONAL] Seed to use for initializing the random gnerator
        --diskbuf	[OPTIONAL] Use a disk buffer to reduce memory footprint. Provide path to a fast partition (ideally NVMe).
        -n		[OPTIONAL] DEPRECATED, no-op.

    Environment variables: 
        CPU parallelism is controlled by OMP_NUM_THREADS. If not defined, all detected core will be used.
        GPU offload can be disabled with UNIFRAC_USE_GPU=N. By default, if a NVIDIA GPU is detected, it will be used.
        A specific GPU can be selected with ACC_DEVICE_NUM. If not defined, the first GPU will be used.

    Citations: 
        For UniFrac, please see:
            Sfiligoi et al. mSystems 2022; DOI: 10.1128/msystems.00028-22
            McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8
            Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005
            Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06
            Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97
            Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133
        For Generalized UniFrac, please see: 
            Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342
        For Variance Adjusted UniFrac, please see: 
            Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118

    $ which faithpd
    /Users/<username>/miniconda3/envs/unifrac/bin/faithpd
    $ faithpd --help
	usage: faithpd -i <biom> -t <newick> -o <out.txt>

		-i          The input BIOM table.
		-t          The input phylogeny in newick.
		-o          The output series.

	Citations: 
		For Faith's PD, please see:
			Faith Biological Conservation 1992; DOI: 10.1016/0006-3207(92)91201-3

            
## Shared library access

In addition to the above methods to access UniFrac, it is also possible to link against the shared library. The C API is described in `src/api.hpp`, and examples of linking against this API can be found in `examples/`. 

## Minor test dataset

A small test `.biom` and `.tre` can be found in `src/`. An example with expected output is below, and should execute in 10s of milliseconds:

    $ ssu -i src/test.biom -t src/test.tre -m unweighted -o test.out
    $ cat test.out
    	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
    Sample1	0	0.2	0.5714285714285714	0.6	0.5	0.2
    Sample2	0.2	0	0.4285714285714285	0.6666666666666666	0.6	0.3333333333333333
    Sample3	0.5714285714285714	0.4285714285714285	0	0.7142857142857143	0.8571428571428571	0.4285714285714285
    Sample4	0.6	0.6666666666666666	0.7142857142857143	0	0.3333333333333333	0.4
    Sample5	0.5	0.6	0.8571428571428571	0.3333333333333333	0	0.6
    Sample6	0.2	0.3333333333333333	0.4285714285714285	0.4	0.6	0
