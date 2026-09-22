# molConv
Monte Carlo simulation of molecule creation from ultracold atoms

## Overview

This C++ code performs a static (i.e. no temporal change) Monte Carlo simulation where one or two species of ultracold atoms with Ferion or Boson statistics (including non-interacting BEC) or thermal gases are created in a harmonic trap and from which molecules are created using different molecule pairing criteria. Varying species, atom number, temperature, T/TF or T/Tc, or trap geometry allows to calculate molecule creation efficiency for different pairing criteria and to compare with experiment data. The Monte Carlo simulation with the different pairing criteria is inspired by [Tomotake Yamakoshi, Shinichi Watanabe, Chen Zhang, and Chris H. Greene, `Stochastic and equilibrium pictures of the ultracold Fano-Feshbach-resonance molecular conversion rate`, Physical Review A87, 053604 (2013)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.87.053604) [Ref. 1].

The present code originates from the version by Florian Schreck (now Amsterdam) used by the LiK experiment in Innsbruck/Austria (around 2010). The improvements include multi-threading on Windows and Linux, the implementation of the Metropolis argorithm for efficient atom generation, and the choice of different 64bit or 128bit random number generators which were extensively tested with PractRand (v0.94). The molecule search algorithm has been extended with more pairing criteria and the nearest-pair option. For future applications (like temporal evolution) a cells grid has been introduced but which is not yet efficiently used. Although this code is tailored for the specific purpose of molecule formation it might serve as a starting point for extension and other studies.

The Monte Carlo simulation shows that for the mass-balanced case only T/TF or T/Tc matters and the difference between the pairing criteria behave is negligible apart of a scaling factor for each cirterium chosen to match the experiment. Below a figure for the mass-balanced Fermi-Fermi case (similar to Fig.2a of Ref.1, MB = thermal gas):

<img src="figures/Grene_Fig2a.png" width="600"/>

For the mass-imbalanced case at sufficient low temperature the different criteria deviate sufficiently from each other such that by comparison with experiment one might be able to pinpoint the microscopic origin of the pairing mechanism. Below a figure for the mass-imbalanced Fermi-Bose case (similar to Fig.2c of Ref. 1, MB = thermal gas):

<img src="figures/Grene_Fig2c.png" width="600"/>

In our experiment we have Li6 and different isotopes of Chromium (50, 52 and 53) and we can create molecules with all combinations of species with a mass imbalance of about 8 and different combinations of statistics. We measure the pairing efficiency at different T/TF or T/TF of these mixtures and the comparison with the results of the Monte Carlo simulation might lead to a better understanding of the microsopic picture.

The project structure is below:

```
├── figures     example figures
├── log         PractRand 0.94 logfiles
├── params      example params files
├── python      python scripts
└── src         C++ source code
```

## Compiling

The [/src/ source folder](/src/) contains the C++ files and a Makefile for easy compilation on Linux systems (tested with GCC on Ubuntu 22.04). In the terminal just type `make` in the folder containing the source files. The compilation requires pthreads (part of the newer GCC versions) and has no other dependencies.

For compilation on Windows a Visual Studio (2019) solution can be found in the the [/src/Visual-Studio/ folder](/src/VS_solution). Copy this into the same folder as the source files, open the solution with Visual Studio (follow instructions in case it wants to upgrade the project) and compile it (tested on Windows 10). The project is however also easily created manually: create an empty new `win32` project, add all existing files to the project and select `x64 Release` and compile. The project has no additional dependencies and should be easily compiled. The Lehmer64/128 random number generators use 128bit integer arithmetics which is not available on the Visual Studio compiler. Therefore a custom implementation (using GCC source code) is provided. It should be however possible to compile the project on Visual Studio with `clang` compiler which is natively supporting 128bit integers. I have not yet tested this however.

A list of compile options can be found at the end of this page.

## How to use

The most basic usage is to execute `./molConv` or `molconv.exe` from the command line where a few command line arguments can be given [see list of command line arguments below][]. 

MolConv accepts many more options which can be specified in the parameter file using the `-p` argument. The [/params/](/params/) folder which has a few example parameter files.

If no seed value is given to MolConv then it uses `std::random_device` to generate two (hardware) seed values which is used to generate a seed sequence for the used random number generator which is different for each thread. These two hardware seed values are printed on the screen and are saved into the result file such that if needed the user can provide them to a subsquent call to MolConv and to repeat the simulation. MolConv ensures that even though the seed values are also used for internal self-tests, the atom creation is always initialize with exactly the provided seed values, such that even if the tests would be changed in the future or different number of tests executed the result is still the same.

The most important output file is the result file (specified by `result_file` parameter in the parameters file). This is a text file which summarizes all input values, calculated T/TF and T/Tc and other critical parameters, provides the molecule conversion efficiency and more information and estimated quantities of the molecule. Its structure is optimized for automatic processing and is not very convenient to read. It also contains all information for each repetition and summaries of atoms and molecules and repetitions. 

For easiest use of MolConv the `molConv_run.py` python script is provided in the [/python/](/python/) folder which not only generates the parameter file, but also launches molConv and collects the results saved into the result file. It generates a result csv file with the averaged result from all repetitions and generates `nice` plots.


## List of command line arguments

<!--
| argument  | description                                                                 | example         | use case
|-----------|-----------------------------------------------------------------------------|-----------------|--------------------|
| -p <file> | read options from parameter file <file name>                                | -s myParams.txt | simulation         |
| -g <RNG>  | select random number generator <RNG name>                                   | -g Lehmer128    | simulation         |
| -s <seed> | give min. 1 comma-separated seed value within {}, no spaces, hex with 0x    | -s {1,2,0xab}   | reproduce old sim. |
| -tI       | stream uniform distributed 64bit integer to stdout                          | -tI             | PractRand          |     
| -tD       | stream uniform distributed double converted to 64bit integer to stdout      | -tD             | PractRand          |
| -tN       | stream normal  distributed double converted to 64bit integer to stdout      | -tN             | PractRand          |
| -ti       | stream uniform distributed 32bit integer to stdout                          | -ti             | PractRand          |
| -td       | stream uniform distributed float converted to 32bit integer to stdout       | -td             | PractRand          |
| -tn       | stream normal  distributed float converted to 32bit integer to stdout       | -tn             | PractRand          |
-->

## List of parameter file entrie

<!--
| parameter               | description                                                | remark                                  |
|-------------------------|------------------------------------------------------------|-----------------------------------------|
| species_0               | name of first atom                                         |                                         |
| species_1               | name of second atom                                        | *                                       |
| num_species             | number of species. must be 1 or 2                          | must be 1 or 2                          |
| statistics_0            | statistics of first atom                                   | see list of statistics                  |
| statistics_1            | statistics of second atom                                  | see list of statistics                  |
|                         |                                                            | *                                       |
| statistics_mol          | statistics of molecules                                    | see list of statitsics                  |
| mass_0                  | mass of first atom in amu                                  |                                         |
| mass_1                  | mass of second atom in amu                                 | *                                       |
| N_0                     | atom number of first atom                                  |                                         |
| N_1                     | atom number of second atom                                 | *                                       |
| T_0                     | temperature of first atom in nK                            |                                         |
| T_1                     | temperature of second atom in nK                           | *                                       |
| fx_0                    | harmonic trapping frequency of first atom on x-axis in Hz  |                                         |
| fy_0                    | harmonic trapping frequency of first atom on y-axis in Hz  |                                         |
| fz_0                    | harmonic trapping frequency of first atom on z-axis in Hz  |                                         |
| fx_1                    | harmonic trapping frequency of second atom on x-axis in Hz | *                                       |
| fy_1                    | harmonic trapping frequency of second atom on y-axis in Hz | *                                       |
| fz_1                    | harmonic trapping frequency of second atom on z-axis in Hz | *                                       |
| cx                      | offset of second atom center x-position in mu              | default = 0, *                          |
| cy                      | offset of second atom center y-position in mu              | default = 0, *                          |
| cz                      | offset of second atom center z-position in mu              | default = 0, *                          |
| size_0                  | calculation size and energy scaling                        | default = 1.0                           |
| size_1                  | calculation size and energy scaling                        | default = 1.0, *                        |
| error_allowed           | chemical potential calc. goal |delta N/N|                  | default = 1e-10                         |
| step_size               | chemical potential calc. initial step size                 | default = 0.1                           |
| repeat                  | repetition of simulation                                   | default = 1                             |
| distance_measure        | molecule pairing criterion                                 | see list of pairing criteria            |
| gamma                   | maximum distance to take for given pairing criterium       | obtained by comparison with experiment  |
| threads                 | number of threads                                          | default = 8                             |
| num_bins                | number of bins for histogram                               | default = 300                           |
| random_number_generator | random number generator                                    | see list of random number generators    |
|                         |                                                            | default = 'Lehmer128'                   |
| distribution_generator  | distribution generator                                     | see list of distribution generators     |
|                         |                                                            | default = Metropolis                    |
| seed                    | seed values given as comma-separated list within {}        | default = "" = hardware seed            |
|                         |                                                            | min. 1 value, no spaces allowed         |
| result_file             | path and filename of result text file                      | default = result.dat                    |
| histogram_file          | path and filename of histogram text file                   | default = "" = no file created, **      |
| export_atom_0           | path and filename of exported first atoms csv file         | default = "" = no file created, **      |
| export_atom_1           | path and filename of exported second atoms csv file        | *, **                                   |
|                         |                                                            | default = "" = no file created, **      |
| export_molecule         | path and filename of exported molecules csv file           | default = "" = no file created, **      |
| neighbors_file          | path and filename of exported neightbors csv file          | at the moment disabled                  |
| neighbors_num           | number of neighbors to export                              | at the moment disabled                  |
| find_nearest            | if nonzero molecule search pairs nearest neighbors (slow)  | default = 0                             |
|                         | otherwise first matching pair is used (faster)             |                                         |
-->

`*` can only provided for num_species == 2
`**` histogram, atoms and molecule file names are appended with species/molecule name and repetition

# List of statistics

These statistics are implemented. For `BoseEinstein` two seperate species are generated (`thermal` and `condensed`) depending on the thermal and BEC fraction. For the BEC part non-interacting (harmonic ground-state wavefunction) is assumed using the truncated Wigner approximation Equ. 4 in Ref. 1.

<!--
| statistics               | description                                                                                         |
|--------------------------|-----------------------------------------------------------------------------------------------------|
| FermiDirac               | ideal Fermi gas, see Ketterle, Durfee, Stamper Kurn,                                                |
|                          | Making, probing and understanding Bose-Einstein condensates,                                        |
|                          | https://arxiv.org/abs/cond-mat/9904034                                                              |
| BoseEinstein             | non-interacting BEC and thermal Bose gas, see Wolfgang Ketterle, Martin W. Zwierlein,               |
|                          | Making, probing and understanding ultracold Fermi gases,                                            |
|                          | https://arxiv.org/abs/0801.2500                                                                     |
| MaxwellBoltzman          | thermal gas, see Refs. for other statistics                                                         |
-->

# List of pairing criteria

These are the possibile pairing criteria which include the ones from Ref. 1.

<!--
| criterium                | description                                                                      | Ref. 1 criterion |
|--------------------------|----------------------------------------------------------------------------------|------------------|
| delta_x                  | real space distance                                                              | -                |
| delta_v                  | velocity difference                                                              | -                |
| delta_x*delta_v          | momentum-space distance in center of mass frame                                  | -                |
| delta_x*delta_p          | momentum-space distance                                                          | 1                |    
| cross(delta_x^,delta_p^) | partial waves                                                                    | 2                |
| <delta_x*delta_p>        | phase-space volume                                                               | 3                |
| max\|delta_xi*delta_pi\| | all individual phase-space coordinates must match                                | 4                |
-->

# List of recommended random number generators

These random number generators (RNGs) are enabled by default and are recommended for use. 

All of these RNGs have been tested extensively with the PractRand (v0.94) test suite (using the molConv_PractRand.py script together with the command line arguments `-t#`. All, except Wyhash64, fail the `-tI` test but pass the `-tD` test (for > 1TB). In the first test the generated uniform-distributed random (64bit) numbers are sent directly to PractRand while in the second test the derived random double-precision floats (in the range 0..1) are converted back into 64bit integers before sending to PractRand. The conversion from double gives only 53bits which requires that two random doubles are combined. This procedure obviously masks the correlations between the bits generated by the original RNG. The Lehmer generators pass however PractRand 0.93 test (see Ref. 6 and comments in Ref. 3) but they intrinsically fail more modern statistical tests. Since the `-tD` test passes which is more close to the use case of the Monte Carlo I beleive all of the provided RNGs are suitable. If you want to be sure you can make a comparison simulation with `Wyhash64` which passes both the `-tI` and `-tD` tests (with more than 1TB size). This is more modern but not as established as the other generators. 

All of tese RNGs fail the PractRand test with the `-tN` option, i.e. where the normal-distributed random doubles are converted back to a uniform distribution. I do not consider this as a serious flaw since most likely the problem is the conversion back into uniform distribution and not of the normal-distributed random numbers. These are used for the Metropolis algorithm. If you want to be sure you can use the older rejection-sampling method for atom generaton - which is much slower but requires only the uniform distributed doubles tested with the `-tD` option.

MersenneTwiswer64 is used to generate a variable size seed sequence for the other random number generators unique for each thread. This overcomes the limitations in the implemenation of `std::seed_seq` which (on Ubuntu 22.04) creates for each call to `generate` always the same sequence of seed values and cannot be reset to a new state.

<!--
| RNG               | description                                                    |speed | Ref. | PractRand -tI test fails at | 
|-------------------|----------------------------------------------------------------|------|------|-----------------------------|
| MersenneTwister64 | from c++11 standard library, 2^19937-1 period                  | slow |      | 512GB                       |
| Lehmer64          | linear congruential generator, state 128bit, unclear period    | fast | 2,3  |  64GB                       |
| Lehmer128         | linear congruential generator, state 128bit, 2^126 period      | fast | 3,4  |  64GB                       |
| Wyhash64          | more modern, state 64bit, unclear period, passes PractRand!    | fast |   5  |   -                         |
-->


Ref. 2: 
D. H. Lehmer, Mathematical methods in large-scale computing units. Proceedings of a Second Symposium on Large Scale Digital Calculating Machinery; Annals of the Computation Laboratory, Harvard Univ. 26 (1951), pp. 141-146.

Ref 3:
https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/
https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h

Ref. 4:
https://en.wikipedia.org/wiki/Lehmer_random_number_generator
P L'Ecuyer,  Tables of linear congruential generators of different sizes and good lattice structure. 
Mathematics of Computation of the American Mathematical Society 68.225 (1999): 249-260.

Ref. 5:
https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
adapted to this project by D. Lemire, from https://github.com/wangyi-fudan/wyhash/blob/master/wyhash.h

Ref. 6:
https://www.pcg-random.org/posts/does-it-beat-the-minimal-standard.html
https://www.pcg-random.org/posts/too-big-to-fail.html

# List of not recommended random number generators

> [!WARNING]
> Do not use these random number generators for the Monte Carlo simulation!

These are for reference and testing purposes only! They have a very limited internal state and short period (`WELL1024` might be an exception and `Ranlux48` is intermediate). These are normally disabled and must be enabled by compiling with the (not recommended) option `REAL_PRECISION` = `REAL_SINGLE` which reduces the floating point precision to single. Using `Lehmer32` with the very inefficient `rejection-sampling` method and the less efficient method to generate random floating-point values I have observed clear `noise` in the y-axis histograms already for 2x15k generated atoms! Changing to `Lehmer64` this noise disappears. The 32bit RNG simply runs out of random numbers and starts to repeat the sequence under these conditions.

<!--
| RNG               | description                                                                                                | 
|-------------------|------------------------------------------------------------------------------------------------------------|
| MinStd            | minimum standard, from c++11 standard library, 32bit state                                                 |
| Ranlux24          | from c++11 standard library, 24bit state, very slow                                                        |
| Ranlux48          | from c++11 standard library, 48bit state, very slow                                                        |
| Lehmer32          | original random number generator used by Florian, fast, similar or even the same as MinStd                 |
| WELL1024          | copied from one of my old codes, supposed to be good at that time, large internal state, hard to seed      |
|                   | passed a small size `-tD` test but needs more testing.                                                     |
-->

# List of distribution generators

There are only two possible distribution generators available. `Metropolis` is the default since it is much faster, but it creates a less precise distribution than the classical rejection-sampling algorithm which is extremely slow even on multiple threads. I keep both here in case one has to check the effect of the more noisy Metropolis algorithm.

<!--
| generator                | description                                                                                         |
|--------------------------|-----------------------------------------------------------------------------------------------------|
| rejection-sampling       | original sampling method, extremely slow, produces a very precise distribution                      |
| Metropolis               | very fast sampling method, produces more variation of the samples                                   |
-->

## List of compile options

These are the most important compile time options found in the `MoleculeConversion.h` header file.

<!--
| option                   | description                                                                                         |
|--------------------------|-----------------------------------------------------------------------------------------------------|
| _DEBUG                   | if defined additional debugging code is executed and more output is generated                       |
|                          | slows down execution                                                                                |
| CLASS_UINT128            | if defined uses uint128.h to implement uint128_t arithmetics.                                       |
| REAL_PRECISION           | if this is set to REAL_DOUBLE double precision is used, otherwise single precision                  |
|                          | double precision is recommended, single precision was used to test if faster calculation is         |
|                          | possible (had no big effect, but should be tested again)                                            |
|                          | and to enable the not recommended 32bit random number generators                                    |
| SEED_VALUES_NUM          | number of hardware seed values generated. default = 2 is a good compromize of large enough          |
|                          | and not loosing too much entropy by excessive use                                                   |
| STREAM_DATA_OUT          | if defined gives the number of samples which are output in parallel on stdout and stderr            |
|                          | by RNG_stream for tests '-tI' and '-ti'. this allows to debug if on Windows the pipe                |
|                          | is in text or binary mode which is needed for PractRand to work properly.                           |
-->



