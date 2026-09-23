// this code is based in parts on the version used by the LiK experiment in Innsbruck/Austria (around 2010) developed by Florian Schreck (now Amsterdam)
// revised, debugged and extended by Andi
// last change 23/09/2026 by Andi

#if defined(_WIN32) || defined(_WIN64) 
// define early for windows. does not work when inside MolecularConversion.h
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <chrono>
#include <ctime>
#include <algorithm>

#include "MolecularConversion.h"    // cells and species
#include "random.hpp"               // random number generators

std::string seed_use = "";

//Physical constants in SI units from NIST/CODATA 
const REAL_TYPE h                  = 6.62607015E-34;           // checked 3/8/2025, exact
const REAL_TYPE hbar               = 1.054571817E-34;          // checked 3/8/2025, exact
const REAL_TYPE kB                 = 1.380649E-23;             // checked 3/8/2025, exact
const REAL_TYPE Pi                 = 3.1415926535897932385;    // checked 3/8/2025 with Mathematica
const REAL_TYPE mn                 = 1.67492750056e-27;        // checked 3/8/2025, error (85), old value 1.67492749804e-27
//const REAL_TYPE amu                = 1.66053906892e-27;        // checked 3/8/2025, error (52), at the moment not used.
const REAL_TYPE zeta_3             = 1.2020569031595942854;    // checked 3/8/2025 with Mathematica, zeta(3) = Li_3(1)

// conversion constants. updated 3/8/2025 with constants above (mn changed slightly)
// TODO: I use here mn but also amu could be used which is slightly (<1%) different. not clear which is more appropriate? 
const REAL_TYPE k_Epot = (1.67492750056/1.380649)*1e-7;        // get [Epot] = amu*(rad/s*mu)^2 -> nK; mn*(1e-6)^2/kB*1e9 = mn/kB*1e-3;
const REAL_TYPE k_Ekin = (1.67492750056/1.380649)*1e-1;        // get [Ekin] = amu*(mm/s)^2     -> nK; mn*(1e-3)^2/kB*1e9 = mn/kB*1e3;
const REAL_TYPE k_nK   = (1.054571817  /1.380649)*1.0e-2;      // get TF and Tc in nK out of hbar*omega: k_nK = (hbar*1e9)/kB;
const REAL_TYPE k_h    = (6.62607015   /1.67492750056)*1e2;    // scale phase space density (gamma): k_h = h/(mn*1e-9)
const REAL_TYPE k_Rc   = std::sqrt(1.054571817/1.67492750056*1e5);  // get BEC size in mu = 1e6*np.sqrt(hbar/mn)

// distance measures string and corresponding distance index
// note: no spaces allowed here.
#define NUM_DISTANCE            7       // user-selections for "distance_measure"
#define DISTANCE_X              0       // real space distance            :    |delta x|              <= gamma
#define DISTANCE_V              1       // real space distance            :    |delta v|              <= gamma
#define DISTANCE_PV             2       // velocity-space distance        :    |delta x| * |delta v|  <= gamma
#define DISTANCE_PP             3       // momentum-space distance        :    |delta x| * |delta p|  <= gamma
#define DISTANCE_CROSS          4       // cross momentum-space distance  :    |delta x  x  delta p|  <= gamma
#define DISTANCE_MEAN           5       // mean momentum-space distance   :   <|delta x| * |delta p|> <= gamma
#define DISTANCE_MAX            6       // maximum momentum-space distance: max|delta xi * delta pi|  <= gamma
std::string distance_all[NUM_DISTANCE]  = {"delta_x",
                                           "delta_v",
                                           "delta_x*delta_v",
                                           "delta_x*delta_p",
                                           "cross(delta_x^,delta_p^)", 
                                           "<delta_x*delta_p>",
                                           "max|delta_xi*delta_pi|"
                                          };

// distance measure used for molecule search
uint8_t     mol_distance_index          = DISTANCE_PP;
std::string mol_distance_use            = distance_all[mol_distance_index];
REAL_TYPE      mol_gamma = 0.38;           // threshold in units of Planck constant or mu or mm/s depending on distance used
//REAL_TYPE      k_gamma   = 0;              // squared and scaled mol_gamma depending on distance used

// statistics 
#define NUM_STAT                5
#define STAT_MB                 0       // Maxwell-Boltzmann = thermal gas at any T
#define STAT_FD                 1       // Fermi-Dirac       = Fermions at any T
#define STAT_BE                 2       // Bose-Einstein     = Bosons at any T: non-condensed + condensed parts
#define STAT_BE_NC              3       // Bosons: only non-condensed part
#define STAT_BE_C               4       // Bosons: only condensed part - non-interating - Gaussian ground state of HO. 
//#define STAT_NONE               0xff    // invalid statistics
// TODO: implement interacting BEC = inverted parabola
std::string stat_all        [NUM_STAT]  = { "MaxwellBoltzmann",
                                            "FermiDirac",
                                            "BoseEinstein", 
                                            "BoseEinstein(thermal)",
                                            "BoseEinstein(condensed)"};
uint8_t distance_default    [NUM_STAT]  = {DISTANCE_X, DISTANCE_PP, DISTANCE_PP, DISTANCE_PP, DISTANCE_PP};

// global double-linked list of shared cells containing atoms and optional cells with different distance measure.
// each list is sorted with individual distance measure.
class cells_list *cells[NUM_CELLS_LISTS] = {NULL};

// global double-linked list of shared nearest neighbors found by FindNearest
class double_linked_list<class species_entry, uint64_t> nearest;

// distribution generator
#define NUM_DSTRB               2
#define DSTRB_BASIC             0
#define DSTRB_METROPOLIS        1
std::string dstrb_all[NUM_DSTRB]    = {"rejection-sampling","Metropolis"};
uint8_t     dstrb_index             = DSTRB_METROPOLIS;      // default distribution generator
std::string dstrb_use               = dstrb_all[dstrb_index];

// random number generator
#define RND_MERSENNE_TWISTER    0
#define RND_LEHMER64            1
#define RND_LEHMER128           2
#define RND_WYHASH64            3
//#define RND_PCG64               4
//#define RND_XOSHIRO256          5
//#define RND_SPLIMIX64           6
const std::string rnd_MT         = "MersenneTwister64";  // old, slow, considered not good anymore, from C++11, period 2^19937-1
const std::string rnd_Lehmer64   = "Lehmer64";           // improved 64bit with 64x128bit multiplication, fast, period?
const std::string rnd_Lehmer128  = "Lehmer128";          // improved 64bit with 128x128bit multiplication, fast, period 2^126
const std::string rnd_WyHash64   = "Wyhash64";           // modern, fast, possibly not yet fully understood?
//const std::string rnd_PCG64      = "PCG64";              // modern, fast, possibly not yet fully understood?
//const std::string rnd_xoshiro256 = "xoshiro256**";       // modern, fast, entire xo/r/oshiro family has flaws. this seems still ok. use only for testing!
//const std::string rnd_SplitMix64 = "Splitmix64";         // modern, fast, possibly not yet fully understood?
#if REAL_PRECISION == REAL_DOUBLE
// default random number generators with long period and at least 64bit internal state
#define NUM_RND                 4
#define RND_DEFAULT             RND_LEHMER128            // default random number generator index
std::string rnd_all[NUM_RND]     = {rnd_MT,rnd_Lehmer64,rnd_Lehmer128,rnd_WyHash64
                                    //rnd_PCG64,rnd_xoshiro256,rnd_SplitMix64
                                   };
#elif REAL_PRECISION == REAL_SINGLE
#define RND_LEHMER32            4
#define RND_WELL1024            5
#define RND_LCG_MINST           6
#define RND_RANLUX24            7
#define RND_RANLUX48            8
#define NUM_RND                 9
#define RND_DEFAULT             RND_LEHMER32             // default random number generator index
// old random number generators are not recommended!
const std::string rnd_Lehmer32   = "Lehmer32";           // original 32bit with short period of 2^32 (most likely same as MinSt), slower than Lehmer64
const std::string rnd_WELL       = "WELL1024";           // better than Mersenne Twister but only 32bit
const std::string rnd_MinSt      = "MinStd";             // 32bit linear congruent generator, newer "Minimum standard" from C++11, short period of 2^32
const std::string rnd_Ranlux24   = "Ranlux24";           // supposed to be very good but very slow, from C++11
const std::string rnd_Ranlux48   = "Ranlux48";           // supposed to be very good but very slow, slower than Ranlux24, from C++11
std::string rnd_all[NUM_RND]     = {rnd_MT,rnd_Lehmer64,rnd_Lehmer128,rnd_WyHash64,             
                                    rnd_Lehmer32,rnd_WELL,rnd_MinSt,rnd_Ranlux24,rnd_Ranlux48
                                    //rnd_PCG64,rnd_xoshiro256,rnd_SplitMix64
                                   };
#endif
std::string rnd_use              = rnd_all[RND_DEFAULT]; // name of current selected random number generator

// atoms
uint8_t     atom_num                            = 2;
uint8_t     atom_stat_index [MAX_NUM_SPECIES]   = {STAT_FD, STAT_FD};
std::string atom_stat_use   [MAX_NUM_SPECIES]   = {stat_all[STAT_FD],stat_all[STAT_FD]};
uint8_t     atom_distance_index[MAX_NUM_SPECIES]= {distance_default[STAT_FD], distance_default[STAT_FD]}; // TODO: not used! check this
REAL_TYPE   atom_gamma      [MAX_NUM_SPECIES]   = {0.38, 0.38}; // TODO: not used! check this
std::string atom_name       [MAX_NUM_SPECIES]   = {"Li6","K40"};
REAL_TYPE   atom_mass       [MAX_NUM_SPECIES]   = {6,40}; 
uint64_t    atom_N          [MAX_NUM_SPECIES]   = {25000,4000};
uint64_t    atom_Nc         [MAX_NUM_SPECIES]   = {0,0}; // condensed fraction for STAT_BE
REAL_TYPE   atom_T          [MAX_NUM_SPECIES]   = {300,300};
REAL_TYPE   atom_wx         [MAX_NUM_SPECIES]   = {600,600};
REAL_TYPE   atom_wy         [MAX_NUM_SPECIES]   = {600,600};
REAL_TYPE   atom_wz         [MAX_NUM_SPECIES]   = {80,80};
REAL_TYPE   atom_cx                             = 0;  // displacement in um 2nd species. TODO: not tested. might be used not everywhere.
REAL_TYPE   atom_cy                             = 0;  
REAL_TYPE   atom_cz                             = 0;

// maximum values per species and phase-space dimension + Etot
//REAL_TYPE atom_max[MAX_NUM_SPECIES][DIM_MAXVAL]    = {{0.0,0.0}};
// maximum values for all species. defines volume of cells.
//REAL_TYPE max_x                    [DIM_MAXVAL]    = {0.0};

//Results
REAL_TYPE      atom_Tcrit   [MAX_NUM_SPECIES]   = {0, 0}; // TF or Tc in nK
//REAL_TYPE      atom_R     [MAX_NUM_SPECIES][DIM_PHASESPACE] = {{0,0}}; // in mu
REAL_TYPE      atom_mu      [MAX_NUM_SPECIES]   = {0.0, 0.0}; // in nK
REAL_TYPE      max_E        [MAX_NUM_SPECIES]   = {0.0, 0.0}; // in nK

// molecule
uint8_t     mol_stat_index                      = STAT_BE;
std::string mol_stat_use                        = stat_all[STAT_BE];
std::string mol_name                            = "";

// calculation time in seconds for atoms and molecules
//REAL_TYPE      atom_t          [MAX_NUM_SPECIES]   = {0.0,0.0};
//REAL_TYPE      mol_t                               = 0.0;

// output file precision
#define OUTPUT_PRECISION                        3

// result file
#define FILE_PRECISION                          3
std::string file_result                         = "result.dat";

// files for exporting of atoms and molecules
// when length >0 these are created. these are big files!
std::string export_atom     [MAX_NUM_SPECIES]   = {"", ""};
std::string export_mol                          = "";

// histogram
uint16_t num_bins                               = 300;
// 0 if Histogram should use calculated maximum values.
// 1 if Histogram of atoms should find maximum values.
// for molecules always maximum values are searched.
// in file_result always calculated maximum values are given.
#define HIST_GET_MAX                            0
// if defined export also Epot and Ekin in histogram
#define HIST_EPOT_EKIN
std::string file_histogram                      = "";
// histogram output for each coordinate + total energy
std::string cname[DIM_MAXVAL]                   = {"x","y","z","vx","vy","vz","Etot"};
std::string cunit[DIM_MAXVAL]                   = {"um","um","um","mm/s","mm/s","mm/s","nK"};

// if != 0 find nearest neighbor (slower), otherwise stop molecule search at first pair
uint8_t find_nearest                            = 0;

// when neighbors_file given
// for each species[0] gets neighbors_num nearest neighbors of pecies[1] 
std::string neighbors_file                      = "";
uint8_t     neighbors_num                       = 1;

// number of repetitions
uint32_t    repetitions                         = 1;

// energy cutoff per species
REAL_TYPE CalculationSize   [MAX_NUM_SPECIES]   = {1.0, 1.0};

// calculation of chemical potential
REAL_TYPE ErrorAllowed                          = 1e-10;    // max. relative error of atom number to calculate chemical potential
REAL_TYPE error_mu                              = 1e-3;     // max error of mu in nK
REAL_TYPE InitialStepSize                       = 0.1;      // initial step size*E_Fermi to calculate chemical potential
uint64_t ESteps             [MAX_NUM_SPECIES]   = {100000,100000}; // TODO: check if reasonable

// number of helper threads
uint8_t num_threads                             = 8;

// thread synchronization
MUTEX   global_mutex;
BARRIER global_barrier;

// console handle for coloring text output (only windows)
GET_STDOUT_HANDLE;
GET_STDERR_HANDLE;

// parameter file names, types and pointers
#define NUM_VAR                 41
#define TYPE_U8                 0
#define TYPE_U16                1
#define TYPE_U32                2
#define TYPE_U64                3
#define TYPE_DOUBLE             8
#define TYPE_STRING             9
#define TYPE_MASK               15
#define FLAG_SP1                16  // set for 2nd species. used to check if atom_num=1 and 2nd species selected.
#define FLAG_SET                32  // set automatically after parameter has been read
struct params {
    std::string name;
    void *ptr;
    unsigned char flag;     // TYPE | FLAG
};
struct params var_params[NUM_VAR] = {
    {"species_0",                   &atom_name[0],          TYPE_STRING             },
    {"species_1",                   &atom_name[1],          TYPE_STRING | FLAG_SP1  },
    {"num_species",                 &atom_num,              TYPE_U8                 },
    {"statistics_0",                &atom_stat_use[0],      TYPE_STRING             },
    {"statistics_1",                &atom_stat_use[1],      TYPE_STRING | FLAG_SP1  },
    {"statistics_mol",              &mol_stat_use,          TYPE_STRING             },
    {"mass_0",                      &atom_mass[0],          TYPE_DOUBLE             },
    {"mass_1",                      &atom_mass[1],          TYPE_DOUBLE | FLAG_SP1  },
    {"N_0",                         &atom_N[0],             TYPE_U32                },
    {"N_1",                         &atom_N[1],             TYPE_U32    | FLAG_SP1  },
    {"T_0",                         &atom_T[0],             TYPE_DOUBLE             },
    {"T_1",                         &atom_T[1],             TYPE_DOUBLE | FLAG_SP1  },
    {"fx_0",                        &atom_wx[0],            TYPE_DOUBLE             },
    {"fy_0",                        &atom_wy[0],            TYPE_DOUBLE             },
    {"fz_0",                        &atom_wz[0],            TYPE_DOUBLE             },
    {"fx_1",                        &atom_wx[1],            TYPE_DOUBLE | FLAG_SP1  },
    {"fy_1",                        &atom_wy[1],            TYPE_DOUBLE | FLAG_SP1  },
    {"fz_1",                        &atom_wz[1],            TYPE_DOUBLE | FLAG_SP1  },
    {"cx",                          &atom_cx,               TYPE_DOUBLE             },
    {"cy",                          &atom_cy,               TYPE_DOUBLE             },
    {"cz",                          &atom_cz,               TYPE_DOUBLE             },
    {"size_0",                      &CalculationSize[0],    TYPE_DOUBLE             },
    {"size_1",                      &CalculationSize[1],    TYPE_DOUBLE | FLAG_SP1  },
    {"error_allowed",               &ErrorAllowed,          TYPE_DOUBLE             },
    {"step_size",                   &InitialStepSize,       TYPE_DOUBLE             },
    {"repeat",                      &repetitions,           TYPE_U32                },
    {"distance_measure",            &mol_distance_use,      TYPE_STRING             },
    {"gamma",                       &mol_gamma,             TYPE_DOUBLE             },
    {"threads",                     &num_threads,           TYPE_U8                 },
    {"num_bins",                    &num_bins,              TYPE_U16                },
    {"random_number_generator",     &rnd_use,               TYPE_STRING             },
    {"distribution_generator",      &dstrb_use,             TYPE_STRING             },
    {"seed",                        &seed_use,              TYPE_STRING             },
    {"result_file",                 &file_result,           TYPE_STRING             },
    {"histogram_file",              &file_histogram,        TYPE_STRING             },
    {"export_atom_0",               &export_atom[0],        TYPE_STRING             },
    {"export_atom_1",               &export_atom[1],        TYPE_STRING | FLAG_SP1  },
    {"export_molecule",             &export_mol,            TYPE_STRING             },    
    {"neighbors_file",              &neighbors_file,        TYPE_STRING             },    
    {"neighbors_num",               &neighbors_num,         TYPE_U8                 },
    {"find_nearest",                &find_nearest,          TYPE_U8                 }
//    {"neighbors_distance_measure",  &neighbors_distance_use,TYPE_STRING             },    
//    {"neighbors_species_0",         &neighbors_species[0],  TYPE_U8                 },    
//    {"neighbors_species_1",         &neighbors_species[1],  TYPE_U8                 },
//    {"neighbors_bothways",          &neighbors_bothways,    TYPE_U8                 }
};

// species data
species_data::species_data() {
    //prev = next = NULL;
    mass = 0.0;
    for (int i = 0; i < DIM_PHASESPACE; ++i) x[i] = 0.0;
    for (int i = 0; i < DIM_ENERGY    ; ++i) E[i] = 0.0;
    for (int i = 0; i < DIM_PHASESPACE; ++i) index[i] = 0;
    SET_STATUS_FREE(status);
#ifdef _DEBUG
    distance = 0.0;
    atom_index[0] = atom_index[1] = (uint64_t)-1;
#endif
}

species_data::~species_data() {
}

/* species list
species_entry::species_entry(species_data *species0, species_data *species1, uint64_t hash_key) {
    this->prev = this->next = NULL;
    this->species[0] = species0;
    this->species[1] = species1;
    this->hash_key = hash_key;
}*/

species_entry::species_entry() {
    this->prev = this->next = NULL;
    this->hash_key = -1;
}

species_entry::~species_entry() {
    //note: we must delete species manually since one species might be in several cells list!
}

/* returns next species_entry entry with distance^2 <= distance_squared starting from first_entry (inclusive)
// returns next species_entry entry with updated distance_squared.
// search includes always first atom, i.e. call first->get_next() if first should be skipped.
// if distance_squared == DBL_MAX returns first atom with distance_squared updated.
// if first_entry == NULL returns NULL
// note: 
// - species can be different: 'this' and 'first_entry' can be from different lists!
//   i.e. they can be species from different cells or even different species with different mass.
class species_entry * species_entry::get_next_distance(
    class species_entry *first_entry, 
    REAL_TYPE              *distance_squared, 
    uint8_t              distance_index
    ) {
    class species_entry *next_entry = first_entry;
    //class species_data *this_species = this->species[0], *next;
    class species_data *this_species = &this->species[0], *next;
    REAL_TYPE d = *distance_squared;
    if (distance_index == DISTANCE_X) {
        // real space distance delta x
        REAL_TYPE r, tmp;
        REAL_TYPE x = this_species->x[0];
        REAL_TYPE y = this_species->x[1];
        REAL_TYPE z = this_species->x[2];
        while (next_entry) {
            next = &next_entry->species[0];
            tmp = next->x[0] - x;
            r = tmp*tmp;
            tmp = next->x[1] - y;
            r += tmp*tmp;
            tmp = next->x[2] - z;
            r += tmp*tmp;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_V) {
        // velocity distance delta v
        REAL_TYPE r, tmp;
        REAL_TYPE vx = this_species->x[3];
        REAL_TYPE vy = this_species->x[4];
        REAL_TYPE vz = this_species->x[5];
        while (next_entry) {
            next = &next_entry->species[0];
            tmp = next->x[3] - vx;
            r = tmp*tmp;
            tmp = next->x[4] - vy;
            r += tmp*tmp;
            tmp = next->x[5] - vz;
            r += tmp*tmp;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_PV) {
        // velocity-space distance delta x * delta v
        // want to substract velocities in center of mass frame
        // this is equivalent to simply substracting the velocities.
        // to convert relative velocity to momentum in center of mass frame, use reduced mass
        REAL_TYPE r, v, tmp;
        REAL_TYPE x  = this_species->x[0];
        REAL_TYPE y  = this_species->x[1];
        REAL_TYPE z  = this_species->x[2];
        REAL_TYPE vx = this_species->x[3];
        REAL_TYPE vy = this_species->x[4];
        REAL_TYPE vz = this_species->x[5];
        while (next_entry) {
            next = &next_entry->species[0];
            // (mu)^2
            tmp = next->x[0] - x;
            r = tmp*tmp;
            tmp = next->x[1] - y;
            r += tmp*tmp;
            tmp = next->x[2] - z;
            r += tmp*tmp;
            // (mm/s)^2
            tmp = next->x[3] - vx;
            v = tmp*tmp;
            tmp = next->x[4] - vy;
            v += tmp*tmp;
            tmp = next->x[5] - vz;
            v += tmp*tmp;
            // um^2*(mm/s)^2 = (m^2/s*1e-9)^2
            // [dx*dp] = kg*m^2/s = Js
            v *= r;
            if (v <= d) {
                *distance_squared = v;
                break;
            }            
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_PP) {
        // phase-space distance delta x * delta p = lab frame momentum-space distance
        REAL_TYPE r, p, tmp;
        REAL_TYPE x  = this_species->x[0];
        REAL_TYPE y  = this_species->x[1];
        REAL_TYPE z  = this_species->x[2];
        REAL_TYPE px = this_species->x[3]*this_species->mass;
        REAL_TYPE py = this_species->x[4]*this_species->mass;
        REAL_TYPE pz = this_species->x[5]*this_species->mass;
        while (next_entry) {
            next = &next_entry->species[0];
            // (mu)^2
            tmp = next->x[0] - x;
            r = tmp*tmp;
            tmp = next->x[1] - y;
            r += tmp*tmp;
            tmp = next->x[2] - z;
            r += tmp*tmp;
            // (mm/s * amu)^2
            tmp = next->x[3]*next->mass - px;
            p = tmp*tmp;

            tmp = next->x[4]*next->mass - py;
            p += tmp*tmp;
            tmp = next->x[5]*next->mass - pz;
            p += tmp*tmp;

            // [dx*dp] = kg*m^2/s = Js
            p *= r;
            if (p <= d) {
                *distance_squared = p;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_CROSS) {
        // cross momentum-space distance: |delta x  x  delta p|
        REAL_TYPE rx, ry, rz, p;
        REAL_TYPE x  = this_species->x[0];
        REAL_TYPE y  = this_species->x[1];
        REAL_TYPE z  = this_species->x[2];
        REAL_TYPE px = this_species->x[3]*this_species->mass;
        REAL_TYPE py = this_species->x[4]*this_species->mass;
        REAL_TYPE pz = this_species->x[5]*this_species->mass;
        while (next_entry) {
            next = &next_entry->species[0];
            // (mm/s * amu) * (mu)
            rx =  (next->x[1] - y)*(next->x[5] - pz) - (next->x[2] - z)*(next->x[4] - py);
            ry = -(next->x[0] - x)*(next->x[5] - pz) + (next->x[2] - z)*(next->x[3] - px);
            rz =  (next->x[0] - x)*(next->x[4] - py) - (next->x[1] - y)*(next->x[3] - px);
            // (kg*m^2/s)^2 = (Js)^2
            p  = rx*rx + ry*ry + rz*rz;
            if (p <= d) {
                *distance_squared = p;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_MEAN) {
        // mean momentum-space distance: <|delta x| * |delta p|>
        REAL_TYPE rx, ry, rz, p;
        REAL_TYPE x  = this_species->x[0];
        REAL_TYPE y  = this_species->x[1];
        REAL_TYPE z  = this_species->x[2];
        REAL_TYPE px = this_species->x[3]*this_species->mass;
        REAL_TYPE py = this_species->x[4]*this_species->mass;
        REAL_TYPE pz = this_species->x[5]*this_species->mass;
        while (next_entry) {
            next = &next_entry->species[0];
            // (mm/s * amu) * (mu)
            rx = (next->x[0] - x)*(next->x[3] - px);
            ry = (next->x[1] - y)*(next->x[4] - py);
            rz = (next->x[2] - z)*(next->x[5] - pz);
            // (kg*m^2/s)^3 = (Js)^3
            p  = abs(rx * ry * rz);
            if (p <= d) {
                *distance_squared = p;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else if (distance_index == DISTANCE_MAX) {
        // maximum momentum-space distance: max|delta xi * delta pi|
        REAL_TYPE rx, ry, rz, p;
        REAL_TYPE x  = this_species->x[0];
        REAL_TYPE y  = this_species->x[1];
        REAL_TYPE z  = this_species->x[2];
        REAL_TYPE px = this_species->x[3]*this_species->mass;
        REAL_TYPE py = this_species->x[4]*this_species->mass;
        REAL_TYPE pz = this_species->x[5]*this_species->mass;
        while (next_entry) {
            next = &next_entry->species[0];
            // (mm/s * amu) * (mu)
            rx = (next->x[0] - x)*(next->x[3] - px);
            ry = (next->x[1] - y)*(next->x[4] - py);
            rz = (next->x[2] - z)*(next->x[5] - pz);
            // (kg*m^2/s)^2 = (Js)^2
            rx *= rx;
            ry *= ry;
            rz *= rz;
            // take maximum (delta_xi*delta_pi)^2
            if (rx >= ry) {
                p = (rx >= rz) ? rx : rz;
            }
            else {
                p = (ry >= rz) ? ry : rz;
            }
            if (p <= d) {
                *distance_squared = p;
                break;
            }
            next_entry = next_entry->get_next();
        }
    }
    else {
        // unknown distance measure
        std::cout << "get_next_distance error: distance measure ' " << distance_all[distance_index] << " ' (" << distance_index << ") not implemented!" << std::endl;
        next_entry = NULL; 
    }
    
    return next_entry;
}

// returns next cell with distance^2 <= distance_squared starting from first
// returns next cell with updated distance_squared.
// if distance_squared == DBL_MAX returns first cell with distance_squared updated.
// TODO: 
// - use neighboring cells instead of plain next cells.
// - maybe move into cells_list instead of species_entry?
// - at the moment uses species_entry::x which is available only during debugging.
// - not all distance measures implemented
class cell_data * species_entry::get_next_distance(
    class cell_data *first, 
    REAL_TYPE *distance_squared, 
    uint8_t distance_index,
    const REAL_TYPE *dx2,
    REAL_TYPE mass2
    ) {
    class cell_data *next = first;
    REAL_TYPE r, d_minus, d_plus, d;
    class species_data * this_species = this->get_species(0);
    REAL_TYPE x  = this_species->x[0];
    REAL_TYPE y  = this_species->x[1];
    REAL_TYPE z  = this_species->x[2];
    REAL_TYPE vx = this_species->x[3];
    REAL_TYPE vy = this_species->x[4];
    REAL_TYPE vz = this_species->x[5];
    d = *distance_squared;
    if (distance_index == DISTANCE_X) {
        // real space distance delta x
        while (next) {
            d_minus = d_plus = next->x[0] - x;
            d_minus -= dx2[0];
            d_plus  += dx2[0];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[1] - y;
            d_minus -= dx2[1];
            d_plus  += dx2[1];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[2] - z;
            d_minus -= dx2[2];
            d_plus  += dx2[2];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next = next->get_next();
        }
    }
    else if (distance_index == DISTANCE_V) {
        // velocity distance delta v
        while (next) {
            d_minus = d_plus = next->x[3] - vx;
            d_minus -= dx2[3];
            d_plus  += dx2[3];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[4] - vy;
            d_minus -= dx2[4];
            d_plus  += dx2[4];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[5] - vz;
            d_minus -= dx2[5];
            d_plus  += dx2[5];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next = next->get_next();
        }
    }
    else if (distance_index == DISTANCE_PV) {
        // velocity-space distance delta x * delta v
        while (next) {
            d_minus = d_plus = next->x[0] - x;
            d_minus -= dx2[0];
            d_plus  += dx2[0];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[1] - y;
            d_minus -= dx2[1];
            d_plus  += dx2[1];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[2] - z;
            d_minus -= dx2[2];
            d_plus  += dx2[2];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[3] - vx;
            d_minus -= dx2[3];
            d_plus  += dx2[3];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            REAL_TYPE p = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[4] - vy;
            d_minus -= dx2[4];
            d_plus  += dx2[4];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            p += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[5] - vz;
            d_minus -= dx2[5];
            d_plus  += dx2[5];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            p += (d_minus < d_plus) ? d_minus : d_plus;
            r *= p;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next = next->get_next();
        }
    }
    else if (distance_index == DISTANCE_PP) {
        // phase-space distance delta x * delta p = lab frame momentum-space distance
        vx *= this_species->mass;
        vy *= this_species->mass;
        vz *= this_species->mass;
        while (next) {
            d_minus = d_plus = next->x[0] - x;
            d_minus -= dx2[0];
            d_plus  += dx2[0];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[1] - y;
            d_minus -= dx2[1];
            d_plus  += dx2[1];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[2] - z;
            d_minus -= dx2[2];
            d_plus  += dx2[2];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            r += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[3]*mass2 - vx;
            d_minus -= dx2[3];
            d_plus  += dx2[3];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            REAL_TYPE p = (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[4]*mass2 - vy;
            d_minus -= dx2[4];
            d_plus  += dx2[4];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            p += (d_minus < d_plus) ? d_minus : d_plus;
            d_minus = d_plus = next->x[5]*mass2 - vz;
            d_minus -= dx2[5];
            d_plus  += dx2[5];
            d_minus *= d_minus;
            d_plus  *= d_plus;
            p += (d_minus < d_plus) ? d_minus : d_plus;
            r *= p;
            if (r <= d) {
                *distance_squared = r;
                break;
            }
            next = next->get_next();
        }
    }
    else if (distance_index == DISTANCE_CROSS) {
        // cross momentum-space distance: |delta x  x  delta p|
    }
    else if (distance_index == DISTANCE_MEAN) {
        // mean momentum-space distance: <|delta x| * |delta p|>
    }
    else if (distance_index == DISTANCE_MAX) {
        // maximum momentum-space distance: max|delta xi * delta pi|
    }
    return next;
}*/


/* returns next atom and cell with distance^2 <= distance_squared starting from input cell (inclusive).
// distance = min(|this->hash_key - cell->d_min|, |this->hash_key - cell->d_max|).
// returns NULL when no further cell found and output cell == NULL, or when input cell == NULL.
// atoms for cell must have been created using the same distance measure as to generate this->hash_key.
class species_entry * species_entry::get_next_distance(class species_entry *entry, REAL_TYPE *distance_squared) {
    REAL_TYPE d0 = this->hash_key;
    REAL_TYPE d2 = (*distance_squared == 0.0) ? DBL_MAX : *distance_squared;
    
    while (entry) {
        REAL_TYPE d = entry->hash_key - d0;
        d *= d;
        if (d <= d2) {
            *distance_squared = d;
            break;
        }
        entry = entry->get_next();
    }

    return entry;
}

// same as before but does not return new distance^2
class species_entry * species_entry::get_next_distance(class species_entry *entry, REAL_TYPE distance_squared) {
    REAL_TYPE d0 = this->hash_key;
    if (distance_squared == 0.0) distance_squared = DBL_MAX;
    
    while (entry) {
        REAL_TYPE d = entry->hash_key - d0;
        if ((d*d) <= distance_squared) break;
        entry = entry->get_next();
    }

    return entry;
}

// returns next cell with distance^2 <= distance_squared starting from input cell (inclusive).
// distance = min(|this->hash_key - cell->d_min|, |this->hash_key - cell->d_max|).
// returns NULL when no further cell found and output cell == NULL, or when input cell == NULL.
// atoms for cell must have been created using the same distance measure as to generate this->hash_key.
// TODO: this is wrong now since hash_key = uint_64(distance/delta)
class cell_data * species_entry::get_next_distance(class cell_data *cell, REAL_TYPE distance_squared) {
    uint64_t d0 = this->hash_key;
    REAL_TYPE delta = 0;
    if (distance_squared == 0.0) distance_squared  = DBL_MAX;
    
    while (cell) {
        uint64_t hash_key = cell->get_hash_key();
        REAL_TYPE d_min      = d0 - (hash_key - delta/2);
        d_min = d_min*d_min;
        REAL_TYPE d_max    = d0 - (hash_key + delta/2);
        d_max = d_max*d_max;
        if (d_max < d_min) d_min = d_max;
        if (d_min <= distance_squared) break;
        cell = cell->get_next();
    }

    return cell;
}*/

// note: pthreads must be the first thread or NULL
cells_list::cells_list(
        struct thread_data *pthreads,
        uint8_t             type, 
        uint64_t            num_cells, 
        uint64_t            num_atoms,
        uint8_t             index0, 
        uint8_t             index1, 
        std::string         species_name, 
        uint8_t             statistics,
        uint8_t             distance_index, 
        REAL_TYPE              gamma, 
        //bool                duplicate,
        REAL_TYPE              *max,
        REAL_TYPE              *R
        ) {

#if NUM_CELLS > 1
    if (type == CELLS_TYPE_ATOMS) {
        // get maximum distance^2 for given species and distance index
        if (distance_index == DISTANCE_X) {
            this->dim    = 3;
            this->offset = 0;
            //mult   = 1.0;
        }
        else if (distance_index == DISTANCE_V) {
            this->dim    = 3;
            this->offset = 3;
            //mult   = 1.0; 
        }
        else if (distance_index == DISTANCE_PV) {
            this->dim    = 6;
            this->offset = 0;
            //mult   = 1.0;
        }
        else if (distance_index == DISTANCE_PP) {
            this->dim    = 6;
            this->offset = 0;
            //mult   = ::atom_mass[species];
        }
        else if (distance_index == DISTANCE_CROSS) {
            this->dim    = 6;
            this->offset = 0;
            //mult   = ::atom_mass[species];
        }
        else if (distance_index == DISTANCE_MEAN) {
            this->dim    = 6;
            this->offset = 0;
            //mult   = ::atom_mass[species];
        }        
        else if (distance_index == DISTANCE_MAX) {
            this->dim    = 6;
            this->offset = 0;
            //mult   = ::atom_mass[species];
        }
        else {
            // unknown distance measure
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << "cells_list init error: distance measure ' " << distance_all[distance_index] << " ' (" << distance_index << ") not implemented!" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            return;
        }
        this->x_min = new REAL_TYPE[dim];
        this->dx    = new REAL_TYPE[dim];
        this->cells_per_dim = (uint8_t)(std::floor(std::pow(num_cells,1.0/dim)));
        if (this->cells_per_dim == 0) {
            this->cells_per_dim = 1;
            std::cout << "cells_list init: num_cells " << num_cells << " gives 0 cells/dim with dim = " << dim << "! increase to 1 cells/dim." << std::endl;
        }
        for (uint8_t i = 0; i < dim; ++i) {
            this->dx[i] = 2 * max[i + offset] / cells_per_dim; // TODO: have removed -1 from cells_per_dim since seems wrong!
            this->x_min[i] = -max[i + offset];
        }
    }
    else if (type == CELLS_TYPE_MOLECULES) {
        // molecule requires 1d list of num_thread cells and hash_key = thread id 
        //this->min        = 0.0;
        //this->max        = 0.0;
        //this->delta      = 0.0;
        this->dim           = 1;
        this->offset        = 0;
        this->cells_per_dim = (uint8_t)num_cells; // TODO: limits to 256 cells!
        this->x_min = new REAL_TYPE[1];
        this->dx    = new REAL_TYPE[1];
        x_min[0] = 0.0;
        dx  [0] = 1.0;
    }
    else {
        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
        std::cout << "cells_list init error: unknown type ' " << type << " '!" << std::endl;
        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
        return;
    }
#else
    // for NUM_CELLS == 1 we create only 1d list of num_thread cells and hash_key = thread id
    //this->min        = 0.0;
    //this->max        = 0.0;
    //this->delta      = 0.0;
    this->dim = 1;
    this->offset = 0;
    this->cells_per_dim = (uint8_t)num_cells;
    this->x_min = new REAL_TYPE[1];
    this->dx = new REAL_TYPE[1];
    x_min[0] = 0.0;
    dx[0] = 1.0;
#endif

    this->type           = type;
    this->index[0]       = index0;
    this->index[1]       = index1;
    this->statistics     = statistics;
    this->distance_index = distance_index;
    this->gamma          = gamma;
    this->species_name   = species_name;
    this->num_atoms      = num_atoms;
    //this->duplicate      = duplicate;
    
    // init max values and size. can be NULL
    this->max = max;
    this->R   = R;
    
    // if owner_thread given and number of cells given, create cells with hash_key = linear cell index
    if ((pthreads != NULL) && (num_cells > 0)) {
        uint64_t i = 0;
        uint8_t ii[DIM_PHASESPACE];
        REAL_TYPE  xc[DIM_PHASESPACE];
        if (dim == 1) {
            xc[0] = x_min[0];
            for (ii[0] = 0; ii[0] < cells_per_dim; ++ii[0], ++i) {
                this->add(new class cell_data(&pthreads[i % num_threads], dim, i, ii, x_min));
                xc[0] += dx[0];
            }
        }
        else if (dim == 3) {
            xc[0] = x_min[0];
            for (ii[0] = 0; ii[0] < cells_per_dim; ++ii[0]) {
                xc[1] = x_min[1];
                for (ii[1] = 0; ii[1] < cells_per_dim; ++ii[1]) {
                    xc[2] = x_min[2];
                    for (ii[2] = 0; ii[2] < cells_per_dim; ++ii[2], ++i) {
                        this->add(new class cell_data(&pthreads[i % num_threads], dim, i, ii, x_min));
                        xc[2] += dx[2];
                    }
                    xc[1] += dx[1];
                }
                xc[0] += dx[0];
            }
        }
        else if (dim == 6) {
            xc[0] = x_min[0];
            for (ii[0] = 0; ii[0] < cells_per_dim; ++ii[0]) {
                xc[1] = x_min[1];
                for (ii[1] = 0; ii[1] < cells_per_dim; ++ii[1]) {
                    xc[2] = x_min[2];
                    for (ii[2] = 0; ii[2] < cells_per_dim; ++ii[2]) {
                        xc[3] = x_min[3];
                        for (ii[3] = 0; ii[3] < cells_per_dim; ++ii[3]) {
                            xc[4] = x_min[4];
                            for (ii[4] = 0; ii[4] < cells_per_dim; ++ii[4]) {
                                xc[5] = x_min[5];
                                for (ii[5] = 0; ii[5] < cells_per_dim; ++ii[5], ++i) {
                                    this->add(new class cell_data(&pthreads[i % num_threads], dim, i, ii, xc));
                                    xc[5] += dx[5];
                                }
                                xc[4] += dx[4];
                            }
                            xc[3] += dx[3];
                        }
                        xc[2] += dx[2];
                    }
                    xc[1] += dx[1];
                }
                xc[0] += dx[0];
            }
        }
    }
} 

cells_list::~cells_list() {
    // delete all cells
    reset(true);
    delete [] x_min;
    delete [] dx;
    if (max) delete [] max;
    if (R)   delete [] R;
}

// delete species and if delete_all delete cells
// note: this does not reset any of the cell properties like num_atoms
int cells_list::reset(bool delete_all) { 
    int error = 0;
    class cell_data *cell = get_first(), *next_cell;
    class species_entry *entry, *next_entry;
    uint64_t count_cells = 0, count_atoms_tot = 0;
#ifdef _DEBUG
    uint64_t num_cells = this->get_num();
#endif
    while (cell) {
        // save next cell
        next_cell = cell->get_next();
        
        // remove all entries in list. 
        // update: duplicate is removed and never delete atoms here, only remove entries from list!
        uint64_t num = cell->species.get_num(), count_atoms = 0;
        entry = cell->species.get_first();
        while (entry) {
            ++count_atoms;
            next_entry = entry->get_next();
            cell->species.remove(entry);
            /*if (!duplicate) {
                // delete species index 0
                class species_data *sp;
                sp = entry->get_species(0);
                if (sp) delete sp;
                // delete species index 1
                sp = entry->get_species(1);
                if (sp) delete sp;
            }*/
            // we always delete entry - update: only thread is allowed to delete its entries 
            //delete entry;
            entry = next_entry;
        }
#ifdef _DEBUG
        if ( (cell->species.get_num() != 0) || (cell->species.get_first() != NULL) || (cell->species.get_last() != NULL) ) {
            std::cout << std::endl << "cells_list::reset " << species_name << ": atom list not empty! remaining " << cell->species.get_num() << " first " << cell->species.get_first() << ", last " << cell->species.get_last() << std::endl << std::endl;
            error = -1;
            break;
        }
        else if (count_atoms != num) {
            std::cout << std::endl << "cells_list::reset " << species_name << ": atoms list corrupt! " << count_atoms << " in list but " << num << " expected!" << std::endl << std::endl;
            error = -2;
            break;
        } 
#endif
        if (delete_all) {
            // remove cell from list and delete
            remove(cell);
            delete cell;
        }

        cell = next_cell;
        ++count_cells;
        count_atoms_tot += num;
    }
#ifdef _DEBUG
    if ((!delete_all) && (count_atoms_tot != this->num_atoms) && (this->type == CELLS_TYPE_ATOMS)) {
        // this code has a problem since num_atoms is not reset to actual number of atoms.
        // when reset_all is used then cells were already reset before and count_atoms_tot = 0 while num_atoms > 0
        std::cout << std::endl << "cells_list::reset " << species_name << ": atoms list corrupt! total " << count_atoms_tot << " in list but " << this->num_atoms << " expected!" << std::endl << std::endl;
        error = -3;
    }
    else if (count_cells != num_cells) {
        std::cout << std::endl << "cells_list::reset " << species_name << ": cells list corrupt! total " << count_cells << " in list but " << num_cells << " expected!" << std::endl << std::endl;
        error = -4;
    } 
    else if (delete_all && (this->get_num() || this->get_first() || this->get_last())) {
        std::cout << std::endl << "cells_list::reset " << species_name << ": cells list not empty! remaining " << this->get_num() << " first " << this->get_first() << ", last " << this->get_last() << std::endl << std::endl;
        error = -4;
    }
#endif
    // TODO: if not delete_all, we do not reset num_atoms here because otherwise we have to recalculate it each repetition.
    //       this is inconsistent with the actual number of atoms in list. 
    if (delete_all) {
        this->num_atoms = 0;
    }
    
    return error;
}


// return hash_key for the given species using distance measure index of cells_list
uint64_t cells_list::get_hash_key(class species_data *species) {
    uint64_t hash_key = 0;
    for (uint8_t i = 0; i < dim; ++i) {
        hash_key = hash_key*cells_per_dim + (uint64_t)((species->x[i+offset]-x_min[i])/dx[i]);
    }
    return hash_key;
}

// cell data
cell_data::cell_data (
    struct thread_data *owner, 
    uint8_t dim, 
    uint64_t hash_key 
//#ifdef _DEBUG
    ,uint8_t *index, REAL_TYPE *x
//#endif
    ) {
    owner_thread = owner;
    prev = next = NULL;
    this->hash_key = hash_key;
    this->nb = new class cell_data*[dim*2];
//#ifdef _DEBUG    
    // 3/6 cell index
    this->index = new uint8_t[dim];
    // 3/6 center coordinates of cell
    this->x = new REAL_TYPE[dim];
    for (uint8_t i = 0; i < dim; ++i) {
        this->index[i] = index[i];
        this->x    [i] = x[i];
        nb[(2*i)+0] = NULL;
        nb[(2*i)+1] = NULL;
    }
//#else
//    for (uint8_t i = 0; i < (dim*2); ++i) {
//        nb[i] = NULL;
//    }
//#endif
    // TODO: assign neighbors (see old code)
}

cell_data::~cell_data() {
    delete [] nb;
#ifdef _DEBUG    
    delete [] index;
    delete [] x;
#endif
}

struct thread_data * cell_data::get_owner(void) { 
    return owner_thread;     
}

uint8_t cell_data::get_owner_id(void) { 
    return owner_thread->id;
}

/* returns coordinates index from absolute index
void cell_data::_get_coords(uint64_t index, uint8_t _index[DIM_CELLS]) {
    for (int d = 0; d < DIM_CELLS; ++d) {
        _index[d] =             index % NUM_CELLS_PER_DIM;
        index     = (uint64_t) (index / NUM_CELLS_PER_DIM);
    }    
}*/

/* returns abolute index of this cell
inline uint64_t cell_data::get_index(void) {
    uint64_t ii = 0;
    for (int8_t d = DIM_CELLS-1; d >= 0; --d) {
        ii = ii*NUM_CELLS_PER_DIM + (uint64_t)this->index[d];
    }
    //printf("get_index [%u,%u,%u] -> %lu\n", +index[0], +index[1], +index[2], ii);
    return ii;
}

// returns absolute index of given dimensional index
// for DIM_CELLS=3 returns ((index[offset+2]*DIM_CELLS_PER_DIM + index[offset+1])*DIM_CELLS_PER_DIM + index[offset])
inline uint64_t cell_data::get_index(uint8_t index[DIM_PHASESPACE], uint8_t offset) {
    uint64_t ii = 0;
    for (uint8_t d = offset+DIM_CELLS-1; ; --d) {
        ii = ii*NUM_CELLS_PER_DIM + (uint64_t)index[d];
        if (d == offset) break;
    }
    //printf("_get_index [%li,%li,%li] -> %li\n", _index[0], _index[1], _index[2], ii);
    return ii;
}*/

// lock and unlock cell
void cell_data::lock(void) { 
    MUTEX_LOCK(owner_thread->mutex); 
}
void cell_data::unlock(void) { 
    MUTEX_UNLOCK(owner_thread->mutex); 
}

/* returns cell of given relative index starting from this cell
// use to find neighbouring cells with _index = cells_offset[#]
class cell_data * cell_data::find_cell(int16_t _index[DIM_CELLS]) {
    class cell_data *cell = this;

    // difference of linear index of cell - this
    int64_t ii = 0, i;
    for (int8_t d = DIM_CELLS-1; d >= 0; --d) {
        i = this->index[d] + _index[d];
        if ((i < 0) || (i >= NUM_CELLS_PER_DIM)) return NULL; // cell coordinate out of range  
        ii = ii*NUM_CELLS_PER_DIM + _index[d];
    }
    //std::cout << "find cell [" << _index[0] << "," << _index[1] << "," << _index[2] << "] from [" << +index[0] << "," << +index[1] << "," << +index[2] << "]: delta = " << ii << std::endl; sleep_ms(10);

    if (ii >= 0) {

        // walk though list towards end
        for (i = 0; (i < ii) && (cell != NULL); ++i) cell = cell->next;

    }

    else {

        // walk though list towards beginning
        for (i = 0; (i > ii) && (cell != NULL); --i) cell = cell->prev;

    }
    
#ifdef _DEBUG

    // check that the found cell is the correct one
    if ( cell ) {
        if ((this->get_index() + ii) != cell->get_index()) {
            printf("\nerror find cell!\n");
            printf("this  cell %8lu [%4u,%4u,%4u]\n", this->get_index(), this->index[0], this->index[1], this->index[2]);
            printf("delta cell %8li [%4u,%4u,%4u]\n", ii, _index[0], _index[1], _index[2]);
            printf("found cell %8lu [%4u,%4u,%4u]\n", cell->get_index(), cell->index[0], cell->index[1], cell->index[2]);
            printf("\n");
            cell = NULL;
        }
        else {
            for (uint8_t d = 0; d < DIM_CELLS; ++d) { 
                if (cell->index[d] != index[d] + _index[d]) {
                    printf("\nerror find cell!\n");
                    printf("this  cell %8lu [%4u,%4u,%4u]\n", this->get_index(), this->index[0], this->index[1], this->index[2]);
                    printf("delta cell %8li [%4u,%4u,%4u]\n", ii, _index[0], _index[1], _index[2]);
                    printf("found cell %8lu [%4u,%4u,%4u]\n", cell->get_index(), cell->index[0], cell->index[1], cell->index[2]);
                    printf("\n");
                    cell = NULL;
                 }
            }
        }
    }

#endif    
    return cell;
}*/

// expands filename with name and number
// name can be NULL
std::string expand_filename(std::string &filename, std::string& name, int number) {
    std::stringstream str_out;
    size_t pos = filename.rfind('.');
    if (pos == std::string::npos) { // no file extension
        if (name.length() != 0) str_out << filename << '_' << name << number;
        else                    str_out << filename << '_' << number;
    }
    else {
        if (name.length() != 0) str_out << filename.substr(0, pos) << '_' << name << '_' << number << filename.substr(pos, filename.length() - pos);
        else                    str_out << filename.substr(0, pos) << '_' << number << filename.substr(pos, filename.length() - pos);
    }
    return str_out.str();
}

// export atoms or molecules coordinates from given cells index
// returns 0 on success.
// note: execute only on main thread.
int ExportAtoms(std::string &filename, uint8_t index) {
    int error = 0;
    std::ofstream out;
    uint32_t t_start = GetTickCount();
    uint32_t tot_atoms = 0, num;
    uint64_t _free = 0, paired = 0; //, unpaired = 0;
    std::string &species_name = ::cells[index]->get_species_name();
    uint8_t distance_index    = ::cells[index]->get_distance_index(); 
    class cell_data *cell     = ::cells[index]->get_first();
    bool is_atom              = (::cells[index]->get_type() == CELLS_TYPE_ATOMS);
    
    out.open(filename);
    if (out.fail()) {
        std::cout << "could not open export file '" << filename << "'!" << std::endl;
        error = -1;
    }
    else {
        if (is_atom) {
            std::cout << "exporting atoms file '" << filename << "' ..." << std::endl;
            out << "index,x(mu),y(mu),z(mu),vx(mm/s),vy(mm/s),vz(mm/s),status,Epot(nK),Ekin(nK),Etot(nK)" << std::endl;
        }
        else {
            std::cout << "exporting molecules file '" << filename << "' ..." << std::endl;
            out << "index,x(mu),y(mu),z(mu),vx(mm/s),vy(mm/s),vz(mm/s),Epot(nK),Ekin(nK),Etot(nK),index0,index1,distance";
            if      (distance_index == DISTANCE_X    ) out << "(mu)^2";
            else if (distance_index == DISTANCE_V    ) out << "(mm/s)^2";
            else if (distance_index == DISTANCE_PV   ) out << "(mu*mm/s)^2";
            else if (distance_index == DISTANCE_PP   ) out << "(mu*mm/s*amu)^2";
            else if (distance_index == DISTANCE_CROSS) out << "(mu*mm/s*amu)^2";
            else if (distance_index == DISTANCE_MEAN ) out << "(mu*mm/s*amu)^3";
            else if (distance_index == DISTANCE_MAX  ) out << "(mu*mm/s*amu)^2";
            out << std::endl;
        }

        // walk through all cells
        // we can walk through cells_x or equivalently cells_v
        //class cell_data *cell = cells_x.get_first();
        class species_entry *entry;
        class species_data *sp;
        for (; cell != NULL; cell = cell->get_next()) {
            //uint64_t num_atoms = list->n[species];
            num = 0;
            entry = cell->species.get_first();

            //REAL_TYPE **x = list->x[species]; 
            //char   *st = is_atom ? list->status[species] : NULL;

            // atoms
            for (; entry != NULL; entry = entry->get_next(), ++num, ++tot_atoms) {
                //sp = entry->get_species(0);
                sp = &entry->species;
                // atom counter.
                // for atoms this is not ordered but unique. 
                // for molecules we give actual counter - its neither ordered nor unique.
                if (is_atom) {
                    out << sp->atom_index[0] << ",";
                }
                else {
                    out << tot_atoms << ",";
                }

                // position in mu
                out << sp->x[0] << ",";
                out << sp->x[1] << ",";
                out << sp->x[2] << ",";
                
                // momentum in mm/s
                out << sp->x[3] << ",";
                out << sp->x[4] << ",";
                out << sp->x[5] << ",";
                
                // for atoms give status string
                if (is_atom) {
                    if   (IS_STATUS_FREE(sp->status)) {
                        ++_free;
                        out << +sp->status << ",";
                        //out << 0 << ",";
                    }
                    else {
                        ++paired;
                        out << +sp->status << ",";
                        //out << 1 << ",";
                    }
                }
                
                // Epot, Ekin and Etot in nK
          	    out << sp->E[E_POT] << ",";
          	    out << sp->E[E_KIN] << ",";
          	    out << (sp->E[E_POT] + sp->E[E_KIN]);

                if (!is_atom) {
              	    // atom index for molecules and distance^2 or distance^3
              	    out << "," << sp->atom_index[0] << "," << sp->atom_index[1] << "," << sp->distance << std::endl;
          	    }
          	    else {
          	        out << std::endl;
          	    }
            } // next entry
            if (cell->species.get_num() != num) {
                std::cout << std::endl << "error " << species_name << ": unexpected number of atoms in lists = " << num << " != " << cell->species.get_num() << std::endl << std::endl;
                error = -430;
            }
        } // next cell
        out.close();
        
    }

    if (!error) {
        std::cout << "exported " << tot_atoms << " " << species_name << " (" << ((REAL_TYPE)get_ticks_delta(t_start)) << " ms)" << std::endl << std::endl;
    }
    
    return error;
}

/* for each species[0] gets num nearest neighbors of species[1] 
// saves absolute coordinates of found species into file.
// notes: 
// - comparison with molecule search agrees for majoirity but in rare cases not!
//   the reason is that atoms might not be mutually nearest and then pairing is arbitrary.
// - finding num >1 is nearly as fast as num = 1.
// TODO: - search with DISTANCE_PV and DISTANCE_PP works but correctness is not tested. 
//       - search with DISTANCE_X and cell1 = cell0->get_next_distance(cell1, &d) is not fully tested.
//         but it would be much more efficient since most pairs are within same cell.
//       - molecule search could be based on this result easily:
//         execute in parallel by all threads without requirement of locking and save num nearest for each atom.
//         num>1 allows to select alternative atom1 when already paired.
//         however, I do not see a big improvement in efficiency especially for DISTANCE_PV/PP.
int ExportNeighbors(const std::string &filename, 
                    uint8_t index0, 
                    uint8_t index1, 
                    const uint8_t num,
                    const uint8_t distance_index) {
    int error = 0;
    uint32_t tot_atoms = 0;
    uint32_t t_start = GetTickCount();
    std::ofstream out;
    REAL_TYPE *distance         = new REAL_TYPE[num];
    species_data **neighbors = new species_data*[num];
    const std::string name0  = ::cells[index0]->get_species_name(); 
    const std::string name1  = ::cells[index1]->get_species_name();
    REAL_TYPE mass1             = 0.0; //::cells[index1]->get_first()->species.get_first()->get_species(0)->mass;       
    REAL_TYPE *dx               = ::cells[index1]->get_dx();
    REAL_TYPE *dx2              = new REAL_TYPE[::cells[index1]->get_dim()];    
    uint8_t offset           = ::cells[index1]->get_offset();

    // init dx/2
    for (uint8_t i = 0; i < ::cells[index1]->get_dim(); ++i) {
        dx2[i] = dx[i+offset]/2.0;
    }

#ifdef _DEBUG
    {
#else    
    if (distance_index == DISTANCE_PP) {
#endif    
        // get mass of 2nd species: take from first 2nd species atom.
        class cell_data *cell1 = ::cells[index1]->get_first();
        while(cell1) {
            class species_entry *entry1 = cell1->species.get_first();
            if(entry1) {
                mass1 = entry1->get_species(0)->mass;
                break;
            }
            if (mass1 != 0.0) break;
            cell1 = cell1->get_next();
        }
        if (mass1 == 0.0) {
            std::cout << "no species 2 atom (could not get mass)?!" << std::endl;
            error = -1;
        }
    }
    
    if (!error) {
        out.open(filename);
        if (out.fail()) {
            std::cout << "could not open neighbors file '" << filename << "'!" << std::endl;
            error = -2;
        }
        else {
            std::cout << "exporting neighbors file '" << filename << "' ..." << std::endl; 
            
            // create header and init neighbors first time
            out << "index0,x0(mu),y0(mu),z0(mu),vx0(mm/s),vy0(mm/s),vz0(mm/s)";
            for (uint8_t i = 0; i < num; ++i) {
                distance [i] = DBL_MAX;
                neighbors[i] = nullptr;
                out << ",index"<<i+1<<",x"<<i+1<<"(mu),y"<<i+1<<"(mu),z"<<i+1<<"(mu),vx"<<i+1<<"(mm/s),vy"<<i+1<<"(mm/s),vz"<<i+1<<"(mm/s),d"<<i+1;
                if      (distance_index == DISTANCE_X    ) out << "(mu)^2";
                else if (distance_index == DISTANCE_V    ) out << "(mm/s)^2";
                else if (distance_index == DISTANCE_PV   ) out << "(mu*mm/s)^2";
                else if (distance_index == DISTANCE_PP   ) out << "(mu*mm/s*amu)^2";
                else if (distance_index == DISTANCE_CROSS) out << "(mu*mm/s*amu)^2";
                else if (distance_index == DISTANCE_MEAN ) out << "(mu*mm/s*amu)^3";
                else if (distance_index == DISTANCE_MAX  ) out << "(mu*mm/s*amu)^2";
            }
            out << std::endl;

            // walk through all cells and entries of species 0 
            class cell_data *cell1;
            class species_entry *entry0, *entry1;
            class species_data *atom0, *atom1;
            class cell_data *cell0 = ::cells[index0]->get_first();
            for (; (!error) && (cell0 != NULL); cell0 = cell0->get_next()) {
                entry0 = cell0->species.get_first();

                // atom 0
                for (; (!error) && (entry0 != NULL); entry0 = entry0->get_next(), ++tot_atoms) {
                    atom0 = entry0->get_species(0);
                    // atom counter
                    //out << tot_atoms << ",";
                    out << atom0->atom_index[0] << ",";

                    // position in mu
                    out << atom0->x[0] << ",";
                    out << atom0->x[1] << ",";
                    out << atom0->x[2] << ",";
                    
                    // momentum in mm/s
                    out << atom0->x[3] << ",";
                    out << atom0->x[4] << ",";
                    out << atom0->x[5];

                    uint32_t n = 0; // loop counter. this is not atom index

                    // walk through all cells and entries of species 1
                    cell1 = ::cells[index1]->get_first();
                    while (cell1 != NULL) {
                        // atom 1
                        // searches atom1 with distance <= largest distance of num nearest neighbors
                        // for the first num nearest neighbors returns consecutive atoms since distance[num-1] == 0.0.
                        // for further searches returns only atoms with smaller distance.
                        REAL_TYPE d = DBL_MAX; // DBL_MAX = infinite
                        entry1 = entry0->get_next_distance(cell1->species.get_first(), &d, distance_index); 
#ifdef _DEBUG
                        if (cell1->species.get_num() == 0) {
                            // empty cell: expect NULL and d == DBL_MAX
                            if ((entry1 != nullptr) || (d != DBL_MAX)) {
                                error = -340;
                                break;
                            }
                        }
                        else {
                            // non-empty cell: expect first atom of cell and finite distance
                            if ((entry1 != cell1->species.get_first()) || (d == DBL_MAX)) {
                                error = -342;
                                break;
                            }
                        }
#endif                        
                        for (; entry1 != NULL; ++n) {
                            if (entry1 != entry0) {
                                atom1 = entry1->get_species(0);
                                //std::cout << "atom " << n << ": " << d << " ... " << std::endl;
                                for (uint8_t i = num-1; ; --i) {
                                    if (distance[i] != DBL_MAX) {
                                        if (d < distance[i]) {
                                            if ((i+1) < num) {
                                                distance [i+1] = distance [i];
                                                neighbors[i+1] = neighbors[i];
                                            }
                                            if (i == 0) {
                                                distance [i] = d;
                                                neighbors[i] = atom1;
                                                break;
                                            }
                                        }
                                        else {
                                            if ((i+1) < num) {
                                                distance [i+1] = d;
                                                neighbors[i+1] = atom1;
                                            }
                                            break;
                                        }
                                    }
                                    else if (i == 0) {
                                        distance [i] = d;
                                        neighbors[i] = atom1;
                                        break;
                                    }
                                    if (i == 0) break;  
                                }
                            }
                            // find next atom with d <= distance[num-1]
                            d = distance[num-1];
                            entry1 = entry0->get_next_distance(entry1->get_next(), &d, distance_index);
                        } // next entry1
                        
                        // find next cell with d <= distance[num-1]
                        d = distance[num-1];
                        cell1 = entry0->get_next_distance(cell1->get_next(), &d, distance_index, dx2, mass1);

                    } // next cell1
                    
                    //std::cout << "atom " << tot_atoms << " found nearest " << num << " / " << n << std::endl;
                    
                    // save absolute coordinates and distance to species[0] of num nearest neighbors of species[1]
                    for (uint8_t i = 0; i < num; ++i) {
                        if (neighbors[i] != nullptr) {
                            atom1 = neighbors[i];
                            // index
                            out << ",";
                            out << atom1->atom_index[0] << ",";
                            // absolute position in mu
                            out << atom1->x[0] << ",";
                            out << atom1->x[1] << ",";
                            out << atom1->x[2] << ",";
                            
                            // absolute velocity in mm/s
                            out << atom1->x[3] << ",";
                            out << atom1->x[4] << ",";
                            out << atom1->x[5] << ",";

                            // save distance^2 or distance^3
                            //out << std::sqrt(distance[i]);
                            out << distance[i];
                     
                            // reset for next atom 0                               
                            distance [i] = DBL_MAX;
                            neighbors[i] = nullptr;
                        }
                        else break;
                    }
                    out << std::endl;
                }// next entry0 (atom 0)
            }// next cell0

            out.close();
        }
    
        std::cout << "exported " << tot_atoms << " " << name0 << " with " << +num << " nearest " << mol_distance_use << " neighbors of " << name1 << " (" << (((REAL_TYPE)get_ticks_delta(t_start))/1e3) << " s)" << std::endl << std::endl;
        
    }
    
    delete[] distance;
    delete[] neighbors;
    delete[] dx2;

    return error;
}*/


// convert string to number of given data type
// returns 0 if ok or -1 if failed
// notes: 
// - other than atof() and std::stod() ensures all characters in string are converted, otherwise returns -1.
//   tellg() could be also applied and would return number of converted chars,
//   but its implemenation is platform dependent. eof() is more reliable.
//   tellg() works on Ubuntu with g++ but caused that fail() returns true after this is called.
// - this works also with strings as type.
// - if allow_hex = true checks if string after optional spaces or tabs starts with NUM_HEX
//                       and reads string as hexadecimal number 
template<typename T>
int to_number(const std::string& str, T *ptr, const std::string& info, bool allow_hex=false, bool output=true) {
    bool is_hex = false;
    std::stringstream stream(str);
    T d;
    if (allow_hex) {
        const auto start = str.find_first_not_of(" \t");
        if ((str.size() - start) >= 3) {
            is_hex = ((str.compare(start, 2, NUM_HEX) == 0) && (std::isxdigit(str[start+2])));
            if (is_hex) stream.ignore(start + 2);
        }
    }
    if (is_hex) stream >> std::hex >> d;
    else        stream >> d;
    if (stream.fail()) {
        std::cout << "cannot convert string '" << str << "'! (failed): " << info << std::endl;
        //std::runtime_error e(numberAsString);
        //throw e;
        return -1;
    }
    else if (stream.bad()) {
        std::cout << "cannot convert string '" << str << "'! (bad): " << info << std::endl;
        //std::runtime_error e(numberAsString);
        //throw e;
        return -1;
    }
    else if (!stream.eof()) {
        std::cout << "cannot convert string '" << str << "'! (partial): " << info << std::endl;
        //std::runtime_error e(numberAsString);
        //throw e;
        return -1;
    }
    else {
        //*static_cast<T*>(ptr) = d;
        *ptr = d;
        if (output) {
            if (is_hex) std::cout << "0x" << std::hex << d << std::dec;
            else        std::cout << d;
        }
    }
    return 0;
}

// read a vector of values between V_OPEN and V_CLOSE separated by V_SEP
// if allow_hex = true hexadecimal values starting with NUM_HEX are allowed.
template<typename T>
int to_vector(const std::string& str, std::vector<T> &vec, const std::string& info, bool allow_hex) {
    int error = 0;
    auto start = str.find_first_not_of(" \t");
    if (str[start] != V_OPEN) error = -1;
    else {
        start += 1;
        const auto end = str.find(V_CLOSE, start);
        if ((end == std::string::npos) || (end == start)) error = -2;
        else {
            std::stringstream stream(str.substr(start, end-start));
            std::string tmp;
            T value;
            //std::cout << "reading vector in '" << str.substr(start, end-start) << "' ( " << start << ", " << end-start << " )" << std::endl;
            while(std::getline(stream, tmp, V_SEP)) {
                error = to_number<T>(tmp, &value, info, allow_hex, false);
                //std::cout << ": to_number ('" << tmp << "')= " << value << " (error " << error << ")" << std::endl;
                if (error) break;
                vec.push_back(value);
            }
        }
    }
    return error;
}


void show_cwd(char* argv[]) {
    std::cout << "executable                '" << argv[0] << "'" << std::endl;
#if defined(_WIN32) || defined(_WIN64)
    TCHAR _cwd[MAX_PATH + 1] = L"";
    DWORD len = GetCurrentDirectory(MAX_PATH, _cwd);
    std::wstring cwd(_cwd);
    std::wcout << L"current working directory '" << cwd << L"'" << std::endl;
#else
    char *cwd = new char[101];
    if (getcwd(cwd, 100)) {
        std::cout << "current working directory '" << cwd << "'" << std::endl;
    }
    delete[] cwd;
#endif

}

int ReadParams(const char *filename, char* argv[]) {
    // read parameters from given file + path
    // TODO: it would be convenient to give also vectors as list of numbers
    int error = 0, line = 1, count = 0;
    std::ifstream in;
    in.open(filename);
    if (!in.is_open()) {
        std::cout << "error open file '" << filename << "'!" << std::endl;
        show_cwd(argv);
        error = -500;
    }
    else {
        std::string buffer;
        bool found;
        size_t len1, len2, i, c;
        while(getline(in, buffer)) {
            found = false;
            len1 = buffer.length();
            if ((len1 > 0) && (buffer[0] != '*') && (buffer[0] != '/') && (buffer[0] != '#')) {
                for (i=0; (i<NUM_VAR) && (!found) && (!error); ++i) {
                    len2 = var_params[i].name.length();
                    if (len1 >= len2) {
                        c = 0;
                        while (true) {
                            if (buffer[c] != var_params[i].name[c]) break;
                            if (++c == len2) {
                                // parameter found
                                // read value after '=' skipping leading spaces.
                                // number must be terminated with new line and nothing else!
                                ++c;
                                while( (c < len1) && ( (buffer[c] == ' ') || (buffer[c] == '=') ) ) ++c; 
                                if (c >= len1) {
                                    std::cout << "line " << line << ": '" << buffer << "' unexpected!" << std::endl;
                                    error = -501;
                                    break;
                                }
                                else if (var_params[i].flag & FLAG_SET) {
                                    std::cout << "line " << line << ": '" << buffer << "' already set!" << std::endl;
                                    error = -502;
                                    break;
                                }
                                else {
                                    std::cout << "line " << line << ": " << var_params[i].name << " = ";
                                    switch(var_params[i].flag & TYPE_MASK) {
                                        case TYPE_U8: 
                                            { 
                                                // uint8_t fails with stream.eof() but uint16_t works.
                                                uint16_t tmp = 0;
                                                error = to_number<uint16_t>(buffer.substr(c), 
                                                    &tmp, "8-bit unsigned integer required"); 
                                                if (!error) {
                                                    if (tmp < 256) *((uint8_t*)var_params[i].ptr) = (uint8_t)tmp;
                                                    else error = -511;
                                                }
                                            }
                                            break;
                                        case TYPE_U16: 
                                            error = to_number<uint16_t>(buffer.substr(c), 
                                                (uint16_t*)var_params[i].ptr, "16-bit unsigned integer required"); 
                                            break;
                                        case TYPE_U32: 
                                            error = to_number<uint32_t>(buffer.substr(c), 
                                                (uint32_t*)var_params[i].ptr, "32-bit unsigned integer required"); 
                                            break;
                                        case TYPE_U64: 
                                            error = to_number<uint64_t>(buffer.substr(c), 
                                                (uint64_t*)var_params[i].ptr, "64-bit unsigned integer required"); 
                                            break;
                                        case TYPE_DOUBLE: 
                                            error = to_number<REAL_TYPE>(buffer.substr(c), 
                                                (REAL_TYPE*)var_params[i].ptr, "real number required"); 
                                            break;
                                        case TYPE_STRING: 
                                            error = to_number<std::string>(buffer.substr(c), 
                                                (std::string*)var_params[i].ptr, "string required (no space allowed!)"); 
                                            break;
                                        default: 
                                            error = -510; 
                                            break;
                                    }
                                    if (error) {
                                       std::cout << "'" << buffer << "' failed with error " << error << std::endl;
                                    }
                                    else {                                    
                                        found = true;
                                        var_params[i].flag |= FLAG_SET;
                                        ++count;
                                        std::cout << " (ok)" << std::endl;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                if ( (!error) && (!found) ) {
                    std::cout << "line " << line << ": '" << buffer << "' unknown parameter!" << std::endl;
                    error = -520;
                }
            }
            if (error) break;
            ++line;
        }
        in.close();
        if (error) std::cout << "'" << filename << "' error " << error << " in line " << line << "!" << std::endl << std::endl;
        else       std::cout << "'" << filename << "' " << line << " lines, " << count << " parameters ok." << std::endl << std::endl;
    }
    return error;
}

// save histogram to given file for given species_name.
// if atom_max != NULL and for parameters with atom_max[i] == 0.0 the maximum values is obtained from the data. 
// length of atom_max = DIM_MAXVAL since includes Emax. original atom_max are not changed.
// atom/molecule data is obtained from list->x[species][coordinate][data].
// species 0 = atom 0, 1 = atom 1, 2 = molecule.
// the histogram for all list->num threads is taken by walking through list->next.
// returns 0 if ok, otherwise error 
// note: excute only on main thread with list from any thread.
int Histogram(std::string &filename, std::string &species_name, std::string &stat,
                uint16_t rep, uint16_t num_bins, //uint8_t species,
                //struct thread_data *list,
                //class cell_data *cells,
                uint8_t index,
                REAL_TYPE T, REAL_TYPE mu,
                REAL_TYPE wx, REAL_TYPE wy, REAL_TYPE wz,
                REAL_TYPE Mass,
                bool use_atom_max) {
    int error = 0;
    class cell_data *cells = ::cells[index]->get_first();
    uint64_t Ntot = 0;
    //uint64_t num_lists = list->num; // number of lists = number of threads
    REAL_TYPE h;
    uint64_t **H = new uint64_t*[DIM_MAXVAL]; // histogram for DIM_PHASESPACE coordinates + Etot
    REAL_TYPE MaxV[DIM_MAXVAL]; // maximum value used for histogram
    bool findMax[DIM_MAXVAL]; // true if have to find maximum value for each dimension
    bool findAny = false; // true if have to find any maximum value
#ifdef HIST_EPOT_EKIN
    uint64_t *HEpot = new uint64_t[num_bins];
    uint64_t *HEkin = new uint64_t[num_bins];
    for (uint16_t k = 0; k < num_bins; ++k) {
        HEpot[k] = HEkin[k] = 0;
    }
#endif    
    
    // allocate histogram and check if maximum values need to be found
    for (uint8_t j = 0; j < DIM_MAXVAL; ++j) {
        H[j] = new uint64_t[num_bins];
        for (uint16_t k = 0; k < num_bins; ++k) {
            H[j][k] = 0;
        }
        // get max values for any zero max value
        MaxV[j] = (use_atom_max) ? ::cells[index]->get_max(j) : 0.0;
        if (MaxV[j] == 0.0) {
            findMax[j] = true;
            findAny = true;
        }
    }
        
    if ( findAny ) {
        // have to find at least one maximum value
        
        // walk through all cells
        // we can walk through cells_x or equivalently cells_v
        //class cell_data *cell = cells_x.get_first();
        class cell_data *cell = cells;
        class species_entry *entry;
        class species_data *sp;
        for (; cell != NULL; cell = cell->get_next()) {
            //uint64_t num_atoms = list->n[species];
            // walk through all atoms
            entry = cell->species.get_first();
            for (; entry != NULL; entry = entry->get_next()) {
                //sp = entry->get_species(0);
                sp = &entry->species;
            
                //for (uint64_t i = 0; i < num_lists; ++i, list = list->next) {
                //uint64_t num_atoms = list->n[species];
                // dimension
                for (uint8_t j = 0; j < DIM_PHASESPACE; ++j) {
                    if (findMax[j]) {
                        h = abs(sp->x[j]);
                        if (h > MaxV[j]) MaxV[j] = h;
                    }
                }
                if (findMax[DIM_PHASESPACE]) {
                    /* calculate Etot
                    for (uint64_t k = 0; k < num_atoms; ++k, ++xk) {
                        REAL_TYPE hx = list->x[species][0][k]*wx;
                        REAL_TYPE hy = list->x[species][1][k]*wy;
                        REAL_TYPE hz = list->x[species][2][k]*wz;        
                        REAL_TYPE Epot = 0.5*Mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
                        hx = list->x[species][3][k];
                        hy = list->x[species][4][k];
                        hz = list->x[species][5][k];
                        REAL_TYPE Ekin = 0.5*Mass*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
                        h = abs(Epot + Ekin);
                        if (h > MaxV[DIM_PHASESPACE]) MaxV[DIM_PHASESPACE] = h;
                    }
                    */
                    h = sp->E[E_POT] + sp->E[E_KIN];
                    if (h > MaxV[DIM_PHASESPACE]) MaxV[DIM_PHASESPACE] = h;
                }
            }
        }
    }
    
    
    // walk through all cells
    // we can walk through cells_x or equivalently cells_v
    //class cell_data *cell = cells_x.get_first();
    class cell_data *cell = cells;
    class species_entry *entry;
    class species_data *sp;
    for (; cell != NULL; cell = cell->get_next()) {
        //uint64_t num_atoms = list->n[species];
        // walk through all atoms
        entry = cell->species.get_first();
        for (; entry != NULL; entry = entry->get_next()) {
            //sp = entry->get_species(0);
            sp = &entry->species;
            //uint64_t num_atoms = list->n[species];
            //Ntot += num_atoms;
            ++Ntot;
            // dimension
            for (uint8_t j = 0; j < DIM_MAXVAL; ++j) {
                //REAL_TYPE *xk = list->x[species][j];
                uint16_t bin;
                if (j < DIM_PHASESPACE) {
                    // position and momentum REAL_TYPE-sided
                    bin = (uint16_t)std::floor(num_bins*0.5*((sp->x[j]/MaxV[j])+1));
                }
                else {
                    /* calculate single-sided Etot
                    REAL_TYPE hx = list->x[species][0][k]*wx;
                    REAL_TYPE hy = list->x[species][1][k]*wy;
                    REAL_TYPE hz = list->x[species][2][k]*wz;        
                    REAL_TYPE Epot = 0.5*Mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK */
#ifdef HIST_EPOT_EKIN
                    // get histogram of Epot
                    bin = (uint16_t)std::floor(num_bins*((sp->E[E_POT]/MaxV[j])));
                    if (bin >= num_bins) bin = num_bins-1;
                    else if (bin < 0)    bin = 0;
                    ++HEpot[bin];
#endif                        
                    /*hx = list->x[species][3][k];
                    hy = list->x[species][4][k];
                    hz = list->x[species][5][k];
                    REAL_TYPE Ekin = 0.5*Mass*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK*/
#ifdef HIST_EPOT_EKIN
                    // get histogram of Ekin
                    bin = (uint16_t)std::floor(num_bins*((sp->E[E_KIN]/MaxV[j])));
                    if (bin >= num_bins) bin = num_bins-1;
                    else if (bin < 0)    bin = 0;
                    ++HEkin[bin];
#endif                        
                    // get histogram of Etot
                    bin = (uint16_t)std::floor(num_bins*(((sp->E[E_POT] + sp->E[E_KIN])/MaxV[j])));
                }
                if (bin >= num_bins) bin = num_bins-1;
                else if (bin < 0)    bin = 0;
                ++H[j][bin];
            }
        }        
    }

    if (true) {
        // Andi: ensure sum of bins == total N
        uint64_t sum;
        // dimension
        for(uint8_t j = 0; j < DIM_MAXVAL; ++ j) { 
            // bins
            sum = 0;
            for(uint16_t k = 0; k < num_bins; ++k) sum += H[j][k];
            if (sum != Ntot) { 
                std::cout << std::endl << "error sum of " << cname[j] << " bins " << sum << " != " << Ntot << std::endl << std::endl;
                error = -600;
            }
        }
    }
    
    std::ofstream out;
    out.open(filename);
    if (out.fail()) {
        std::cout << "could not open Histogram file '" << filename << "'!" << std::endl;
        error = -4;
    }
    else {
        out << "// "          << species_name << " histogram" << std::endl;
        out << "// stat   = " << stat           << std::endl;
        out << "// rep    = " << rep            << std::endl;
        out << "// N      = " << Ntot           << std::endl;
        out << "// T      = " << T     << " nK" << std::endl;
        out << "// mu     = " << mu    << " nK" << std::endl;
        for (uint8_t j = 0; j < DIM_MAXVAL; ++j) {
            out << "// max " << cname[j] << " = " << MaxV[j] << " " << cunit[j] << std::endl;
        }
        std::string header = "// #";
        for (uint8_t j = 0; j < DIM_MAXVAL; ++j) {
            header += " " + cname[j];
        }
#ifdef HIST_EPOT_EKIN
        header += " Epot Ekin";
#endif
        out << header << std::endl;
        for (uint8_t j = 0; j < header.length(); ++j) {
            out << "/";
        }
        out << std::endl;
        for (uint16_t k = 0; k < num_bins; ++k) {
            out << k; 
            for (uint8_t j = 0; j < DIM_MAXVAL; ++j) {
                out << " " << H[j][k];
            }
#ifdef HIST_EPOT_EKIN
            out << " " << HEpot[k] << " " << HEkin[k] << std::endl;
#else
            out << std::endl;
#endif
        }
        out << std::endl;
        out.close();
    }
        
#ifdef HIST_EPOT_EKIN
    delete [] HEpot;
    delete [] HEkin;
#endif    

    for (uint8_t j = 0; j < DIM_MAXVAL; ++j) delete [] H[j];
    delete [] H;

    return error;
}

// get calculation volume
// TODO: this could be done more dynamic and could be part of cell generation: 
//       select Emax from distribution p(E) where N*p(Emax) << 1
//       this gives a radius for arb. potential U(rmax) = Emax with k~1
//       the cells could be adjusted to be within this radius, i.e. instead of a cube would be an ellipsoid.
//       the cell size could be also adjusted: 
//       instead of constant cartesian volume one could do constant dR*dPhi*dTheta in spherical or cylindrical coords.
//       i.e. cells farther out are larger but might still have fewer atoms
//       maybe one can also split cells when they have more than a threshold of atoms (say 1e4)
//       and combine cells when they have fewer than a threshold atoms (say 1e2)
void get_size(uint8_t species, uint8_t statistics, REAL_TYPE **max, REAL_TYPE **R) {
    // allocate vectors
    *max = new REAL_TYPE[DIM_MAXVAL];
    *R   = new REAL_TYPE[DIM_SPACE];
    std::cout << "max_E = " << ::max_E[species] << " nK, Tc = " << atom_Tcrit[species] << " nK, T = " << atom_T[species] << " nK, size = " << CalculationSize[species] << " * " << EMAX_SCALING << std::endl;
    (*max)[DIM_MAXVAL - 1] = ::max_E[species]; // copy Emax into vector

    if (statistics == STAT_BE_C) {
        // recalculate the max_ values for BEC, otherwise it is extremely slow
        //exp(-mass/k_nK*((x[0]*x[0]*omega[0] + x[1]*x[1]*omega[1] + x[2]*x[2]*omega[2])*k_Epot + 
        //                (x[3]*x[3]/omega[0] + x[4]*x[4]/omega[1] + x[5]*x[5]/omega[2])*k_Ekin ))
        REAL_TYPE fmax = 0.7; // fmax^6 = probability of seeing 1 out of N atoms at max_..: 0.7^6 = 0.1, 0.5^6 = 0.015.
        (*max)[0] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Epot)/::atom_wx[species]);
        (*max)[1] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Epot)/::atom_wy[species]);
        (*max)[2] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Epot)/::atom_wz[species]);

        (*max)[3] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Ekin)*::atom_wx[species]);
        (*max)[4] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Ekin)*::atom_wy[species]);
        (*max)[5] = CalculationSize[species]*EMAX_SCALING*std::sqrt(std::log(::atom_N[species]/fmax)*k_nK/(::atom_mass[species]*k_Ekin)*::atom_wz[species]);
        
        // Radius of condensed part = ground state of non-interacting HO in um = Gaussian.
        // see Ketterle, Making-Probing-Understanding, Equ. 38. [k_Rc = 1e6*np.sqrt(hbar/mn)] 
        // note: factor 0.5 is not in Equ. 38 but accounts for definition of Gauss sigma.
        // TODO: in Thomas-Fermi limit - strong interactions gives inverted parabola with R = sqrt(2*mu/m*omega^2)
        (*R)[0] = std::sqrt(0.5/(atom_mass[species]*atom_wx[species]))*k_Rc; // mu
        (*R)[1] = std::sqrt(0.5/(atom_mass[species]*atom_wy[species]))*k_Rc;
        (*R)[2] = std::sqrt(0.5/(atom_mass[species]*atom_wz[species]))*k_Rc;
    }
    else {
        // Epot = 1/2*m*omega^2*x^2 * k_Epot
        (*max)[0] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*atom_wx[species]*atom_wx[species]*k_Epot)); // mu
        (*max)[1] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*atom_wy[species]*atom_wy[species]*k_Epot));
        (*max)[2] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*atom_wz[species]*atom_wz[species]*k_Epot));

        // E_kin = 1/2*m*v^2 * k_Ekin 
        (*max)[3] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*k_Ekin)); // mm/s
        (*max)[4] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*k_Ekin));
        (*max)[5] = std::sqrt(2.0*max_E[species]/(atom_mass[species]*k_Ekin));

        if (statistics == STAT_FD) {
            // radius used in Polylog
            (*R)[0] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wx[species]*atom_wx[species]*k_Epot)); // mu
            (*R)[1] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wy[species]*atom_wy[species]*k_Epot));
            (*R)[2] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wz[species]*atom_wz[species]*k_Epot));
        }
        else if ((statistics == STAT_BE   ) || 
                 (statistics == STAT_BE_NC)) {
            // radius of non-condensed part used in PolyLog
            (*R)[0] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wx[species]*atom_wx[species]*k_Epot)); // mu
            (*R)[1] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wy[species]*atom_wy[species]*k_Epot));
            (*R)[2] = std::sqrt(2.0*atom_T[species]/(atom_mass[species]*atom_wz[species]*atom_wz[species]*k_Epot));
        }
        else if (statistics == STAT_MB) {
            // Gauss sigma
            (*R)[0] = std::sqrt(atom_T[species]/(atom_mass[species]*atom_wx[species]*atom_wx[species]*k_Epot)); // mu
            (*R)[1] = std::sqrt(atom_T[species]/(atom_mass[species]*atom_wy[species]*atom_wy[species]*k_Epot));
            (*R)[2] = std::sqrt(atom_T[species]/(atom_mass[species]*atom_wz[species]*atom_wz[species]*k_Epot));
        }
        else {
            // unknown distribution (should not happen)
            //return std::numeric_limits<REAL_TYPE>::quiet_NaN();
            throw ERROR_UNKNOWN_DIST;
        }

    }
}

// calculate chemical potential in nK for non-interacting Fermions or Bososons or thermal gas.
// integrate FD or BE distribution over phase space gives atom number
// in a loop adjust chemical potential to get AtomNumber == N[species]
// throws ERROR_UNKNOWN_DIST if species distribution is unknown
// note: this assumes harmonic trap
// TODO: I have a newer code on python which works pretty well and adapts maximum energy and energy steps automatically.
//       however, present code here is most likely faster but uses fixed size which might be a limit.
REAL_TYPE get_mu(uint8_t species) {
    uint32_t t_start = GetTickCount();
    REAL_TYPE AtomNumber = 0.0;
    REAL_TYPE who = std::pow(atom_wx[species]*atom_wy[species]*atom_wz[species], 1.0/3.0);
    //REAL_TYPE PreFactor=1.0/(2*std::pow(hbar*who,3.0)); // this is 5e91!!!
    REAL_TYPE PreFactor = who*k_nK; // k_nK = (hbar*1e9)/kB; omega -> nK
    //REAL_TYPE beta=1.0/(kB*T[species]);
    //REAL_TYPE FermiEnergy=pow(6*N[species],1.0/3.0)*hbar*who;
    //REAL_TYPE ChemicalPotential=FermiEnergy;
    //REAL_TYPE StepSize=InitialStepSize*FermiEnergy;
    REAL_TYPE beta = 1.0/atom_T[species];
    REAL_TYPE ChemicalPotential;
    REAL_TYPE plus_minus_one;
    REAL_TYPE limit;
    bool   check_N = true, limit_reached = false;

    if (atom_stat_index[species] == STAT_FD) {
        // Fermi-Dirac statistics: Tcrit = TFermi
        ::atom_Tcrit[species] = std::pow(6 * atom_N[species], 1.0 / 3.0) * who * k_nK; // k_nK = (hbar*1e9)/kB; omega -> nK
        // cutoff energy in nK. TODO: check if reasonable [we already check E > Emax]
        ::max_E[species] = CalculationSize[species] * EMAX_SCALING * (::atom_Tcrit[species] + 4 * ::atom_T[species]);
        // starting energy = TF and limit mu <= TF
        ChemicalPotential = atom_Tcrit[species];
        // limit mu <= TF. we give 5% margin to take into account numerical errors. check_N = true even when this is reached.
        limit = ::atom_Tcrit[species]*1.05;
        // statistics contains +1
        plus_minus_one = +1.0;
    }
    else if ((atom_stat_index[species] == STAT_BE   ) || 
             (atom_stat_index[species] == STAT_BE_NC) || 
             (atom_stat_index[species] == STAT_BE_C )) {
        // Bose-Einstein statistics: Tcrit = Tc
        // see: Ketterle MakingProbingUnderstanding, Equ. 51 (note: Equ. on top of page 40 is wrong)
        //      Zwierlein MakingProbingUnderstanding, Equ. 29
        //      zeta_3=Li_3(1) must be is inside sqrt.
        ::atom_Tcrit[species] = std::pow(::atom_N[species] / zeta_3, 1.0 / 3.0) * who * k_nK; // k_nK = (hbar*1e9)/kB; omega -> nK.
        REAL_TYPE TTc = ::atom_T[species] / ::atom_Tcrit[species];
        REAL_TYPE fc = 1.0 - TTc * TTc * TTc;
        if (fc < 1e-10) fc = 0.0;
        ::atom_Nc[species] = uint64_t(std::round(fc * ::atom_N[species]));
        // cutoff energy in nK. TODO: check if reasonable [we already check E > Emax]
        ::max_E[species] = CalculationSize[species] * EMAX_SCALING * (::atom_Tcrit[species] + 4 * ::atom_T[species]);
        // starting energy = -T
        ChemicalPotential = -::atom_T[species];
        // limit mu <= 0.0. this limit is sharp to avoid singularity. check_N = false when this is reached.
        limit = 0.0;
        // statistics contains -1
        plus_minus_one = -1.0;
    }
    else if (atom_stat_index[species] == STAT_MB) {
        // for Maxwell-Boltzmann distribution Tcrit = 0
        ::atom_Tcrit[species] = 0.0;
        // cutoff energy in nK. TODO: check if reasonable [we already check E > Emax]
        ::max_E[species] = CalculationSize[species] * EMAX_SCALING * 4.0 * ::atom_T[species];
        return 0.0;
    }
    else {
        // unknown distribution (should not happen)
        //return std::numeric_limits<REAL_TYPE>::quiet_NaN();
        throw ERROR_UNKNOWN_DIST;
    }

    REAL_TYPE StepSize = InitialStepSize*fabs(ChemicalPotential);
    bool Direction = true;
    //REAL_TYPE de = atom_max[species][DIM_MAXVAL-1]/ESteps[species];
    REAL_TYPE de = max_E[species]/ESteps[species]; // TODO: max_E and delta_E (de) could be determined here in a loop - see my Python code.
    uint64_t loops = 0;
    do {
        AtomNumber = 0;
        uint64_t e = (limit_reached && (!check_N)) ? 1 : 0; // for mu = 0 and Energy = 0 avoid division by 0 
        for (; e < ESteps[species]; e++) {
            REAL_TYPE Energy = e * de;
            AtomNumber += de*Energy*Energy/(exp(beta*(Energy-ChemicalPotential))+plus_minus_one);    
        }
        //AtomNumber=AtomNumber*PreFactor;
        AtomNumber /= 2*PreFactor*PreFactor*PreFactor;
        if (AtomNumber > atom_N[species]) {
            if (!Direction) {
                StepSize *= 0.8;
                std::cout << "+";
                Direction = true;
                limit_reached = false;
            } 
            ChemicalPotential -= StepSize;
        } else {
            if (Direction) {
                StepSize *= 0.8;
                std::cout << "-";
            } 
            ChemicalPotential += StepSize;
            if (ChemicalPotential > limit) {
                // limit mu but do not abort, just reduce StepSize when not already done
                // for Bosons with T <= Tc: AtomNumber = non-condensed fraction <= N[species]. atom number cannot be used as abort criterion.
                // for Fermions we limit to mu <= TF*1.05, atom number should be still fulfilled within error. 
                ChemicalPotential = limit;
                if (limit_reached) {
                    StepSize *= 0.8;
                    std::cout << "^";
                }
                else {
                    limit_reached = true; 
                    if ((atom_stat_index[species] == STAT_BE   ) || 
                        (atom_stat_index[species] == STAT_BE_NC) || 
                        (atom_stat_index[species] == STAT_BE_C )) check_N = false;
                    SET_CONSOLE_COLOR_NOTE(STDOUT_HANDLE);
                    std::cout << std::endl << "note: chemical potential " << ::atom_name[species] << " limit mu to " << limit << "!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    if (!Direction) {
                        StepSize *= 0.8;
                        std::cout << "^";
                    }
                }
            }
            Direction = false;
        }

        if (++loops >= 1000) {
            if ((atom_stat_index[species] == STAT_FD) && (ChemicalPotential == limit)) {
                REAL_TYPE TTF = atom_T[species] / atom_Tcrit[species];
                if (TTF <= 0.025) {
                    if (loops == 1000) {
                        // at low T/TF the loop is stuck at limit of mu without improving error. we can savely assume mu = TF.
                        // continue one loop to recalculate N error.
                        SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                        std::cout << std::endl << std::endl << "warning: maximum number of loops " << loops << " reached with mu limited!" << std::endl << "This is at low T/TF = " << TTF << " and we can assume mu = TFermi = " << atom_Tcrit[species] << " nK." << std::endl;
                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                        ChemicalPotential = atom_Tcrit[species];
                        continue;
                    }
                }
            }
            break;
        }

    }
    while ( (check_N && (fabs((AtomNumber-atom_N[species])/atom_N[species]) > ErrorAllowed)) || (StepSize > error_mu) );

    std::cout<<std::endl;
    //std::cout << "PreFactor = " << PreFactor << std::endl;
    std::cout << "ChemicalPotential " << ::atom_name[species] << " (" << atom_stat_use[species] << ") = " << ChemicalPotential << " nK +/- " 
        << std::scientific << std::setprecision(OUTPUT_PRECISION) << StepSize << " nK, N error = " 
        << fabs(AtomNumber-atom_N[species])
        << std::fixed << std::setprecision(OUTPUT_PRECISION) << " (" 
        << loops << " loops, " << get_ticks_delta(t_start) << " ms)" << std::endl << std::endl;

    if (loops > 1000) {
    }

    return ChemicalPotential;
}

// probability function of thermal Maxwell-Boltzmann gas in a harmonic trap with given parameters
// this is in units of (hbar*2*pi)^3 [e.g. Equ. 20 in Zwierlein, Ketterle, MPU Fermi Gases]
// saves Epot and Ekin into given pointers
REAL_TYPE f_MB(REAL_TYPE mass, REAL_TYPE T, REAL_TYPE mu, REAL_TYPE omega[DIM_SPACE], REAL_TYPE vx[DIM_PHASESPACE], REAL_TYPE *Epot, REAL_TYPE *Ekin, REAL_TYPE Emin) {
    REAL_TYPE hx = omega[0]*vx[0]; // rad/s*mu
    REAL_TYPE hy = omega[1]*vx[1];
    REAL_TYPE hz = omega[2]*vx[2];
    *Epot = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
    hx = vx[3]; // mm/s
    hy = vx[4];
    hz = vx[5];
    *Ekin = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
    return exp(-((*Epot) + (*Ekin))/T);
}

// probability function of Fermions in a harmonic trap
// this is in units of (hbar*2*pi)^3
// saves Epot and Ekin into given pointers
REAL_TYPE f_FD(REAL_TYPE mass, REAL_TYPE T, REAL_TYPE mu, REAL_TYPE omega[DIM_SPACE], REAL_TYPE x[DIM_PHASESPACE], REAL_TYPE *Epot, REAL_TYPE *Ekin, REAL_TYPE Emin) {
    REAL_TYPE hx = omega[0]*x[0]; // rad/s*mu
    REAL_TYPE hy = omega[1]*x[1];
    REAL_TYPE hz = omega[2]*x[2];
    *Epot = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
    hx = x[3]; // mm/s
    hy = x[4];
    hz = x[5];
    *Ekin = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
    return (1.0/(exp(((*Epot) + (*Ekin) - mu)/T) + 1.0));
}

// probability function of non-condensed Bosons in a harmonic trap
// this is in units of (hbar*2*pi)^3
// saves Epot and Ekin into given pointers
// for T <= Tc: mu = 0 and we have to limit Emin to ground-state energy, otherwise function diverges.
// in this case f_BE returns 0.0 and no particle should be created. 
// TODO: maybe all the limited atoms = BEC fraction? -> no, this is much smaller!
REAL_TYPE f_BE(REAL_TYPE mass, REAL_TYPE T, REAL_TYPE mu, REAL_TYPE omega[DIM_SPACE], REAL_TYPE x[DIM_PHASESPACE], REAL_TYPE *Epot, REAL_TYPE *Ekin, REAL_TYPE Emin) {
    REAL_TYPE hx = omega[0]*x[0]; // rad/s*mu
    REAL_TYPE hy = omega[1]*x[1];
    REAL_TYPE hz = omega[2]*x[2];
    *Epot = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
    hx = x[3]; // mm/s
    hy = x[4];
    hz = x[5];
    *Ekin = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
    REAL_TYPE Etot = (*Epot) + (*Ekin) - mu;
    //REAL_TYPE Emin = 0.5*(omega[0] + omega[1] + omega[2])*k_nK; // ground-state energy in nK
    //REAL_TYPE Emin = 0.5*std::pow(omega[0]*omega[1]*omega[2],1.0/3.0)*k_nK; // ground-state energy in nK
    if (Etot < Emin) return 0.0;
    else             return (1.0/(exp(Etot/T) - 1.0));
}

// probability function of condensed Bosons in a harmonic trap
// chemical potential mu=0 for BEC. 
// BEC fraction = max(1-(T/Tc)^3,0) = 0..1
// this is in units of (hbar*2*pi)^3
// saves Epot and Ekin into given pointers
// adapted from: T. Yamakoshi, S. Watanabe, Ch. Zhang, and Ch. Greene, "Stochastic and equilibrium pictures of the ultracold Fano-Feshbach-resonance molecular conversion rate", Physical Review A 87, 053604 (2013)
REAL_TYPE f_BEC(REAL_TYPE mass, REAL_TYPE T, REAL_TYPE mu, REAL_TYPE omega[DIM_SPACE], REAL_TYPE x[DIM_PHASESPACE], REAL_TYPE *Epot, REAL_TYPE *Ekin, REAL_TYPE Emin) {
    REAL_TYPE hx = omega[0]*x[0]; // rad/s*mu
    REAL_TYPE hy = omega[1]*x[1];
    REAL_TYPE hz = omega[2]*x[2];
    *Epot = 0.5*mass*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
    *Ekin = 0.5*mass*( x[3]*x[3] + x[4]*x[4] + x[5]*x[5] ) * k_Ekin; // nK
    // check units: inner brackets with mass give [energy/omega] = nK/(rad/s)
    // k_nK = (hbar*1e9)/kB
    // exp gives: nK/((rad/s)*k_nK) = nK/((rad/s)*(hbar*1e9)/kB) = nK/nK = 1
    return exp(-mass/k_nK*((x[0]*x[0]*omega[0] + x[1]*x[1]*omega[1] + x[2]*x[2]*omega[2])*k_Epot + 
                           (x[3]*x[3]/omega[0] + x[4]*x[4]/omega[1] + x[5]*x[5]/omega[2])*k_Ekin ));
}

// calculate chemical potential in nK
// integrate distribution over phase space gives atom number
// in a loop adjust chemical potential to get AtomNumber == N[species]
// this assumes harmonic trap but can be extended to arbitrary potential as well
// TODO: 
// 1. find coordinates giving highest probability. energy = start value.
// 2. find atom_max of each coordinate where 0.1-1.0 atoms are likely for given N. maybe using peak values in other dims.
// 3. allow arbitrary potentials
// 4. use PolyLog

#define TEST_STEPS  250

/* TODO: unfinised! use CalculateChemicalPotential (adapt as in v1.2)

REAL_TYPE CalculateChemicalPotential_adv(uint8_t species, uint8_t stat_index) {
    REAL_TYPE plus_minus_one;
    if      ( stat_index[species] == STAT_FD ) {
        std::cout << "species " << name[species] << " calculating mu (FD) for N = " << N[species] << " at T = " << T[species] << " nK ..." << std::endl;
        plus_minus_one = +1.0;
    }
    else if ((atom_stat_index[species] == STAT_BE   ) || 
             (atom_stat_index[species] == STAT_BE_NC) || 
             (atom_stat_index[species] == STAT_BE_C )) {
        std::cout << "species " << name[species] << " calculating mu (BE) for N = " << N[species] << " at T = " << T[species] << " nK ..." << std::endl;
        plus_minus_one = -1.0;
    }
    else {
        std::cout << "species " << name[species] << " calculating thermal (MB) gas for N = " << N[species] << " at T = " << T[species] << " nK ..." << std::endl;
        return 0.0;
    }
    uint32_t t_start = GetTickCount(), t_next, t_last;
    REAL_TYPE AtomNumber;
    REAL_TYPE x[DIM_PHASESPACE], dx[DIM_PHASESPACE];
    REAL_TYPE dxdp = Mass[species]/k_h; // k_h = h/(mn*1e-9)
    dxdp = 64*dxdp*dxdp*dxdp; // factor 2^6 since all dimensions are symmetric around 0
    for (uint8_t j = 0; j < DIM_PHASESPACE; ++j) {
        dx[j] = atom_max[species][j]/TEST_STEPS;
        dxdp *= dx[j];
    }
    // starting energy = Fermi energy. 
    REAL_TYPE who=std::pow(wx[species]*wy[species]*wz[species],1.0/3.0);
    REAL_TYPE mu = 413.428; //std::pow(6*N[species],1.0/3.0)*who*k_nK; // k_nK = (hbar*1e9)/kB; omega -> nK
    REAL_TYPE StepSize = InitialStepSize*mu;
    REAL_TYPE MaxP = atom_max[species][3]*std::sqrt(3.0);
    dx[3] = MaxP/TEST_STEPS;
    //std::cout << std::scientific << std::setprecision(6); 
    std::cout << "start mu = " << mu << " nK " << dxdp << std::endl;
    bool Direction = true;
    AtomNumber = 0;
    while (fabs((AtomNumber-N[species])/N[species]) > ErrorAllowed) {
        // integrate distribution over all dimensions
        // goal is here not to have the best performing function but to use this to verify PolyLog
        // but this is too slow to get good resolution. 
        // maybe solution: adapt step size for each dimension depending on df/dx.
        //                 to find derivative use peak value in all other dimensions
        //                 using derivative of integral is not symmetric vs. other dimensions.
        // integral can be done unitless
        t_last = GetTickCount();
        AtomNumber = 0.0;
        for (x[0] = 0.0; x[0] < atom_max[species][0]; x[0] += dx[0]) {
            for (x[1] = 0.0; x[1] < atom_max[species][1]; x[1] += dx[1]) {
                for (x[2] = 0.0; x[2] < atom_max[species][2]; x[2] += dx[2]) {
                    for (x[3] = 0.0; x[3] < MaxP; x[3] += dx[3]) {
                        //for (x[4] = 0.0; x[4] < atom_max[species][4]; x[4] += dx[4]) {
                            //for (x[5] = 0.0; x[5] < atom_max[species][5]; x[5] += dx[5]) {
                                REAL_TYPE hx = wx[species]*x[0]; // rad/s*mu
                                REAL_TYPE hy = wy[species]*x[1];
                                REAL_TYPE hz = wz[species]*x[2];
                                REAL_TYPE Epot = 0.5*Mass[species]*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
                                hx = x[3]; // mm/s
                                //hy = x[4];
                                //hz = x[5];
                                REAL_TYPE Ekin = 0.5*Mass[species]*( hx*hx ) * k_Ekin; // nK
                                AtomNumber += 4*Pi/8*hx*hx*1.0/(exp((Epot + Ekin - mu)/T[species]) + plus_minus_one);
                            //}
                        //}
                    }
                    if (((t_next=GetTickCount())-t_last) > 1000000) {
                        std::cout << AtomNumber << std::endl; 
                        t_last = t_next;
                    }
                }
                //std::cout << AtomNumber << std::endl; 
            }
        }
        AtomNumber *= dxdp;
        std::cout << "N(mu = " << mu << ") = " << AtomNumber << std::endl;
        return 0.0;
        // adjust mu to get closer to atom number
        if (AtomNumber > N[species]) {
            if (!Direction) {
                StepSize *= 0.8;
                std::cout<<"*";
            } 
            Direction = true;
            mu -= StepSize;
        } else {
            if (Direction) {
                StepSize *= 0.8;
                std::cout<<"/";
            } 
            Direction = false;
            mu += StepSize;
        }
    }
    std::cout<<std::endl;
    //std::cout << "PreFactor = " << PreFactor << std::endl;
    std::cout << "ChemicalPotential " << name[species] <<" = "<< mu <<" nK (" 
         << get_ticks_delta(t_start) << " ms)" << std::endl;
    return mu;
}
*/

// insert num atoms into shared lists of cells.
// atoms must be sorted by increasing hash_key.
// this way we just need to go through atoms and cells one time simultaneously.
// inserting is done by all threads in parallel, i.e. cells will be locked/unlocked.
// atoms must correspond to cells and must be NULL if nothing should be inserted.
int insert_atoms(class double_linked_list<class species_entry, uint64_t> *atoms, uint8_t cells_index, uint64_t num_atoms, std::string info) {
    int error = 0;
    //for(uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
        if (atoms) {
            if (num_atoms != atoms->get_num()) {
                //std::cout << std::endl << info << " insert_atoms error (1) inconsistent atom number created " << num_atoms <<  " != " << atoms->get_num() << " in list # " << cells_index << " !?" << std::endl << std::endl;
                error = -500;
                //break;
            }
            else {
                //std::cout << info << " insert list of " << atoms->get_num() << " atoms into cell " << +cells_index << " ..." << std::endl;
                class cell_data *cell = ::cells[cells_index]->get_first();
                class species_entry *next = atoms->get_first(), *remove;
                uint64_t num = 0;
                while(next) {

                    // find all atoms with same hash_key
                    uint64_t hash_key = next->get_hash_key(), count = 1;
                    next = next->get_next();
                    while(next && (next->get_hash_key() == hash_key)) {
                        next = next->get_next();
                        ++count;
                    }
                    num += count;
                    
                    if (next != NULL) {
                        // check incrementing hash_key
                        if (next->get_hash_key() < hash_key) {
                            std::cout << info << " insert_atoms error atom " << num << " hash_key " << next->get_hash_key() << " < previous " << hash_key << std::endl;
                            error = -501;
                            break;
                        }
                        next = next->get_prev();
                    }
                    else {
                        next = atoms->get_last();
                    }
                    // remove all atoms from list with same hash_key
                    // returns removed elements starting from first
                    //class species_entry * first = atoms->get_first(), * last = atoms->get_last();
                    remove = atoms->remove_until(next);
                                        
                    // find cell with same hash_key
                    // note: during search cells are not locked although all threads are searching all cells in parallel.
                    //       this is safe since cells and hash_key do not change!
                    //       only the cells->species are changed, which we access only when cell is locked.
                    //       the locking time should be kept at a minimum for best performance.
                    //       add_sorted() cannot be used since: 
                    //       1. is only sorting atoms within cell, and not atoms into different cells.
                    //       2. is not thread-safe and would require to lock cell during search.
                    while(cell) {
                        if (hash_key == cell->get_hash_key()) {
                            // cell found: add count atoms with same hash_key to cell
                            // Attention: add(first, last, num) is adding elements without counting!
                            cell->lock();
                            //cell->species.add(next);
                            cell->species.add(remove, next, count);
                            cell->unlock();
                            break;
                        }
                        cell = cell->get_next();
                    }
                    if (cell == NULL) {
                        // cell not found
                        std::cout << std::endl << info << "insert_atoms error list " << +cells_index << ", atom " << num << ", hash_key " << hash_key << " no cell out of " << ::cells[cells_index]->get_num() <<" cells found! distance measure '" << ::distance_all[cells[cells_index]->get_distance_index()] << "'" << std::endl << "first cell hash_key " << (::cells[cells_index]->get_first() ? ::cells[cells_index]->get_first()->get_hash_key() : 0) << ", last cell hash_key " << (::cells[cells_index]->get_last() ? ::cells[cells_index]->get_last()->get_hash_key() : 0) << std::endl << std::endl;
                        error = -502;
                        // atoms starting from 'next' are already removed from atoms[i] list 
                        // we insert them back into list in order they can be deleted below.
                        atoms->add(next);
                        break;            
                    }
                    
                    // next atom and cell
                    next = atoms->get_first();
                    cell = cell->get_next();
                }
                if (!error) {
                    if (num != num_atoms) {
                        std::cout << std::endl << info << " insert_atoms error (2) inconsistent atom number created " << num <<  " != " << num_atoms << " expected!" << std::endl << std::endl;
                        error = -503;
                    }
                    // ensure all atoms are in cells, i.e. atoms list is empty
                    else if ( (atoms->get_num() != 0) || (atoms->get_first() != NULL) || (atoms->get_last() != NULL) ) {
                        std::cout << std::endl << info << " insert_atoms error not all atoms inserted into cells! remaining " << atoms->get_num() << ", first " << atoms->get_first() << ", last " << atoms->get_last() << std::endl << std::endl;
                        error = -504;
                    }
                }
            }
        }
    //}
    return error;
}


// create distributions
// position in mu, energy in nK, momentum in mm/s
// notes: 
// - this uses rejection-sampling to generate the atoms from the FD distribution.
//   this is extremely inefficient since many atoms (10^3 - 10^5) need to be rejected!
//   to improve this one could use linear or polynomial upper boundaries of the distribution for sampling.
//   CreateAtoms_Metropolis uses the more efficient Metropolis algorithm.
// - in the original implementation the used random number generator (Lehmer32) was inadequate for rejection-sampling
//   due to its short periodicity of 2^32. in addition, to increase the resolution two 32-bit random integers
//   were used to generate one double, which not only requires twice the random numbers but still
//   does not produce doubles equally distributed.  
int CreateAtoms(struct thread_data *pdata, struct cmd_data *cdata) {
    int error = 0;
    uint32_t t_start = GetTickCount();
    uint8_t  id                             = pdata->id;
    class random_generator<REAL_TYPE> *uniform = pdata->gen_uniform;
    uint8_t  species                        = cdata->in_index[0];
    uint8_t  statistics                     = cdata->in_index[1];
    uint64_t num_atoms                      = cdata->total;         // atoms to create
    uint64_t atom_index                     = cdata->count[0];      // atoms index to start
    REAL_TYPE p;
    REAL_TYPE p_max;
    uint64_t count_zero = 0, count_max = 0;
    REAL_TYPE max_x[DIM_PHASESPACE];
    REAL_TYPE Epot, Ekin;
    REAL_TYPE rnd[DIM_PHASESPACE+1];
    //class random_generator<REAL_TYPE> *gen_uniform = pdata->gen_uniform;
    uint64_t loops; // efficiency = number of (created + rejected) atoms / created atoms 
    class species_entry *atom;
    prob_func func;
    std::string info;
    REAL_TYPE omega[DIM_SPACE] = {atom_wx[species],atom_wy[species],atom_wz[species]};
    REAL_TYPE Emin       = 0.5*std::pow(omega[0]*omega[1]*omega[2],1.0/3.0)*k_nK; // ground-state energy in nK. k_nK = (hbar*1e9)/kB
    REAL_TYPE Emax       = 0.0; // maximum allowed energy in nK 
    REAL_TYPE Emax_atoms = 0.0; // maximum energy in nK of generated atoms

    // reset counters    
    cdata->total = 0;
    cdata->count[0] = cdata->count[1] = cdata->count[2] = loops = 0;

    // create double linked-list to collect created atoms for cells.
    // note: at the moment atoms can be inserted only into one list. this might be changed later.
    class double_linked_list<class species_entry, uint64_t> *atoms = nullptr;
    //class species_entry *patoms = nullptr;
    uint8_t cells_index = 0;
    for(; cells_index < NUM_CELLS_LISTS; ++cells_index) {
        if (::cells[cells_index]) {
            if ((::cells[cells_index]->get_type()       == CELLS_TYPE_ATOMS) && 
                (::cells[cells_index]->get_index(0)     == species         ) && 
                (::cells[cells_index]->get_index(1)     == SPECIES_NONE    ) && 
                (::cells[cells_index]->get_statistics() == statistics      )) {
                // allocate memory for all atoms for this thread
                // atoms must be in first 4 lists corresponding to pdata->species index
                if (cells_index >= NUM_ATOMS_LISTS) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << info << "error: cell " << +cells_index << " set for atoms! atoms must be in first " << NUM_ATOMS_LISTS << " cells!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -401;
                    break;
                }
                else if ((pdata->num[cells_index] != 0) || (pdata->mass[cells_index] != 0.0)) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << info << "error: cells " << +cells_index << " already " << pdata->num[cells_index] << " atoms allocated!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -402;
                    break;
                }
                /*patoms = pdata->species[cells_index] = new class species_entry[num_atoms];
                if (patoms == nullptr) {
                    std::cout << info << "error: cells " << +cells_index << " allocation of " << num_atoms << " failed!" << std::endl;
                    error = -403;
                    break;
                }*/
                pdata->num [cells_index] = num_atoms;
                pdata->mass[cells_index] = atom_mass[species];
                // create empty species list for this atom
                atoms = new class double_linked_list<class species_entry, uint64_t>;
                if (atoms == nullptr) {
                    std::cout << info << "error: cells " << +cells_index << " allocation failed!" << std::endl;
                    error = -404;
                    break;
                }
                // use atom name from first cells
                if (info.length() == 0) {
                    statistics = ::cells[cells_index]->get_statistics();
                    for (uint8_t j = 0; j < DIM_PHASESPACE; ++j) max_x[j] = ::cells[cells_index]->get_max(j);
                    Emax = ::cells[cells_index]->get_max(DIM_MAXVAL - 1);
                    info = "thread " + std::to_string(+id) + " rejection-sampling " + ::cells[cells_index]->get_species_name() + " - " + stat_all[statistics] + ": ";
                }
                else {
                    if (id == 0) {
                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                        std::cout << info << "error: cells lists not unique! cells index " << +cells_index << std::endl;
                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    }
                    error = -405;
                    break;
                }
                break;
            }
        }
    }
    
    // check if at least one cell found
    if ((!error) && ((info.length() == 0) || (atoms == nullptr)) ) {
        if (id == 0) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << info << "error: no cell found for species & distribution!" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
        }
        error = -406;
    }
    else {
        // assign probability function
        // note: for STAT_BE we create non-condensed and condensed atoms (STAT_BE_NC + STAT_BE_C)
        if      (statistics == STAT_MB   ) func = f_MB;
        else if (statistics == STAT_FD   ) func = f_FD;
        else if (statistics == STAT_BE_NC) func = f_BE;
        else if (statistics == STAT_BE_C ) func = f_BEC;
        else {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << info << "error: statistics " << stat_all[statistics] << " ( " << +statistics << " ) not implemented!" << std::endl << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            throw ERROR_UNKNOWN_DIST;
        }

        // maximum population is at energy = 0
        for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) rnd[i] = 0.0; 
        p_max = func(::atom_mass[species], atom_T[species], atom_mu[species], omega, rnd, &Epot, &Ekin, Emin);
        if (p_max == 0.0) {
            if (statistics == STAT_BE_NC) {
                // non-condensed bosons must be limited to Emin otherwise f_BE diverges
                p_max = (1.0/(exp(Emin/atom_T[species]) - 1.0));
            }
            else {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << info << "error: distribution function for " << atom_stat_use[species] << " returned zero!" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -407;
            }
        }
    
        if (!error) {
        
            if (id == 0) {
                MUTEX_LOCK(::global_mutex);
                std::cout << info << "creating " << num_atoms << " atoms ..." << std::endl; 
                std::cout << "max population " << p_max << ", Emin = " << Emin << " nK, Emax = " << Emax << " nK" << std::endl;
                MUTEX_UNLOCK(::global_mutex);
            }

            //if (MUTEX.LOCK(::global_mutex)) std::cout << "error lock!";
            //else                std::cout << "thread " << +id << "/" << +pdata->num_threads << " species " << name[species] << " create " << cdata->num[species] << " atoms ... " << std::endl;
            //MUTEX_UNLOCK(::global_mutex);
                
            if (id == 0) {
                std::cout << "max [x, y, z ] = [" << max_x[0] << "," << max_x[1] << "," << max_x[2] << "] mu" << std::endl;
                std::cout << "max [vx,vy,vz] = [" << max_x[3] << "," << max_x[4] << "," << max_x[5] << "] mm/s" << std::endl;
            }

            for (uint64_t num = 0; num < num_atoms; ++num) {
                bool AtomCreated = false;
                
                // allocate atom. update: get pointer to atom species_data
                // TODO: this should throw an exception when out of memory!?
                atom = new class species_entry;
                //atom = &patoms[num].species;
                            
                while (!AtomCreated) {
                    ++loops;
                    
                    (*uniform)(rnd, DIM_PHASESPACE+1);
                    //Random_thread_safe(pdata->seed, &rnd[6], 1, true);
                    rnd[0] = (1.0-2.0*rnd[0])*max_x[0]; // um
                    rnd[1] = (1.0-2.0*rnd[1])*max_x[1];
                    rnd[2] = (1.0-2.0*rnd[2])*max_x[2];
                    rnd[3] = (1.0-2.0*rnd[3])*max_x[3]; // mm/s
                    rnd[4] = (1.0-2.0*rnd[4])*max_x[4];
                    rnd[5] = (1.0-2.0*rnd[5])*max_x[5];
                    
                    p = func(::atom_mass[species], atom_T[species], atom_mu[species], omega, rnd, &Epot, &Ekin, Emin)/p_max;
                    if (p <= 0.0) ++count_zero;
                    else {
                        AtomCreated=(rnd[6] < p);
                        //if ((AtomCreated) && (id==0) && (statistics==STAT_BE_C)) std::cout << atom_index << " / " << num_atoms << std::endl;
                    }
                }
                atom->species.mass = ::atom_mass[species];
                for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                    atom->species.x[i] = rnd[i];
                }
                SET_STATUS_FREE(atom->species.status); 
                atom->species.E[E_POT] = Epot;
                atom->species.E[E_KIN] = Ekin;
                atom->species.E[E_RED] = 0.0;
                atom->species.distance = 0.0;
                atom->species.atom_index[0] = atom_index++;

                // E > Emax when max_ values are too small!
                // in this case calculation is not reliable and you must increase calc_size variables!
                REAL_TYPE tmp = Epot + Ekin;
                if (tmp > Emax_atoms) Emax_atoms = tmp;
                if ((p >= 1.0) || (tmp > Emax)) ++count_max;
                            
                // add atom to atoms lists using cells distance measure to generate hash_key
                //for(uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                //    if (atoms[i]) {
#if NUM_CELLS > 1
                        //atoms[i]->add_sorted(new class species_entry(atom, NULL, ::cells[i]->get_hash_key(atom))); 
                        atom->hash_key = ::cells[cells_index]->get_hash_key(atom);
                        atoms->add_sorted(atom); 
#else
                        //atoms[i]->add(new class species_entry(atom, NULL, id));
                        atom->hash_key = id;
                        atoms->add(atom); 
#endif
                //    }
                //}
            } // next atom
            
            if (!error) {
#if NUM_CELLS > 1
                // insert atoms into shared lists of cells. this clears atoms list.
                error = insert_atoms(atoms, cells_index, num_atoms, info);
#else
                // insert atoms into linear arrays of thread data structure of owner thread
                // molecule search is significantly faster with this array
                REAL_TYPE *x[DIM_PHASESPACE];
                REAL_TYPE *E[DIM_ENERGY];
                uint8_t *status = pdata->status[cells_index] = new uint8_t[num_atoms];
                pdata->x[cells_index] = new REAL_TYPE*[DIM_PHASESPACE];
                for(uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
                    x[d] = pdata->x[cells_index][d] = new REAL_TYPE[num_atoms];
                }
                pdata->E[cells_index] = new REAL_TYPE*[DIM_ENERGY];
                for(uint8_t d = 0; d < DIM_ENERGY; ++d) {
                    E[d] = pdata->E[cells_index][d] = new REAL_TYPE[num_atoms];
                }
                class species_entry *atom = atoms->get_first();
                uint64_t c = 0;
                while(atom) {
                    for(uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
                        *x[d]++ = atom->species.x[d];
                    } 
                    for(uint8_t d = 0; d < DIM_ENERGY; ++d) {
                        *E[d]++ = atom->species.E[d];
                    } 
                    SET_STATUS_FREE(*status++);
                    atom = atom->get_next();
                    ++c;
                }
                if (c != num_atoms) {
                    std::cout << "error num_atoms " << c << " != " << num_atoms << " inserted into linear array!" << std::endl;
                    error = -475;
                }
                
                // insert atoms into shared lists of cells. this clears atoms list.
                // TODO: might not be needed here. maybe just need num_atoms, owner and mass
                error = insert_atoms(atoms, cells_index, num_atoms, info);
#endif 
                if (!error) {
                    MUTEX_LOCK(::global_mutex);
                    std::cout << info << num_atoms << "/" << loops << " atoms/loops = 1:" << (uint64_t)std::ceil(((REAL_TYPE)loops)/num_atoms) << " (" << (((REAL_TYPE)get_ticks_delta(t_start)) / 1e3) << " s)" << std::endl;
                    MUTEX_UNLOCK(::global_mutex);

                    // sum up total number of created atoms
                    cdata->total = num_atoms;
                    // save loop counter for efficiency calculation
                    cdata->count[0] = loops;
                    // save E < Emin counter for non-condensed bosons
                    cdata->count[1] = count_zero;
                    // save E > Emax counter. this indicates that max_ values are too small!
                    cdata->count[2] = count_max;
                    // save maximum energy of generated atoms
                    cdata->E[E_MAX] = Emax_atoms;
                }
            }
        }
        
        // now atoms list should be empty and can be deleted.
        // on error there might be remaining atoms entries which we remove first.
        // atoms data is never deleted since atoms might be already inserted into shared cells.
        // they can only be deleted at end of each repetition after cells have been deleted.
        //for(int8_t i = NUM_CELLS_LISTS - 1; i >= 0; --i) {
        //    if (atoms[i]) {
                if (atoms->get_num() > 0) {
                    //std::cout << (::cells[i]->get_duplicate() ? "clearing " : "deleting ") << ::atom_name[species] << " list (" << +cells_index << ") " << atoms->get_num() << " entries ..." << std::endl;
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: not all atoms inserted into cells " << ::atom_name[species] << " list (" << +cells_index << ")! remaining " << atoms->get_num() << " entries. is there an error (" << error << ") ?" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    class species_entry *next = atoms->get_first();
                    while(next) {
                        //std::cout << "first " << atoms[i]->get_first() << " next " << next << std::endl; sleep_ms(10);
                        atoms->remove(next);
                        //std::cout << "first " << atoms[i]->get_first() << " next " << next << std::endl; sleep_ms(10);
                        /*
                        if (!::cells[i]->get_duplicate()) {
                            class species_data *sp = next->get_species(0);
                            if (sp) delete sp;
                            sp = next->get_species(1); // should be NULL
                            if (sp) delete sp;
                        }
                        delete next;
                        */
                        next = atoms->get_first();
                    }
                }
                delete atoms;
            //}
        //}
        
        //if(error && patoms) delete[] patoms; 
    }        
    
    return error;
}

// number of samples to skip at the beginning and from sample to sample to avoid correlations
// the efficiency is given mainly by METROPOLIS_SKIP
// note: to keep every sample for testing it is fine to set these values to 0 for testing.
#define METROPOLIS_SKIP_START   1000
#define METROPOLIS_SKIP         100

// initial starting sigma in units of maximum values
REAL_TYPE sigma_start = 0.1;

// ideal rejection rate is 0.25-0.50. sigma is adjusted to match this.
// the rejection rate does not affect the sampling efficiency but defines how fast the state changes.
// larger values smear out distribution (longer high-Etot tail)
// too small values prevent from exploring low-probability regions and 
// to recover from bad starting values or outlying states (with 0.05 saw one time that distribution in p is shifted)
#define METROPOLIS_REJECT_IDEAL     0.25

// number of loops such that sigma changes from 0 to 1 for a large error
// gain is inverse proportional to this
// for METROPOLIS_REJECT_IDEAL = 0.2 - 0.8 (for 0.5 is always closer to ideal)
// 10k is not sufficient
// 1k  is about 1% within ideal
// 100 is about 0.5% within ideal
// 10  is about 0.2% within ideal
#define METROPOLIS_GAIN_LOOPS       100

// create distributions with Metropolis algorithm which is much more efficienct than CreateAtoms above
// uses statistics as specified for species in stat_index and stat_use 
int CreateAtoms_Metropolis(struct thread_data *pdata, struct cmd_data *cdata) {
    int error = 0;
    uint32_t t_start = GetTickCount();
    if      ((cdata == NULL) || (pdata == NULL)   ) return -450;
    else if (cdata->in_index[0] >= MAX_NUM_SPECIES) return -451;
    else if (cdata->total == 0                    ) return 0; // no atoms to create (happens when Nc or Nthermal < num_threads)
    uint8_t id                                  = pdata->id;
    class random_generator<REAL_TYPE> *normal   = pdata->gen_normal;
    class random_generator<REAL_TYPE> *uniform  = pdata->gen_uniform;
    uint8_t  species            = cdata->in_index[0];
    uint8_t  statistics         = cdata->in_index[1];
    uint64_t num_atoms          = cdata->total;
    uint64_t atom_index         = cdata->count[0];
    REAL_TYPE mass              = ::atom_mass[species];
    REAL_TYPE T                 = ::atom_T   [species];
    REAL_TYPE mu                = ::atom_mu  [species];
    REAL_TYPE omega[DIM_SPACE]  = {::atom_wx[species],::atom_wy[species],::atom_wz[species]};
    uint64_t num                = 0;
    uint64_t loops              = 0;
    uint64_t count_zero         = 0;
    uint64_t count_max          = 0;
    uint64_t accept[DIM_PHASESPACE];
    REAL_TYPE rnd[DIM_PHASESPACE];

    // minimum and (maximum-minimum) values for each axis
    //REAL_TYPE min  [DIM_PHASESPACE];
    //REAL_TYPE delta[DIM_PHASESPACE];

    // sigma as a fraction from maximum
    REAL_TYPE sigma[DIM_PHASESPACE];

    // state vector of each coordinate    
    REAL_TYPE vx[DIM_PHASESPACE], vx_test[DIM_PHASESPACE];
    
    REAL_TYPE Epot, Ekin;
    std::string info;
    REAL_TYPE Emin       = 0.5*std::pow(omega[0]*omega[1]*omega[2],1.0/3.0)*k_nK; // ground-state energy in nK
    REAL_TYPE Emax       = 0.0; // maximum energy in nK used to check if size is ok.
    REAL_TYPE Emax_atoms = 0.0; // generated atoms maximum energy in nK
        
    //std::cout << "CreateAtoms_Metropolis species " << ::atom_name[species] << " creating " << num_atoms << " atoms ... " << std::endl;

    // reset counters    
    cdata->total    = 0;
    cdata->count[0] = cdata->count[1] = cdata->count[2] = 0;

    // create double linked-list to collect created atoms for cells.
    // TODO: cdata->out_index[0] = index into cells for output not used now.
    class double_linked_list<class species_entry, uint64_t> *atoms = nullptr;
    //class species_entry *patoms = nullptr;
    uint8_t cells_index = 0;
    for(; cells_index < NUM_CELLS_LISTS; ++cells_index) {
        if (::cells[cells_index]) {
            if ((::cells[cells_index]->get_type()       == CELLS_TYPE_ATOMS) &&
                (::cells[cells_index]->get_index(0)     == species         ) &&
                (::cells[cells_index]->get_index(1)     == SPECIES_NONE    ) && 
                (::cells[cells_index]->get_statistics() == statistics      )) {
                // allocate memory for all atoms for this thread
                // atoms must be in first 4 lists corresponding to pdata->species index
                if (cells_index >= NUM_ATOMS_LISTS) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << info << "error: cells " << +cells_index << " set for atoms! atoms must be in first 4 cells!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -451;
                    break;
                }
                else if ((pdata->num[cells_index] != 0) || (pdata->mass[cells_index] != 0.0)) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << info << "error: cells " << +cells_index << " already " << pdata->num[cells_index] << " atoms allocated!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -452;
                    break;
                }
                pdata->num [cells_index] = num_atoms;
                pdata->mass[cells_index] = mass;
                // create empty species list for this atom
                atoms = new class double_linked_list<class species_entry, uint64_t>;
                if (atoms == nullptr) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << info << "error: cells " << +cells_index << " allocation failed!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -454;
                    break;
                }
                // use atom name from first cells
                if (info.length() == 0) {
                    statistics = ::cells[cells_index]->get_statistics();
                    info = "thread " + std::to_string(+id) + " Metropolis " + ::cells[cells_index]->get_species_name() + " - " + stat_all[statistics] + ": ";
                    // init all vectors
                    for (uint8_t j = 0; j < DIM_PHASESPACE; ++j) {
                        accept[j] = 0;
                        vx    [j] = 0;
                        //min   [j] = ::cells[cells_index]->get_max(j) * (-1);
                        //delta [j] = ::cells[cells_index]->get_max(j) * 2;
                        sigma [j] = ::cells[cells_index]->get_max(j) * sigma_start;
                    }
                    // get maximum energy
                    Emax = ::cells[cells_index]->get_max(DIM_MAXVAL-1);
                }
                else {
                    if (id == 0) {
                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                        std::cout << info << "error: cells lists not unique! cells index " << +cells_index << std::endl;
                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    }
                    error = -455;
                    break;
                }
#ifdef _DEBUG
                // check if cells are empty: num_threads must be 1 otherwise other threads fill already cells!
                if (num_threads == 1) {
                    num = 0;
                    class cell_data *cell = ::cells[cells_index]->get_first();
                    while(cell) {
                        num += cell->species.get_num();
                        cell = cell->get_next();
                    }
                    if (num != 0) {
                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                        std::cout << info << "error: cells index " << +cells_index << " not empty! num = " << num << std::endl << std::endl;
                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                        error = -456;
                    }
                }
#endif          
                break;  
            }
        }
    }
    // check if at least one cell found
    if ((!error) && ((atoms==nullptr) || (info.length() == 0))) {
        if (id == 0) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << info << "error: no cells found for species & istribution!" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
        }
        error = -460;
    }
    else {
        // gain for adjustment of sigma
        // sigma changes completely within 100 x METROPOLIS_SKIP loops
        REAL_TYPE gain = 1.0/METROPOLIS_GAIN_LOOPS;
        
        // get probability of start distribution
        // note: for STAT_BE we create non-condensed and condensed atoms (STAT_BE_NC + STAT_BE_C)
        REAL_TYPE p;
        prob_func func;
        if      (statistics == STAT_MB   ) func = f_MB;
        else if (statistics == STAT_FD   ) func = f_FD;
        else if (statistics == STAT_BE_NC) func = f_BE;
        else if (statistics == STAT_BE_C ) func = f_BEC;
        else {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << info << "error: statistics not implemented! (index " << +statistics << " )" << std::endl << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            throw ERROR_UNKNOWN_DIST;
        }

        // start with zero energy probability
        p = func(mass, T, mu, omega, vx, &Epot, &Ekin, Emin);
        if (p < 1e-10) {
            if (statistics == STAT_BE_NC) {
                // non-condensed bosons must be limited to Emin otherwise f_BE diverges
                p = (1.0/(exp(Emin/atom_T[species]) - 1.0));
            }
            else {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << info << "error: distribution function returned zero!" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -461;
            }
        }
    
        if (!error) {
            // loop until we have all atoms generated
            uint32_t cskip = METROPOLIS_SKIP_START + 1;
            while(num < num_atoms) {
                ++loops;
                for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                    // calculate new state and map inside [min,max]
                    //vx_test[i] = min[i] + std::fmod(((*normal)()*sigma[i] - min[i]), delta[i]);
                    vx_test[i] = vx[i] + (*normal)()*sigma[i];
                }
                // calculate probability of state
                REAL_TYPE pt = func(mass, T, mu, omega, vx_test, &Epot, &Ekin, Emin);
                if (pt < 0.0) ++count_zero;
                else {
                    REAL_TYPE ratio = pt/p;
                    // accept or reject new state on each axis
                    // if new state has higher probability than older it will be always accepted
                    // if not, it will be still sometimes accepted when its not too unlikely
                    // note: regardless of state acceptance or rejectance in this step the state might be saved anyway.
                    (*uniform)(rnd, DIM_PHASESPACE);
                    for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                        if (rnd[i] <= ratio) {
                            ++accept[i];
                            vx[i] = vx_test[i];
                        }
                    }
                    //std::cout << info << "loop " << loops << " x " << vx_test[0] << " pt " << pt << " ratio " << ratio << " rnd " << rnd[0] << std::endl; 
                    // recalculate probability of new state
                    p = func(mass, T, mu, omega, vx, &Epot, &Ekin, Emin);
                    if (p == 0.0) ++count_zero;
                    else if (--cskip == 0) {
                        // save sample
                        cskip = METROPOLIS_SKIP + 1;
                        class species_entry *atom = new class species_entry;
                        //class species_data *atom = &patoms[num].species;
                        atom->species.mass = mass;
                        for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                            atom->species.x[i] = vx[i];
                        }
                        SET_STATUS_FREE(atom->species.status); 
                        atom->species.E[E_POT]      = Epot;
                        atom->species.E[E_KIN]      = Ekin;
                        atom->species.E[E_RED]      = 0.0;
                        atom->species.distance      = 0.0;
                        atom->species.atom_index[0] = atom_index++;
                        
                        // this can happen when max_ values are too small
                        // note: Metropolis neither needs max_ values nor normalization of p.
                        REAL_TYPE tmp = Epot + Ekin;
                        if (tmp > Emax) ++count_max;
                        if (tmp > Emax_atoms) Emax_atoms = tmp;

                        // add atom to atoms lists using requested distance measure to generate hash_key
                        //for(uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                        //    if (atoms[i]) {
#if NUM_CELLS > 1                        
                                //atoms[i]->add_sorted(new class species_entry(atom, NULL, ::cells[i]->get_hash_key(atom)));
                                atom->hash_key = ::cells[cells_index]->get_hash_key(atom);
                                atoms->add_sorted(atom);
#else
                                //atoms[i]->add(new class species_entry(atom, NULL, id));
                                atom->hash_key = id;
                                atoms->add(atom);
#endif
                        //    }
                        //}
                        
                        // adjust sigma of generator
                        for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                            sigma[i] += gain*sigma[i]*(METROPOLIS_REJECT_IDEAL - (1.0 - ((REAL_TYPE)accept[i])/loops));
                        }

                        // next atom
                        ++num;
                    }
                }
            }
        }
        
        if (!error) {    
            if (num != num_atoms) {
                std::cout << info << "error (1) inconsistent atom number created " << num <<  " != " << num_atoms << " expected!" << std::endl;
                error = -470;
            }
        }

        if (error) {
            cdata->total = 0;
        }
        else {
            // output of efficiency and rejection rate
            // efficiency     = 1 : (loops/N) = 1 : ( METROPOLIS_SKIP + METROPOLIS_SKIP_START/N )
            // rejection rate = sum rejected/( DIM_PHASESPACE*loops )
            // efficiency is indepenent of rejection rate since any atom is taken regardless if was accepted or rejected.
            // rejection rate gives how fast state changes
            uint64_t rejected = 0;
            for (uint8_t i = 0; i < DIM_PHASESPACE; ++i) {
                rejected += loops-accept[i];
                //std::cout << "axis " << +i << " : sigma " << sigma[i] << " rejected " << loops-accept[i] << " = " << (1.0-((REAL_TYPE)accept[i])/loops) << std::endl;
            }
            
#ifdef _DEBUG
            // check list of atoms
            uint64_t c = 0;
            class species_entry *next = atoms->get_first(), *last = nullptr; 
            while(next) {
                if (next->get_prev() != last) {
                    std::cout << "error prev " << last << " != " << next->get_prev() << std::endl;
                }
                last = next;
                next = next->get_next();
                ++c;
            }
            if ((c != num_atoms) || (atoms->get_num() != num_atoms)) {
                std::cout << "error num_atoms " << c << " != " << num_atoms << ", " << atoms->get_num() << std::endl;
                error = -471;
            }
            else if (last != atoms->get_last()) {
                std::cout << "error last " << last << " != " << atoms->get_last() << std::endl;
                error = -472;
            }
#endif

            if (!error) {
#if NUM_CELLS > 1
                // insert atoms into shared lists of cells. this clears atoms list.
                error = insert_atoms(atoms, cells_index, num_atoms, info);
#else
                // insert atoms into linear arrays of thread data structure of owner thread
                // molecule search is significantly faster with this array
                REAL_TYPE *x[DIM_PHASESPACE];
                REAL_TYPE *E[DIM_ENERGY];
                uint8_t *status = pdata->status[cells_index] = new uint8_t[num_atoms];
                pdata->x[cells_index] = new REAL_TYPE*[DIM_PHASESPACE];
                for(uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
                    x[d] = pdata->x[cells_index][d] = new REAL_TYPE[num_atoms];
                }
                pdata->E[cells_index] = new REAL_TYPE*[DIM_ENERGY];
                for(uint8_t d = 0; d < DIM_ENERGY; ++d) {
                    E[d] = pdata->E[cells_index][d] = new REAL_TYPE[num_atoms];
                }
                class species_entry *atom = atoms->get_first();
                uint64_t c = 0;
                while(atom) {
                    for(uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
                        *x[d]++ = atom->species.x[d];
                    } 
                    for(uint8_t d = 0; d < DIM_ENERGY; ++d) {
                        *E[d]++ = atom->species.E[d];
                    } 
                    SET_STATUS_FREE(*status);
                    atom = atom->get_next();
                    ++status;
                    ++c;
                }
                if (c != num_atoms) {
                    std::cout << "error num_atoms " << c << " != " << num_atoms << " inserted into linear array!" << std::endl;
                    error = -475;
                }
                
                // insert atoms into shared lists of cells. this clears atoms list.
                // TODO: might not be needed here. maybe just need num_atoms, owner and mass
                error = insert_atoms(atoms, cells_index, num_atoms, info);
#endif 
            }
                        
            if (!error) {
                MUTEX_LOCK(::global_mutex);
                //const auto default_precision{std::std::cout.precision()};
                std::cout << std::fixed << std::setprecision(1);
                std::cout << info << num_atoms
                     << " created 1:" << (int)std::ceil(((REAL_TYPE)loops)/num_atoms) 
                     << ", rejected " << ((REAL_TYPE)rejected)*100.0/(DIM_PHASESPACE*loops) << " % "; 
                if (error) std::cout << "error " << error << " !" << std::endl;
                else       std::cout << "ok" << " (" << (((REAL_TYPE)get_ticks_delta(t_start)) / 1e3) << " s)" << std::endl;
                std::cout << std::fixed << std::setprecision(OUTPUT_PRECISION);
                MUTEX_UNLOCK(::global_mutex);
                
                // total number of created atoms
                cdata->total = num_atoms;
                // save loop counter for efficiency calculation
                cdata->count[0] = loops;
                // save E < Emin counter for non-condensed bosons
                cdata->count[1] = count_zero;
                // save E > Emax counter
                cdata->count[2] = count_max;
                // save maximum energy of created atoms
                cdata->E[E_MAX] = Emax_atoms;
            }            
        }

        // now atoms list should be empty and can be deleted.
        // on error there might be remaining atoms entries which we remove first.
        // atoms data is never deleted since atoms might be already inserted into shared cells.
        // they can only be deleted at end of each repetition after cells have been deleted.
        //for(int8_t i = NUM_CELLS_LISTS - 1; i >= 0; --i) {
        //    if (atoms[i]) {
                if (atoms->get_num() > 0) {
                    //std::cout << (::cells[i]->get_duplicate() ? "clearing " : "deleting ") << ::atom_name[species] << " list (" << +cells_index << ") " << atoms->get_num() << " entries ..." << std::endl;
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: not all atoms inserted into cells " << ::atom_name[species] << " list (" << +cells_index << ")! remaining " << atoms->get_num() << " entries. is there an error (" << error << ") ?" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    class species_entry *next = atoms->get_first();
                    while(next) {
                        //std::cout << "first " << atoms[i]->get_first() << " next " << next << std::endl; sleep_ms(10);
                        atoms->remove(next);
                        //std::cout << "first " << atoms[i]->get_first() << " next " << next << std::endl; sleep_ms(10);
                        /*
                        if (!::cells[i]->get_duplicate()) {
                            class species_data *sp = next->get_species(0);
                            if (sp) delete sp;
                            sp = next->get_species(1); // should be NULL
                            if (sp) delete sp;
                        }
                        delete next;
                        */
                        next = atoms->get_first();
                    }
                }
                delete atoms;
            //}
        //}
        
        //if(error && patoms) delete[] patoms; 

    }    
    
    return error;
}

// search molecules in atoms in a linear fashion
// this is the older version v1.2 code which is very fast! 
// TODO: new code needs to be tested against this code to ensure not getting utterly slow (using cells is bad).
// stops the search when it finds the first matching atom 1 for atom 0 and does not look for nearest neighbors.
// gives in average the same result for any number of threads, but different solutions might be found from shot to shot,
// even for the same seed value. for large atom number can take very long time!
// note: this uses global mutex! not pdata->mutex.
int CreateMolecules_v1_2(struct thread_data *pdata, struct cmd_data *cdata) {
    int error = 0;
    REAL_TYPE hx, hy, hz;    
    uint8_t id        = pdata->id;           // thread id
    uint8_t sp0       = cdata->in_index [0]; // cell index of atom 0
    uint8_t sp1       = cdata->in_index [1]; // cell index of atom 1
    uint8_t mol_index = cdata->out_index[0]; // cell index of molecule. will be allocated here.
    long num_lists    = pdata->num_threads;
    uint64_t molecules = 0;
    long failed[3] = {0, 0, 0}; // count how often a possible pair is already paired by another thread. rarely happens.

    REAL_TYPE Mass[2], Mmol;
    
    REAL_TYPE *x0    = pdata->x     [sp0][0];
    REAL_TYPE *y0    = pdata->x     [sp0][1];
    REAL_TYPE *z0    = pdata->x     [sp0][2];
    REAL_TYPE *px0   = pdata->x     [sp0][3];
    REAL_TYPE *py0   = pdata->x     [sp0][4];
    REAL_TYPE *pz0   = pdata->x     [sp0][5];
    uint8_t *status0 = pdata->status[sp0];
    uint64_t num0    = pdata->num   [sp0];
    Mass[0]          = pdata->mass  [sp0]; 
        
    REAL_TYPE PotentialEnergyOfMolecules = 0.0;
    REAL_TYPE KineticEnergyOfMolecules   = 0.0;
    REAL_TYPE ReductionOfEnergy          = 0.0;

    char old_status[2];
    bool not_found;
    
    REAL_TYPE gamma = ::cells[sp0]->get_gamma();

    // pre-allocate maximum number of molecules == num0
    // note: for single species it would be min(num0,atom_N[0]/2) but num0=atom_N[0]/number_threads.
    //       for simplicity we take num0 which is a safe values in all conditions.
    pdata->x     [mol_index] = new REAL_TYPE*[DIM_PHASESPACE];
    for (uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
        pdata->x [mol_index][d] = new REAL_TYPE[num0];
    }
    // neither energy nor status are used!
    pdata->E     [mol_index] = nullptr;
    pdata->status[mol_index] = nullptr;
    REAL_TYPE *xmol  = pdata->x[mol_index][0];
    REAL_TYPE *ymol  = pdata->x[mol_index][1];
    REAL_TYPE *zmol  = pdata->x[mol_index][2];
    REAL_TYPE *pxmol = pdata->x[mol_index][3];
    REAL_TYPE *pymol = pdata->x[mol_index][4];
    REAL_TYPE *pzmol = pdata->x[mol_index][5];
    REAL_TYPE EpotMol, EkinMol;
    
    std::string info = "thread " + std::to_string(+id) + " create molecules (v1.2): ";
                
    MUTEX_LOCK(::global_mutex);
    std::cout << info << "searching " << num0 << " atoms, " << ((sp0 == sp1) ? "1" : "2") << " species, gamma " << gamma << ", mol index " << +mol_index << std::endl;
    MUTEX_UNLOCK(::global_mutex);

    // status indicating atom 0 is unpaired. this must be different for each combination of species
    //const uint8_t checked_code = MAKE_CHECKED_CODE(cdata->out_index[1]);
    
    // go through all atoms of species 0
    for (uint64_t atom0 = 0; (atom0 < num0) && (!error); ++atom0) {
    
        // for num_species == 1 we need to check if atom0 is already paired
        if ( IS_STATUS_FREE(*status0) ) {
                
            // go through all lists
            not_found = true;
            struct thread_data *list = pdata;
            for (long list_count = 0; not_found && (list_count < num_lists); ++list_count, list = list->next) {
            
                // for same species take same pool of atoms species 0, otherwise species 1
                //int offset = (num_species == 1) ? 0 : 1;

                // assign new atom1 pointers from new list

                // initialize atom1 pointers
                REAL_TYPE *x1    = list->x     [sp1][0];
                REAL_TYPE *y1    = list->x     [sp1][1];
                REAL_TYPE *z1    = list->x     [sp1][2];
                REAL_TYPE *px1   = list->x     [sp1][3];
                REAL_TYPE *py1   = list->x     [sp1][4];
                REAL_TYPE *pz1   = list->x     [sp1][5];
                uint8_t *status1 = list->status[sp1];
                uint64_t num1    = list->num   [sp1];
                Mass[1]          = list->mass  [sp1];
                Mmol = Mass[0] + Mass[1]; 

                // go through all atoms of species 1
                for (uint64_t atom1 = 0; atom1 < num1; ++atom1) {
                    // take free atoms unless its the same (if num_species == 1)
                    if ( ( IS_STATUS_FREE(*status1) ) && (status0 != status1) ) {
                        // calculate phase space distance
                        // note: for improved efficiency we use squared psd to avoid sqrt. result is the same.

                        hx = (*x0) - (*x1)- atom_cx; // um
                        hy = (*y0) - (*y1)- atom_cy;
                        hz = (*z0) - (*z1)- atom_cz;
                        REAL_TYPE RealSpaceDistance = hx*hx + hy*hy + hz*hz; // um^2
                        REAL_TYPE PhaseSpaceDistance;
                        
                        if (mol_distance_index == DISTANCE_X) {
                            // real-space distance
                            PhaseSpaceDistance = RealSpaceDistance;
                        }
                        else if (mol_distance_index == DISTANCE_PP) {
                            // lab frame momentum-space distance
                            hx = (*px0)*Mass[0] - (*px1)*Mass[1]; // mm/s * amu
                            hy = (*py0)*Mass[0] - (*py1)*Mass[1];
                            hz = (*pz0)*Mass[0] - (*pz1)*Mass[1];
                            REAL_TYPE MomentumSpaceDistance = hx*hx + hy*hy + hz*hz; // (mm/s * amu)^2
                            // [dx*dp] = kg*m^2/s = Js
                            PhaseSpaceDistance = RealSpaceDistance*MomentumSpaceDistance; // um^2*(mm/s)^2 = (m^2/s*1e-9)^2
                        }
                        else if (mol_distance_index == DISTANCE_PV ){
                            // use velocity-space distance
                            //Want to substract velocities in center of mass frame
                            //This is equivalent to simply substracting the velocities
                            hx = (*px0) - (*px1); // mm/s
                            hy = (*py0) - (*py1);
                            hz = (*pz0) - (*pz1);
                            //To convert relative velocity to momentum in center of mass frame, use reduced mass
                            REAL_TYPE MomentumSpaceDistance = hx*hx + hy*hy + hz*hz; // (mm/s)^2
                            // [dx*dp] = kg*m^2/s = Js
                            PhaseSpaceDistance = RealSpaceDistance*MomentumSpaceDistance; // um^2*(mm/s)^2 = (m^2/s*1e-9)^2
                        }
                        else {
                            // error: illegal distance measure
                            // note : this code does not have all distance measures implemented!
                            error = -710;
                            not_found = false;
                            break;
                        }
                        
                        if (PhaseSpaceDistance < gamma) {
                            // possible partner found: if both atoms are still free then pair them.
                            
                            MUTEX_LOCK(::global_mutex);
                            old_status[0] = *status0;
                            old_status[1] = *status1;
                            if ( IS_STATUS_FREE( old_status[0] ) && IS_STATUS_FREE( old_status[1] ) ) {
                                SET_STATUS_PAIRED(*status0);
                                SET_STATUS_PAIRED(*status1);
                            }
                            MUTEX_UNLOCK(::global_mutex);
                            
                            if ( IS_STATUS_FREE( old_status[0] ) && IS_STATUS_FREE( old_status[1] ) ) {
                                // pair found
                                not_found = false;
                                if (++molecules > num0) {
                                    error = -711;
                                    break;
                                }
                                
                                //Statistics on molecules
                                //std::cout << "pair found " << atom0 << " with atom " << atom1 << ": psd = " << PhaseSpaceDistance << std::endl;

                                //Center of mass of molecule: 
                                *xmol = (Mass[0]*(*x0)+Mass[1]*((*x1)+atom_cx))/Mmol; // um
                                *ymol = (Mass[0]*(*y0)+Mass[1]*((*y1)+atom_cy))/Mmol;
                                *zmol = (Mass[0]*(*z0)+Mass[1]*((*z1)+atom_cz))/Mmol;

                                REAL_TYPE Epot[2];
                                REAL_TYPE Ekin[2];
                                hx = atom_wx[0]*(*x0); // rad/s*mu
                                hy = atom_wy[0]*(*y0);
                                hz = atom_wz[0]*(*z0);
                                Epot[0] = 0.5*Mass[0]*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
                                hx = atom_wx[1]*(*x1); // rad/s*mu
                                hy = atom_wy[1]*(*y1);
                                hz = atom_wz[1]*(*z1);
                                Epot[1] = 0.5*Mass[1]*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
                                hx = (*px0);
                                hy = (*py0);
                                hz = (*pz0);
                                Ekin[0] = 0.5*Mass[0]*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
                                hx = (*px1);
                                hy = (*py1);
                                hz = (*pz1);
                                Ekin[1] = 0.5*Mass[1]*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
                                
                                //Momentum conservation gives momentum of molecule after association
                                *pxmol = hx = ((*px0)*Mass[0]+(*px1)*Mass[1])/Mmol; // mm/s
                                *pymol = hy = ((*py0)*Mass[0]+(*py1)*Mass[1])/Mmol;
                                *pzmol = hz = ((*pz0)*Mass[0]+(*pz1)*Mass[1])/Mmol;

                                EpotMol = Epot[0] + Epot[1];
                                EkinMol = 0.5*Mmol*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK

                                ReductionOfEnergy += Ekin[0] + Ekin[1] - EkinMol; // nK
                                PotentialEnergyOfMolecules += EpotMol; // nK
                                KineticEnergyOfMolecules   += EkinMol;
                                ++xmol;  ++ymol;  ++zmol;
                                ++pxmol; ++pymol; ++pzmol;
                                
                                break;
                            }
                            else if (IS_STATUS_PAIRED(old_status[0])) {
                                // atom 0 paired by another thread: stop searching atom 0 & 1
                                ++failed[0]; // can happen
                                not_found = false;
                                break;
                            }
                            else if (IS_STATUS_PAIRED(old_status[1])) {
                                // atom 1 paired by another thread: continue searching atom 1
                                ++failed[1]; // never seen?
                            }
                            else {
                                // error cases:
                                // 1. atom 0 unpaired can only be set by this thread.
                                // 2. atom 0 free and atom 1 unpaired although they could be paired.
                                error = -712;
                                not_found = false;
                                break;
                            }
                        }
                    }
                    
                    // increment atom 1 pointers
                    ++x1; ++y1; ++z1;
                    ++px1; ++py1; ++pz1;
                    ++status1;
                } // next atom species 1
            } // next list
            
            if ( not_found && (!error) ) {
                // mark atom 0 if no pair found
                MUTEX_LOCK(::global_mutex);
                old_status[0] = *status0;
                if (IS_STATUS_FREE(old_status[0])) SET_STATUS_UNPAIRED(*status0);
                MUTEX_UNLOCK(::global_mutex);

                if (IS_STATUS_PAIRED(old_status[0])) {
                    // no pair found by this thread but another thread could pair it in the mean time
                    ++failed[2]; // can happen
                }
                else if (IS_STATUS_UNPAIRED(old_status[0])) {
                    // error: atom 0 unpaired can only be set by this thread.
                    error = -713;
                    not_found = false;
                    break;
                }
            }
        }
        
        // increment atom 0 pointers
        ++x0; ++y0; ++z0;
        ++px0; ++py0; ++pz0;
        ++status0;
        
    } // next atom species 0
    
    // save number of molecules into pdata
    pdata->num[mol_index] = molecules;

    MUTEX_LOCK(::global_mutex);
    if ( error ) {
        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
        std::cout << info << "error " << error << " !" << std::endl;
        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
    }
    else if (failed[0] || failed[1] || failed[2])
        std::cout << info << molecules << " found (failed " << failed[0] << " " << failed[1] << " " << failed[2] << ")" << std::endl;
    else
        std::cout << info << molecules << " found" << std::endl;
    MUTEX_UNLOCK(::global_mutex);
    
    cdata->count[0] = failed[0];
    cdata->count[1] = failed[1];
    cdata->count[2] = failed[2];
    cdata->total    = molecules;
    cdata->E[E_RED] = ReductionOfEnergy;
    cdata->E[E_POT] = PotentialEnergyOfMolecules;
    cdata->E[E_KIN] = KineticEnergyOfMolecules;    
    
    return error;
}

// molecule conversion efficiency = number of pairs / number of possible pairs
// returns N0 = min(N0, N1) for two species, or N0/2 for single species (N1=0)
inline REAL_TYPE mol_efficiency(uint64_t Nmol, uint64_t *N0, uint64_t N1) {
    *N0 = (N1    == 0 ) ? ((*N0)>>1) :      // single species. this is floor(N0/2)
          ((*N0) <= N1) ? ((*N0))    : N1;  // two species
    return ((REAL_TYPE)Nmol)/(*N0);
}

// possible partner found: if both atoms are still free then pair them.
// returned values: 
//  STATUS_BOTH_FREE   = both atoms were free, i.e. molecule was created
//  STATUS_A0_PAIRED   = atom 0 already paired
//  STATUS_A1_PAIRED   = atom 1 already paired
//  STATUS_BOTH_PAIRED = both atoms already paired
//  < 0                = error
// returns molecule in mol_entry when return value = 0. input mol_entry must be nullptr.
#define STATUS_BOTH_FREE        0
#define STATUS_FIRST_PAIRED     1
#define STATUS_SECOND_PAIRED    2
#define STATUS_BOTH_PAIRED      3
inline int pair_atoms(
          struct thread_data *p0,
          struct thread_data *p1,
          REAL_TYPE *x0,
          REAL_TYPE *y0,
          REAL_TYPE *z0,
          REAL_TYPE *px0,
          REAL_TYPE *py0,
          REAL_TYPE *pz0,
          REAL_TYPE *x1,
          REAL_TYPE *y1,
          REAL_TYPE *z1,
          REAL_TYPE *px1,
          REAL_TYPE *py1,
          REAL_TYPE *pz1,
          REAL_TYPE m0,
          REAL_TYPE m1,
          uint8_t *st0, 
          uint8_t *st1,
          class species_entry *&mol_entry,
          uint64_t index0,
          uint64_t index1,
          REAL_TYPE distance
          //uint8_t checked_code
         ) {
    int error = -1;

#ifdef _DEBUG
    if (mol_entry != nullptr) return -700;
    // even for same species same atom should not be found! this is a bug!
    if (abs(distance) <= 1e-12) return -701;
    else if (( x0 ==  x1) || ( y0 ==  y1) || ( z0 ==  z1) || 
             (px0 == px1) || (py0 == py1) || (pz0 == pz1) || (st0 == st1)) return -702;
#endif        

    if (p0 == p1) {
        // both cells are owned by same thread
        // pair atoms 0 and 1 when both are free
        MUTEX_LOCK(p0->mutex);
        uint8_t status0 = *st0;
        uint8_t status1 = *st1;
        if ((IS_STATUS_FREE(status0)) && (IS_STATUS_FREE(status1))) {
            SET_STATUS_PAIRED(*st0);
            SET_STATUS_PAIRED(*st1);
            //a0->distance = a1->distance = distance;
        }
        MUTEX_UNLOCK(p0->mutex);
        error = (IS_STATUS_FREE(status0)) ? ((IS_STATUS_FREE(status1)) ? STATUS_BOTH_FREE    : STATUS_SECOND_PAIRED) :
                                            ((IS_STATUS_FREE(status1)) ? STATUS_FIRST_PAIRED : STATUS_BOTH_PAIRED  );
    }
    else {
        // different threads own cells: 
        // possible deadlock when waiting for one lock with holding another lock!
        // to avoid this we lock only one by one. 
        // in rare cases this might cause that another thread sees a0 paired although we set it unpaired afterwards.
        // in the worst case no other pair is found and we lost one molecule, 
        // but it should be a very rare case and this case is counted with failed0-2, so can be tracked easily. 
        // pair atom 0 when free
        MUTEX_LOCK(p0->mutex);
        uint8_t status0 = *st0;
        if (IS_STATUS_FREE(status0)) {
            SET_STATUS_PAIRED(*st0);
        }
        MUTEX_UNLOCK(p0->mutex);
        if (IS_STATUS_FREE(status0)) {
            // pair atom 1 when both free
            MUTEX_LOCK(p1->mutex);
            uint8_t status1 = *st1;
            if (IS_STATUS_FREE(status1)) {
                SET_STATUS_PAIRED(*st1);
            }
            MUTEX_UNLOCK(p1->mutex);
            if (IS_STATUS_FREE(status1)) {
                // success
                error = STATUS_BOTH_FREE;
            }
            else {
                // free atom 0 when atom 1 is not free
                MUTEX_LOCK(p0->mutex);
                status0 = *st0;
                if (IS_STATUS_PAIRED(status0)) {
                    SET_STATUS_FREE(*st0);
                }
                MUTEX_UNLOCK(p0->mutex);
                if (IS_STATUS_PAIRED(status0)) {
                    error = STATUS_SECOND_PAIRED;
                }
                else {
                    error = -703; // not expected to happen!
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << "error: " +p0->id << " status 0x" << std::hex << +status0 << std::endl; //", checked code 0x" << +checked_code << std::dec << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                }
            }
        }
        else {
            // note: no lock needed here.
            error = (IS_STATUS_FREE(*st1)) ? STATUS_FIRST_PAIRED : STATUS_BOTH_PAIRED;
        }
    }

    if (error == STATUS_BOTH_FREE) {
        // both atoms free: pairing successful
        
        //Statistics on molecules
        //std::cout << "pair found " << atom0 << " with atom " << atom1 << ": psd = " << PhaseSpaceDistance << std::endl;
        //mol                      = new class species_data;
        //class species_entry *emol = new class species_entry(mol, NULL, 0);
        mol_entry               = new class species_entry;
        //class species_data *mol = mol_entry->get_species(0);
        class species_data *mol = &mol_entry->species;
        
        //Center of mass of molecule: 
        REAL_TYPE Mmol = m0 + m1;
        mol->x[0] = (m0*(*x0)+m1*((*x1)+::atom_cx))/Mmol; // um
        mol->x[1] = (m0*(*y0)+m1*((*y1)+::atom_cy))/Mmol;
        mol->x[2] = (m0*(*z0)+m1*((*z1)+::atom_cz))/Mmol;
        
        // atom potential and kinetic energy
        REAL_TYPE Epot0, Epot1, Ekin0, Ekin1;
        REAL_TYPE hx, hy, hz;
        hx = ::atom_wx[0]*(*x0); // rad/s*mu
        hy = ::atom_wy[0]*(*y0);
        hz = ::atom_wz[0]*(*z0);
        Epot0 = 0.5*m0*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
        hx = ::atom_wx[1]*(*x1); // rad/s*mu
        hy = ::atom_wy[1]*(*y1);
        hz = ::atom_wz[1]*(*z1);
        Epot1 = 0.5*m1*( hx*hx + hy*hy + hz*hz ) * k_Epot; // nK
        hx = (*px0);
        hy = (*py0);
        hz = (*pz0);
        Ekin0 = 0.5*m0*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
        hx = (*px1);
        hy = (*py1);
        hz = (*pz1);
        Ekin1 = 0.5*m1*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK

        //Momentum conservation gives momentum of molecule after association
        hx = mol->x[3] = ((*px0)*m0+(*px1)*m1)/Mmol; // mm/s
        hy = mol->x[4] = ((*py0)*m0+(*py1)*m1)/Mmol;
        hz = mol->x[5] = ((*pz0)*m0+(*pz1)*m1)/Mmol;

        mol->E[E_POT] = Epot0 + Epot1;
        mol->E[E_KIN] = 0.5*Mmol*( hx*hx + hy*hy + hz*hz ) * k_Ekin; // nK
        mol->E[E_RED] = Ekin0 + Ekin1 - mol->E[E_KIN];
        
        // save atom index
        mol->atom_index[0] = index0;
        mol->atom_index[1] = index1;
        
        // save distance^2
        mol->distance = distance;
    }
    
    return error;
}

// search molecules in atoms in a linear fashion without taking into account cells.
// stops the search when it finds the first matching atom 1 for atom 0 and does not look for nearest neighbors.
// gives in average the same result for any number of threads, but different solutions might be found from shot to shot,
// even for the same seed value.
// for large atom number can take very long time!
// - this creates molecules only in one cell = cells.get_first(id)
int CreateMolecules_linear(struct thread_data *pdata, struct cmd_data *cdata) {
    int error = 0;
    uint32_t t_start = GetTickCount();
    REAL_TYPE hx, hy, hz;    
    //uint64_t num_lists = pdata->num_threads;
    //uint64_t num_atoms0 = pdata->num[0];
    uint64_t mol_num = 0;
    // count how often a possible pair is already paired by another thread. rarely happens.
    uint64_t failed0 = 0, failed1 = 0, failed2 = 0;
    // thread id
    uint8_t id = pdata->id;
    // index into cells for both species. for single species sp0 == sp1.
    uint8_t sp0 = cdata->in_index[0];
    uint8_t sp1 = cdata->in_index[1];
    uint8_t mol_index = cdata->out_index[0]; // cell index of molecule. will be allocated here.
    // scaled gamma^2 for given distance threshold
    REAL_TYPE gamma = ::cells[sp0]->get_gamma();
    
    //REAL_TYPE EpotMol, EkinMol, EredMol;
        
    REAL_TYPE PotentialEnergyOfMolecules = 0.0;
    REAL_TYPE KineticEnergyOfMolecules   = 0.0;
    REAL_TYPE ReductionOfEnergy          = 0.0;

    struct thread_data *pnext;
    bool                not_found;

    // status indicating atom 0 is unpaired. this must be different for each combination of species
    //const uint8_t checked_code = MAKE_CHECKED_CODE(cdata->out_index[1]);
    
    // list of molecules. will be inserted into ::cells[cdata->cells_index] thread-owned cell.
    class double_linked_list<class species_entry, uint64_t> molecules;

    std::string info = "thread " + std::to_string(+id) + " molecule search linear ";
    info += (find_nearest) ? "(nearest): " : "(stop at first): ";
    
    //printf("CreateMol thread %li: %li + %li atoms...\n", +id, cdata->cell->species[0].get_num(), cdata->cell->species[1].get_num());
    //uint64_t num_atoms1 = 0, num_cells1 = 0;
    uint64_t num0, max0, num1, max1, tot_num1 = 0;
    uint64_t num_threads = pdata->num_threads;
    uint64_t num_lists;
    
    REAL_TYPE *x0    = pdata->x     [sp0][0];
    REAL_TYPE *y0    = pdata->x     [sp0][1];
    REAL_TYPE *z0    = pdata->x     [sp0][2];
    REAL_TYPE *px0   = pdata->x     [sp0][3];
    REAL_TYPE *py0   = pdata->x     [sp0][4];
    REAL_TYPE *pz0   = pdata->x     [sp0][5];
    REAL_TYPE m0     = pdata->mass  [sp0];
    uint8_t *status0  = pdata->status[sp0];
    max0              = pdata->num[sp0];

    //MUTEX_LOCK(::global_mutex);
    //std::cout << info << "searching " << max0 << " atoms 0, " << ((sp0 == sp1) ? "1" : "2") << " species, gamma " << gamma << ", mol index " << +mol_index << std::endl;
    //MUTEX_UNLOCK(::global_mutex);

    if (m0 <= 0.0) {
        error = -720;
    }

    // go through all atoms0 which this thread owns
    for (num0 = 0; (!error) && (num0 < max0); ++num0) {

        // for single species (sp0 == sp1) we need to check if atom0 is already paired
        // note: we check this without locking. this is safe but status might be changing afterwards.
        //       therefore, we have to check this later again in pair_atoms within lock.
        if (IS_STATUS_FREE(*status0)) {
        
            struct thread_data *nearest_owner = nullptr;
            REAL_TYPE nearest_distance = std::numeric_limits<REAL_TYPE>::max();
            REAL_TYPE nearest_mass = 0.0;
            uint64_t nearest_index = 0;
            not_found = true;

            pnext = pdata;
            // go through all thread lists
            // notes:
            // - pnext wraps around at last thread. we count num_lists to know when to stop.
            // - for single species this code tests all atoms 2x! 
            //   at least for thread-owned atoms we could start from next atom but the difference is possibly not big.
            for (num_lists = 0; not_found && (num_lists < num_threads); pnext = pnext->next, ++num_lists) {
                // go through all atoms entries of species 1
                // for same species we start from next atom, otherwise we go through all atoms.
                REAL_TYPE *x1    = pnext->x     [sp1][0];
                REAL_TYPE *y1    = pnext->x     [sp1][1];
                REAL_TYPE *z1    = pnext->x     [sp1][2];
                REAL_TYPE *px1   = pnext->x     [sp1][3];
                REAL_TYPE *py1   = pnext->x     [sp1][4];
                REAL_TYPE *pz1   = pnext->x     [sp1][5];
                REAL_TYPE  m1    = pnext->mass  [sp1];
                uint8_t *status1 = pnext->status[sp1];
                max1             = pnext->num   [sp1];
                
                /*if (single_species && (pnext == pdata)) {
                    e1 = e0 + 1;
                    num_atoms1 = num_atoms0;
                }
                else {
                    e1 = pnext->species[sp1];
                    num_atoms1 = 0; 
                }*/
                for (num1 = 0; num1 < max1; ++num1) {

                    // take free atoms unless its the same (if same species)
                    if (IS_STATUS_FREE(*status1) && (status0 != status1)) {
                        // calculate phase space distance
                        // note: for improved efficiency we use squared psd to avoid sqrt. result is the same.

                        // TODO: calculate kinetic and potential energy per atom

                        REAL_TYPE distance;
                        if (mol_distance_index == DISTANCE_X) {
                            // real-space distance
                            hx = (*x0) - (*x1) - ::atom_cx; // um
                            hy = (*y0) - (*y1) - ::atom_cy;
                            hz = (*z0) - (*z1) - ::atom_cz;
                            distance = hx * hx + hy * hy + hz * hz; // um^2
                        }
                        else if (mol_distance_index == DISTANCE_V) {
                            // velocity distance
                            hx = (*px0) - (*px1); // mm/s
                            hy = (*py0) - (*py1);
                            hz = (*pz0) - (*pz1);
                            distance = hx * hx + hy * hy + hz * hz; // um^2
                        }
                        else if (mol_distance_index == DISTANCE_PV) {
                            // velocity-space distance
                            //Want to substract velocities in center of mass frame
                            //This is equivalent to simply substracting the velocities
                            hx = (*x0) - (*x1) - ::atom_cx; // um
                            hy = (*y0) - (*y1) - ::atom_cy;
                            hz = (*z0) - (*z1) - ::atom_cz;
                            REAL_TYPE RealSpaceDistance = hx * hx + hy * hy + hz * hz; // um^2
                            hx = (*px0) - (*px1); // mm/s
                            hy = (*py0) - (*py1);
                            hz = (*pz0) - (*pz1);
                            //To convert relative velocity to momentum in center of mass frame, use reduced mass
                            REAL_TYPE MomentumSpaceDistance = hx * hx + hy * hy + hz * hz; // (mm/s)^2
                            // [dx*dp] = kg*m^2/s = Js
                            distance = RealSpaceDistance * MomentumSpaceDistance; // um^2*(mm/s)^2 = (m^2/s*1e-9)^2
                        }
                        else if (mol_distance_index == DISTANCE_PP) {
                            // lab frame momentum-space distance
                            hx = (*x0) - (*x1) - ::atom_cx; // um
                            hy = (*y0) - (*y1) - ::atom_cy;
                            hz = (*z0) - (*z1) - ::atom_cz;
                            REAL_TYPE RealSpaceDistance = hx * hx + hy * hy + hz * hz; // um^2
                            hx = (*px0) * m0 - (*px1) * m1; // mm/s * amu
                            hy = (*py0) * m0 - (*py1) * m1;
                            hz = (*pz0) * m0 - (*pz1) * m1;
                            REAL_TYPE MomentumSpaceDistance = hx * hx + hy * hy + hz * hz; // (mm/s * amu)^2
                            // [dx*dp] = kg*m^2/s = Js
                            distance = RealSpaceDistance * MomentumSpaceDistance; // um^2*(mm/s)^2 = (m^2/s*1e-9)^2
                        }
                        else if (mol_distance_index == DISTANCE_CROSS) {
                            // cross momentum-space distance: |delta x  x  delta p|
                            // (mm/s * amu) * (mu)
                            hx =  ((*y0) - (*y1) - ::atom_cy)*((*pz0)*m0 - (*pz1)*m1) - 
                                  ((*z0) - (*z1) - ::atom_cz)*((*py0)*m0 - (*py1)*m1);
                            hy = -((*x0) - (*x1) - ::atom_cx)*((*pz0)*m0 - (*pz1)*m1) + 
                                  ((*z0) - (*z1) - ::atom_cz)*((*px0)*m0 - (*px1)*m1);
                            hz =  ((*x0) - (*x1) - ::atom_cx)*((*py0)*m0 - (*py1)*m1) - 
                                  ((*y0) - (*y1) - ::atom_cy)*((*px0)*m0 - (*px1)*m1);
                            distance = hx*hx + hy*hy + hz*hz; // amu^2*um^2*(mm/s)^2 = (amu*m^2/s*1e-9)^2
                        }
                        else if (mol_distance_index == DISTANCE_MEAN) {
                            // mean momentum-space distance: <|delta x| * |delta p|>
                            // (mm/s * amu) * (mu)
                            hx = ((*x0) - (*x1) - ::atom_cx)*((*px0)*m0 - (*px1)*m1);
                            hy = ((*y0) - (*y1) - ::atom_cy)*((*py0)*m0 - (*py1)*m1);
                            hz = ((*z0) - (*z1) - ::atom_cz)*((*pz0)*m0 - (*pz1)*m1);
                            distance = abs(hx*hy*hz); // amu^3*um^3*(mm/s)^3 = (amu*m^2/s*1e-9)^3
                        }
                        else if (mol_distance_index == DISTANCE_MAX) {
                            // maximum momentum-space distance: max|delta xi * delta pi|
                            // (mm/s * amu) * (mu)
                            hx = ((*x0) - (*x1) - ::atom_cx)*((*px0)*m0 - (*px1)*m1);
                            hy = ((*y0) - (*y1) - ::atom_cy)*((*py0)*m0 - (*py1)*m1);
                            hz = ((*z0) - (*z1) - ::atom_cz)*((*pz0)*m0 - (*pz1)*m1);
                            // amu^2*um^2*(mm/s)^2 = (amu*m^2/s*1e-9)^2
                            hx *= hx;
                            hy *= hy;
                            hz *= hz;
                            if (hx >= hy) distance = (hx >= hz) ? hx : hz;
                            else          distance = (hy >= hz) ? hy : hz;  
                        }
                        else {
                            // unknown distance measure
                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                            std::cout << "CreateMolecules_linear error: distance measure ' " << distance_all[mol_distance_index] << " ' (" << mol_distance_index << " not implemented!" << std::endl;
                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                            return -1;
                        }

                        if (distance < gamma) {
                            if (find_nearest) {
                                if (distance < nearest_distance) {
                                    // possible pair found: find nearest distance of pairs                            
                                    nearest_distance = distance;
                                    nearest_index    = num1;
                                    nearest_owner    = pnext;
                                    nearest_mass     = m1;
                                }
                            }
                            else {
                                // possible pair found: try to pair   
                                //class species_data* mol = NULL;
                                class species_entry* mol_entry = NULL;
                                //error = pair_atoms(c0, c1, a0, a1, mol_entry, distance, checked_code);
                                error = pair_atoms(pdata, pnext,
                                                   x0, y0, z0, px0, py0, pz0, 
                                                   x1, y1, z1, px1, py1, pz1,
                                                   m0, m1,
                                                   status0, status1,
                                                   mol_entry, 
                                                   num0, num1, 
                                                   distance);
                                if (error == STATUS_BOTH_FREE) {
                                    // pairing successful: save into list of molecules
                                    error = 0;
                                    //molecules.add(new species_entry(mol, NULL, 0));
                                    molecules.add(mol_entry);
                                    ++mol_num;
                                    //class species_data *mol = mol_entry->get_species(0);
                                    class species_data *mol = &mol_entry->species;
                                    ReductionOfEnergy          += mol->E[E_RED]; // nK
                                    PotentialEnergyOfMolecules += mol->E[E_POT]; // nK
                                    KineticEnergyOfMolecules   += mol->E[E_KIN]; // nK
                                    not_found = false;  // stop cell 1 search
                                    break;              // stop atom 1 search
                                }
                                else if (error == STATUS_FIRST_PAIRED) {
                                    // atom 0 paired by another thread: stop searching atom 0 & 1
                                    error = 0;
                                    ++failed0;
                                    not_found = false;  // stop cell 1 search
                                    break;              // stop atom 1 search
                                }
                                else if (error == STATUS_SECOND_PAIRED) {
                                    // atom 1 paired by another thread: continue searching atom 1
                                    // never seen. unprobable since we already check status1 for atom1.
                                    error = 0;
                                    ++failed1;
                                }
                                else if (error == STATUS_BOTH_PAIRED) {
                                    // both atoms paired by another thread: stop searching atom 0 & 1
                                    error = 0;
                                    ++failed2;
                                    not_found = false;  // stop cell 1 search
                                    break;              // stop atom 1 search
                                }
                                else {
                                    // error!
                                    not_found = false;  // stop cell 1 search
                                    break;              // stop atom 1 search
                                }
                            }
                        } // distance < gamma
                    } // status1 free

                    // increment atom 1 pointers
                    ++x1; ++y1; ++z1;
                    ++px1; ++py1; ++pz1;
                    ++status1;

                } // next atom1
                tot_num1 += max1;
                
            } // next thread list

            if (find_nearest) {
                if (nearest_owner != nullptr) {
                    // possible pair found
                    //class species_data* mol = NULL;
                    class species_entry* mol_entry = NULL;
                    //error = pair_atoms(c0, nearest_cell, a0, nearest_atom, mol_entry, nearest_distance, checked_code);
                    error = pair_atoms(pdata, nearest_owner, 
                                        x0, y0, z0, px0, py0, pz0,
                                        nearest_owner->x[sp1][0] + nearest_index, 
                                        nearest_owner->x[sp1][1] + nearest_index, 
                                        nearest_owner->x[sp1][2] + nearest_index, 
                                        nearest_owner->x[sp1][3] + nearest_index,
                                        nearest_owner->x[sp1][4] + nearest_index,
                                        nearest_owner->x[sp1][5] + nearest_index, 
                                        //Epot0, 
                                        //nearest_owner->E[sp1][E_POT] + nearest_atom, 
                                        //Ekin0, 
                                        //nearest_owner->E[sp1][E_KIN] + nearest_atom,
                                        m0, nearest_mass, 
                                        status0, 
                                        nearest_owner->status[sp1] + nearest_index, 
                                        mol_entry, 
                                        num0, nearest_index,
                                        nearest_distance);
                    if (error == STATUS_BOTH_FREE) {
                        // pairing successful: save into list of molecules. hash_key = 0.
                        //molecules.add(new species_entry(mol, NULL, 0));
                        molecules.add(mol_entry);
                        not_found = false; // do not mark atom 0 as unpaired.
                        ++mol_num;
                        //class species_data *mol = mol_entry->get_species(0);
                        class species_data *mol = &mol_entry->species;
                        ReductionOfEnergy          += mol->E[E_RED]; // nK
                        PotentialEnergyOfMolecules += mol->E[E_POT]; // nK
                        KineticEnergyOfMolecules   += mol->E[E_KIN]; // nK
                    }
                    else if (error == STATUS_FIRST_PAIRED) {
                        // atom 0 is paired by another thread: but it might be different atoms!
                        // maybe I have seen this case? but much less frequent than error == 2.
                        error = 0;
                        ++failed0;
                        not_found = false; // do not mark atom 0 as unpaired
                    }
                    else if (error == STATUS_SECOND_PAIRED) {
                        // atom 1 is paired by another thread: but it might be different atoms!
                        // this happens rarely, but happens.
                        error = 0;
                        ++failed1;
                        // not_found = true already set
                    }
                    else if (error == STATUS_BOTH_PAIRED) {
                        // atom 0 and 1 are paired by another thread: could be even different atoms!
                        // this happens extremely rarely, but happens.
                        error = 0;
                        ++failed2;
                        not_found = false; // do not mark atom 0 as unpaired
                    }
                    else {
                        // error!
                        break;
                    }
                } // end pair found
            } // find_nearest
            
            if (not_found) {
                // mark atom 0 if no pair found
                uint8_t old_status;
                MUTEX_LOCK(pdata->mutex);
                old_status = *status0;
                if (IS_STATUS_FREE(old_status)) SET_STATUS_UNPAIRED(*status0);
                MUTEX_UNLOCK(pdata->mutex);
                if (!IS_STATUS_FREE(old_status)) ++failed2; // can happen
            }
        } // atom 0 is free

        // increment atom 0 pointers
        ++x0; ++y0; ++z0;
        ++px0; ++py0; ++pz0;
        ++status0;
                
    } // next atom 0
                    
    if (!error) {
        // insert molecules into first thread-owned cell.
        // we do not need to lock cells
        class cell_data *cell = ::cells[mol_index]->get_first();
        while(cell) {
            if (cell->get_owner_id() == id) {
                class species_entry *first, *last;
                last  = molecules.get_last();
                first = molecules.remove_all();
                cell->species.add(first, last, mol_num);
                break;
            }
            cell = cell->get_next();
        }
        if (cell == NULL) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << "CreateMolecules_linear error: thread owned cell not found?" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            error = -2;
        }
    }
        
    if (error) {
        // on error delete already created molecules
        class species_entry *entry = molecules.remove_all();
        while (entry) {
            class species_entry *next = entry->get_next();
            delete entry;
            entry = next;
        }
    }

    MUTEX_LOCK(::global_mutex);
    if (error) {
        std::cout << info << "error " << error << " !" << std::endl;
    }
    else {
        std::cout << info << mol_num << " molecules found, " << max0 << " x " << std::setprecision(1) << ((max0==0) ? 0.0 : (REAL_TYPE)tot_num1/max0) << " atoms, " << 
#if NUM_CELLS > 1
            num_cells0 << " x " << ((max0==0) ? 0.0 : (REAL_TYPE)num_cells1/max0) << " cells (" 
#else   
            num_threads << " lists ("
#endif 
            << (((REAL_TYPE)get_ticks_delta(t_start)) / 1e3) << " s)";
        if (failed0 || failed1 || failed2)
            std::cout << ", failed {" << failed0 << ", " << failed1 << ", " << failed2 << "} (ok)" << std::endl;
        else
            std::cout << std::endl;
    }
    MUTEX_UNLOCK(::global_mutex);

    cdata->count[0] = failed0;
    cdata->count[1] = failed1;
    cdata->count[2] = failed2;
    cdata->total    = mol_num;
    cdata->E[E_RED] = ReductionOfEnergy;
    cdata->E[E_POT] = PotentialEnergyOfMolecules;
    cdata->E[E_KIN] = KineticEnergyOfMolecules;
                 
    return error;
}

/* Andi: not used so far.
// find nearest neighbors between two species using distance measure and gamma.
// the resulting neighbors are saved into shared list 'nearest' sorted by increasing distance.
// only neighbors are returned where distance <= gamma for the given distance index.
// function can be called for all threads in parallel.
// species can be the same or different with different mass.
// walks through cells_x or cells_v to efficiently find nearest neighbors.
// note: when species[0] == species[1] and neither DISTANCE_X nor DISTANCE_V 
//       the function will find each pair twice and do twice the required work!
//       in this case the function works as for two different species. 
//       the problem is that searching atom1 > atom0 (in terms of atom index) is difficult 
//       when working with two different unsorted lists cells_x and cells_v.
//       one could sort atom list in each cell by increasing atom index, 
//       when searching atoms in cell, go through atom list in reverse order,
//       and when atom1 index <= atom0 index stop search the cell.
int FindNearest(const struct thread_data *pdata, struct cmd_data *cdata) {
    int error = 0;
    
    class cell_data    *cell0 , *cell1;
    class species_entry *entry0, *entry1;

    uint8_t i0 = cdata->in_index[0];
    uint8_t i1 = cdata->in_index[1];
    bool same_species = (i0 == i1); 
    //uint8_t distance_index = cdata->distance_index;
    //REAL_TYPE  gamma = (cdata->gamma == 0) ? DBL_MAX : cdata->gamma;
    uint8_t distance_index = ::cells[i0]->get_distance_index();
    REAL_TYPE gamma           = ::cells[i0]->get_gamma();

    cell0 = ::cells[i0]->get_first();
    for (; cell0 != NULL; cell0 = cell0->get_next()) {
        // get all atoms in cell0
        entry0 = cell0->species.get_first();
        for (; entry0 != NULL; entry0 = entry0->get_next()) {
            //atom0 = entry0->get_species(0);
            
            // find all cells and atoms where distance measure can be fulfilled.
            //cell1 = entry0->get_next_distance(cells1, gamma_squared);
            cell1 = (same_species) ? cell0->get_next() : ::cells[i1]->get_first();
            while (cell1) {
                entry1 = cell1->species.get_first();
                while (entry1) {

                    if (entry1 != entry0) {
                        // save pair of atoms with hash_key = distance^2
                        // TODO: save into out_index[0]
                        uint64_t d = entry0->get_hash_key() - entry1->get_hash_key();
                        ::nearest.add_sorted(new class species_entry(entry0->get_species(0), 
                                                                    entry1->get_species(0), 
                                                                    d*d));
                    }

                    // next entry within cell1
                    entry1 = entry1->get_next();
                }
                // next cell1
                cell1 = cell1->get_next();
            }
        } // next entry0
    } // next cell0 

    return error;
}
*/


/* gray codes for seeding with thread_id
static const uint8_t gray[THREADS_MAX_NUM] = {
    0x00, 0x01, 0x03, 0x02, 0x06, 0x07, 0x05, 0x04,  
    0x0c, 0x0d, 0x0f, 0x0e, 0x0a, 0x0b, 0x09, 0x08,
    0x18, 0x19, 0x1b, 0x1a, 0x1e, 0xff, 0x1d, 0x1c,
    0x14, 0x15, 0x17, 0x16, 0x12, 0x13, 0x11, 0x10,
    0x30, 0x31, 0x33, 0x32, 0x36, 0x37, 0x35, 0x34,
    0x3c, 0x3d, 0x3f, 0x3e, 0x3a, 0x3b, 0x39, 0x38,
    0x28, 0x29, 0x2b, 0x2a, 0x2e, 0x2f, 0x2d, 0x2c,
    0x24, 0x25, 0x27, 0x26, 0x22, 0x23, 0x21, 0x20
};

// mixing codes for seeding of WELL1024. this assumes WELL1024a_R = 32.
static const uint32_t mix[8] = {
    0x00000000, 
    0x0f0f0f0f, 
    0xf0f0f0f0,
    0x00ff00ff, 
    0xff00ff00, 
    0xff0000ff, 
    0x00ffff00, 
    0xffffffff
};*/

// create or reset seed generator to content of seed_use string or random device.
// if seed == nullptr uses ::seed_use or std::random_device when seed_use is empty to get initial seed values.
// if seed != nullptr uses seed_values to reset seed to initial values.
// this allows to call function several times (with returned seed_values) to reset seed generator.
// returns 0 if ok, otherwise error code.
// notes:
// - std::seed_seq creates the same sequence of numbers for each call to generate!
//   this might depend on OS and might change in the future,
//   but it makes it complicate to init all threads with different seed values.
// - std::seed_seq cannot be reset to initial value, which would be convenient to reset RNGs
//   this way regardless of number of random numbers generated in tests before calculation
//   after reset the actual calculation is always starting with the same random numbers.
//   this allows to reproduce old results even if code/tests were changed in the meantime.
// - to overcome these limits we implement here the Mersenne Twister random number generator
//   which will generate in each call new random numbers, and can be reset at any time.
int get_seed(random_generator<uint32_t> *& seed, std::vector<uint32_t> & seed_values) {
    int error = 0;
    if (seed == nullptr) {
        // initial call
        if (::seed_use.size() == 0) {
            std::cerr << "info: using " << SEED_VALUES_NUM << " random hardware seed values" << std::endl;        
            // non-deterministic random number generator using hardware entropy source (ideally).
            // this is a very precious resource and should not be used excessive. we call it only once!
            std::random_device rd;
            for(uint8_t i = 0; i < SEED_VALUES_NUM; ++i) {
                seed_values.insert(seed_values.end(), rd());
            }
        }
        else {
            std::cerr << "info: reading seed values: " << seed_use << " ... " << std::endl;
            error = to_vector<uint32_t>(
                seed_use, 
                seed_values, 
                std::string("unsigned integer, decimal or hex (with " NUM_HEX ") needed!"), 
                true);
            if (error || (seed_values.size() == 0)) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << std::endl << "error: reading seeding values (" << error << ")!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -141;
            }
            else if ((seed_values.size() == 1) && (seed_values[0] == 0)) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << std::endl << "error: seed values must not be zero!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -142;
            }

        }
        if(!error) {
            // create seed generator seeded with seed_values
            std::seed_seq seq(seed_values.begin(), seed_values.end());
            seed = new RNG_uni_direct_int<uint32_t, std::mt19937_64>(seq);
            if (seed == nullptr) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: creating random number generator!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -143;
            }
            else {
                std::cerr << "info: got " << std::setw(3) << seed_values.size() << " seed values: ";
            }
        }
    }
    else if (seed_values.size() == 0) {
        SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
        std::cerr << std::endl << "error: reset seed generator with empty seed_values!" << std::endl;
        RESET_CONSOLE_COLOR(STDERR_HANDLE);
        error = -144;
    }
    else {
        // reset seed generator to previous seed_values
        std::cerr << "info: reset " << seed_values.size() << " seed values: ";
        seed->seed(seed_values);
    }
    
    if(!error) {
        // output seed values
        std::cerr << V_OPEN << std::hex;
        for(uint8_t i = 0; i < seed_values.size(); ++i) {
            if (i != 0) std::cerr << V_SEP;
            std::cerr << "0x" << seed_values[i];
        }
        std::cerr << V_CLOSE << std::dec << std::endl << std::endl;
    }
    
    return error;
}

// generates the given random number generator and seeds it for given helper thread id
// for any non NULL pointers returns uniform, normal or basic random number generator
// returned generators must be deleted after use
// thread id is used to initialize each of them individually even if seed is the same
// returns 0 if ok, otherwise error
template <typename real_type, typename int_type>
int get_RNG(uint8_t rnd_index, 
            random_generator<uint32_t> & seed, 
            random_generator<real_type> **uniform, 
            random_generator<real_type> **normal,
            random_generator<int_type > **basic
            ) { 
    switch ( rnd_index ) {
        // random number generators with at least 64bit intenal state can be always used
        case RND_MERSENNE_TWISTER: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        std::mt19937_64
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        std::mt19937_64,             
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type, 
                                                        std::mt19937_64 
                                                       >(seq);
                }
            }
            break;
        case RND_LEHMER64: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        Lehmer64
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        Lehmer64,
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type,
                                                        Lehmer64
                                                       >(seq);
                }
            }
            break;
        case RND_LEHMER128: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        Lehmer128
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                                        std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        Lehmer128, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type, 
                                                        Lehmer128
                                                       >(seq);
                }
            }
            break;
        case RND_WYHASH64: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type,
                                                        int_type, 
                                                        Wyhash64
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        Wyhash64, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(4);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type,
                                                        Wyhash64
                                                       >(seq);
                }
            }
            break;
        //case RND_PCG64: should be good but not implemented.
        //    break;
        //case RND_XOSHIRO256: modern but not implemented. not all of this type are good ones however!
        //    break;
        //case RND_SPLIMIX64: # not implemented.
        //    break;
#if REAL_PRECISION == REAL_SINGLE
        // DO NOT USE! except for testing!
        // random number generators with only 32bit internal state can be used only with float!
        // these have way too small period and state space and should never be used!
        // all of these fail miserably and very fast PractRand test!!!
        case RND_LEHMER32: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        Lehmer32
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        Lehmer32, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type, 
                                                        Lehmer32
                                                       >(seq);
                }
            }
            break;
        case RND_WELL1024: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(WELL1024a_R);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        WELL1024
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(WELL1024a_R);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        WELL1024, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(WELL1024a_R);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type, 
                                                        WELL1024
                                                       >(seq);
                }
            }
            break;
        case RND_LCG_MINST: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        std::minstd_rand
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type,
                                                        std::minstd_rand,
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type,
                                                        std::minstd_rand
                                                       >(seq);
                }
            }
            break;   
        case RND_RANLUX24: {
                if (uniform)
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type,
                                                        int_type, 
                                                        std::ranlux24
                                                       >(seq);
                }
                if (normal) { 
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        std::ranlux24, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(1);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type,
                                                        std::ranlux24
                                                       >(seq);
                }
            }
            break;
        case RND_RANLUX48: {
                if (uniform) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *uniform  = new RNG_uni_direct_real<real_type, 
                                                        int_type, 
                                                        std::ranlux48
                                                       >(seq);
                }
                if (normal) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *normal   = new RNG_generic        <real_type, 
                                                        std::ranlux48, 
                                                        std::normal_distribution<real_type>
                                                       >(seq);
                }
                if (basic) {
                    std::vector<uint32_t> seed_values(2);
                    seed.generate(seed_values);
                    std::seed_seq seq(seed_values.begin(), seed_values.end());
                    *basic    = new RNG_uni_direct_int <int_type,
                                                        std::ranlux48
                                                       >(seq);
                }
            }
            break;
#endif
        default:
            std::cout << "random number generator " << rnd_use << " (index " << +rnd_index << ") unknown!? older ones are only available with single precision floating point arithmetics but are not recommended!" << std::endl;
            return -1;
    }
    return 0;
}

// check num random numbers generated by random number generator rng
// we count how often each bit changes. we expect each bit to change half of the time.
// this is intended for integer, not for double or float output!
// TODO: check each bit individually and correlation with other bits and previous bits (requires matrix)
template <typename T>
int test_RNG(random_generator<T> *rng, uint64_t num, uint64_t show, std::string info) {
    int error = 0;
    const uint8_t bits = sizeof(T)<<3;
    std::cout << info << ": counting bit changes of " << num << " random numbers (" << +bits << "bits) ..." << std::endl;

    T value, old = 0, shift;
    uint64_t* count = new uint64_t[bits];
    for (uint8_t i = 0; i < bits; ++i) count[i] = 0;
    for (uint64_t k = 0; k < num; ++k) {
        value = (*rng)();
        shift = value;
        if (k < show) std::cout << std::setw(3) << k << " : " << "0x" << std::hex << std::setw(8) << value << std::dec << std::endl;
        if (k > 0) {
            for (uint8_t i = 0; i < bits; ++i) {
                if ((shift & 1) ^ (old & 1)) ++count[i];
                shift >>= 1;
                old   >>= 1;
            }
        }
        old = value;
    }
    uint64_t sum = 0, sum_sqr = 0;
    for (uint8_t i = 0; i < bits; ++i) {
        sum     += count[i];
        sum_sqr += count[i] * count[i];
    }            
    double mean  = ((double)sum) / bits;
    double stdev = std::sqrt((((double)sum_sqr) - bits*mean * mean)/(bits-1));
    double test = 0.0;
    for (uint8_t i = 0; i < bits; ++i) {
        test += (count[i] - mean)*(count[i] - mean);
    }
    mean  /= num;
    stdev /= num;
    test   = std::sqrt(test / (bits-1)) / num;

    if ((mean < (0.5-stdev)) || (mean > (0.5+stdev)) || (stdev > 0.1)) {
        error = -1154;
    }
    else if (abs(stdev - test) > 1e-10) {
        error = -1155;
    }

    if (error) SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
    else       SET_CONSOLE_COLOR_INFO (STDOUT_HANDLE);
    std::cout << "\tprobability of bit change: " << std::fixed << std::setprecision(3) << mean << " +/- " << stdev;

    if (error) {
        if (error == -1154) 
            std::cout << " error! p = 0.5 expected!" << std::endl;
        else
            std::cout << std::fixed << std::setprecision(3) << " error check " << stdev << " != " << test << " (error = " << std::scientific << abs(stdev - test) << ")" << std::endl;
    }
    else {
        std::cout << " (ok)" << std::endl;
    }
    std::cout << std::endl;
    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
    
    delete[] count;

    return error;
}

template <typename T>
class bit_convert {
private: 
    typedef int (*callback_func)(T);
    T reg_bits;
    T bits;
    T reg;
    callback_func callback;
public:
    bit_convert(T reg_bits, callback_func callback) { 
        this->reg_bits = reg_bits;
        this->callback = callback;
        reg = 0;
        bits = 0; 
    }
    ~bit_convert() {}
    
    // add add_bits to register
    // if reg_bits reached in register call callback function and return result
    // otherwise return 0
    int add(T add_bits, T data) {
        int result = 0;
        data &= (((T)1<<add_bits)-1);
        T num = bits + add_bits;
        if (num < reg_bits) {
            bits += add_bits;
            reg = (reg << add_bits) | data; 
        }
        else {
            T rem = reg_bits - bits;
            bits = add_bits - rem;
            reg = (reg << rem) | (data >> bits);
            result = callback(reg);
            reg = data & (((T)1<<bits)-1);
        }
        return result;
    }
    
    // return number of bits
    T get_bits(void) { return bits; }
};

// is called whenever 64/32bits are ready
// should return always 1
int RNG_write_64bits(uint64_t value) {
    return (int) std::fwrite((void*) &value, sizeof(value), 1, stdout);
}
int RNG_write_32bits(uint32_t value) {
    return (int) std::fwrite((void*) &value, sizeof(value), 1, stdout);
}

// streams binary random numbers on stdout generated with the random number generator index.
// this data can be read and analyzed by PractRand.
// mode gives which random numbers should be used, see RNG_TEST_ constants.
// notes: 
// - RNG_stream returns when stdin is closed.
// - ensure that nothing is written to std::stdout before this test starts! std::cerr is ok.
// - on windows console must be set into binary mode!
#define RNG_TEST_NONE           0       // do not stream binary random numbers
#define RNG_TEST_INT64          1       // directly output random uint64_t
#define RNG_TEST_INT32          2       // directly output random uint32_t
#define RNG_TEST_DOUBLE         3       // get random double, output converted uint64_t
#define RNG_TEST_FLOAT          4       // get random float , output converted uint32_t
#define RNG_TEST_NORMAL_DOUBLE  5       // get normal-distributed double, output converted uint64_t
#define RNG_TEST_NORMAL_FLOAT   6       // get normal-distributed float, output converted uint32_t
#define STREAMING               "STREAMING on stdout"   // keyword on stderr indicating when starting to stream
int RNG_stream(uint8_t rnd_index, random_generator<uint32_t> & seed, uint8_t mode) {
    int error = 0;
    uint64_t loops = 0;
    std::string dist = "uniform distributed";

#if defined(_WIN32) || defined(_WIN64)
    // windows makes big problems!
    int result = _setmode(_fileno(stdout), _O_BINARY);
    if (result == -1) {
        std::cerr << "error: could not set stdout to binary mode! expect PractRand to fail!" << std::endl;
    }
    else {
        std::cerr << "note: stdout set to binary mode!" << std::endl;
    }
#endif

    switch (mode) {
        case RNG_TEST_INT64:
            {
                std::cerr << STREAMING << " using " << rnd_use << " (" << dist << " int64) ... " << std::endl;

                // use integer 64bit random number generator
                // this is less strict since in simulation we use random doubles generated from integers
                random_generator<uint64_t> *rng = nullptr;
                error = get_RNG<double, uint64_t>(rnd_index, seed, NULL, NULL, &rng);

                if (!error) {

                    while (1) {
                        uint64_t value = (*rng)();
#ifdef STREAM_DATA_OUT
                        if (loops < STREAM_DATA_OUT) {
                            std::cerr << "0x" << std::hex << std::setw(16) << std::setfill('0') << value << std::dec << std::endl;
                            if (loops == (STREAM_DATA_OUT-1)) std::cerr << std::endl << std::flush;
                        }
#endif
                        if (std::fwrite((void*) &value, sizeof(value), 1, stdout) != 1) break;
                        ++loops;
                    }
                    delete rng;
                    std::cerr << "RNG stream " << loops << " loops ok" << std::endl;
                }
            }
            break;
        case RNG_TEST_INT32:
            {
                std::cerr << STREAMING << " with " << rnd_use << " (" << dist << " int32) ... " << std::endl;

                // use integer 32bit random number generator
                // this is less strict since in simulation we use random floats generated from integers
                random_generator<uint32_t> *rng = nullptr;
                error = get_RNG<float, uint32_t>(rnd_index, seed, NULL, NULL, &rng);

                if (!error) {
                    while (1) {
                        uint32_t value = (*rng)();
#ifdef STREAM_DATA_OUT
                        if (loops < STREAM_DATA_OUT) {
                            std::cerr << "0x" << std::hex << std::setw(8) << std::setfill('0') << value << std::dec << std::endl;
                            if (loops == (STREAM_DATA_OUT-1)) std::cerr << std::endl << std::flush;
                        }
#endif
                        if (std::fwrite((void*) &value, sizeof(value), 1, stdout) != 1) break;
                        ++loops;
                    }
                    delete rng;
                    std::cerr << "RNG stream " << loops << " loops ok" << std::endl;
                }
            }
            break;
        case RNG_TEST_NORMAL_DOUBLE:
            dist = "normal distributed";
        case RNG_TEST_DOUBLE:
            {
                std::cerr << STREAMING << " with " << rnd_use << " (" << dist << " double -> int64) ... " << std::endl;
                // convert double to uint64_t
                // since we can use only mantissa with 52bits we have to combine several rng's for 64bits.
                // from_double returns mantissa MSB aligned.
#if ( TO_DOUBLE == MANT_EXP )
                const uint8_t reg_bits = 52;
#elif ( TO_DOUBLE == MULT )                
                const uint8_t reg_bits = 53;
#endif
                class bit_convert<uint64_t> bc(64, RNG_write_64bits);
                random_generator<double> *rng = nullptr;
                if (mode == RNG_TEST_DOUBLE)
                    error = get_RNG<double, uint64_t>(rnd_index, seed, &rng, NULL, NULL);
                else
                    error = get_RNG<double, uint64_t>(rnd_index, seed, NULL, &rng, NULL);
                if (!error) {
                    double num64 = 0.0, old = (*rng)();
                    while (1) {
                        double value = (*rng)();
                        if (mode == RNG_TEST_NORMAL_DOUBLE) {
                            // inverse Box-Muller transform. we use only one number calculated with previous value.
                            double U1 = std::exp(-(old*old+value*value)/2.0);
                            //double U2 = std::atan2(value,old)/(2*pi);
                            old = value;
                            value = U1;
                        }
                        int result = bc.add(reg_bits, from_double(value)>>(64-reg_bits));
                        if (result == 1) ++num64;
                        else if (result != 0) break;
                        //if (loops == 200) break; // for testing
                        ++loops;
                    }
                    delete rng;
                    uint64_t rem = bc.get_bits();
                    if (loops*reg_bits == (num64*64 + rem)) {
                        std::cerr << "RNG stream " << loops << " x " << reg_bits << "bits = " << num64 << " x 64bits + " << bc.get_bits() << " remaining ok" << std::endl;
                    }
                    else {
                        std::cerr << "RNG stream " << loops << " x " << reg_bits << "bits = " << num64 << " x 64bits + " << bc.get_bits() << " remaining" << std::endl;
                        error = -1000;
                    }
                }
            }
            break;
        case RNG_TEST_NORMAL_FLOAT:
            dist = "normal distributed";
        case RNG_TEST_FLOAT:
            {
                std::cerr << STREAMING << " with " << rnd_use << " (" << dist << " float -> int32) ... " << std::endl;
                // convert float to uint32_t
                // since we can use only mantissa with 23bits we have to combine several rng's for 32bits.
                // from_double returns mantissa MSB aligned.
#if ( TO_DOUBLE == MANT_EXP )
                const uint8_t reg_bits = 23;
#elif ( TO_DOUBLE == MULT )                
                const uint8_t reg_bits = 24;
#endif
                class bit_convert<uint32_t> bc(32, RNG_write_32bits);
                random_generator<float> *rng = nullptr;
                if (mode == RNG_TEST_FLOAT)
                    error = get_RNG<float, uint32_t>(rnd_index, seed, &rng, NULL, NULL);
                else
                    error = get_RNG<float, uint32_t>(rnd_index, seed, NULL, &rng, NULL);
                if (!error) {
                    float num32 = 0.0, old = (*rng)();
                    while (1) {
                        float value = (*rng)();
                        if (mode == RNG_TEST_NORMAL_FLOAT) {
                            // inverse Box-Muller transform. we use only one number calculated with previous value.
                            float U1 = std::exp(-(old*old+value*value)/2.0f);
                            //float U2 = std::atan2(value,old)/(2*pi);
                            old = value;
                            value = U1;
                        }                        
                        int result = bc.add(reg_bits, from_double(value)>>(32-reg_bits));
                        if (result == 1) ++num32;
                        else if (result != 0) break;
                        //if (loops == 200) break; // for testing
                        ++loops;
                    }
                    delete rng;
                    uint64_t rem = bc.get_bits();
                    if (loops*reg_bits == (num32*32 + rem)) {
                        std::cout << "RNG stream " << loops << " x " << reg_bits << "bits = " << num32 << " x 32bits + " << bc.get_bits() << " remaining ok" << std::endl;
                    }
                    else {
                        std::cout << "RNG stream " << loops << " x " << reg_bits << "bits = " << num32 << " x 32bits + " << bc.get_bits() << " remaining" << std::endl;
                        error = -1001;
                    }
                }
            }
            break;
        default:
            std::cerr << "RNG stream unknown mode " << mode << " !" << std::endl;
            error = -1002;
            break;
    }

    return error;
}

// test real-valued distribution and plot histogram. 
// can be uniform or normal distributed.
#define STDDEV_UNIFORM  std::sqrt(1.0/12.0)      // standard deviation of uniform distribution from 0 to 1
int test_distribution(
        random_generator<REAL_TYPE> *rng, 
        const uint64_t num, 
        bool normal, 
        const REAL_TYPE bin_add, 
        const REAL_TYPE bin_sigma, 
        const uint16_t num_bins, 
        const uint8_t draw_rows) {
    int error = 0;
    REAL_TYPE min, max, bin_min, bin_max, sum, sum_sqr;
    uint64_t below, above, count_max;

#if REAL_PRECISION == REAL_DOUBLE
    std::string real_type("double");
#else
    std::string real_type("float");
#endif
    std::string info(std::to_string(num) + " " + (normal ? "normal" : "uniform") + " distributed " + real_type);

    uint64_t *count = new uint64_t[num_bins];
    if (count == nullptr) {
        std::cout << info << ": could not allocate " << num_bins << " bins for histogram!" << std::endl;        
    }
    else {
        REAL_TYPE old = 0; 
        for (uint8_t rep = 0; rep < 2; ++rep) {
            std::cout << "test " << info << " ... " << std::endl;

            sum = 0.0;
            sum_sqr = 0.0;
            min = std::numeric_limits<REAL_TYPE>::max(); 
            max = std::numeric_limits<REAL_TYPE>::min(); 
            if (normal) {
                bin_min = -bin_sigma;
                bin_max = +bin_sigma;
            }
            else {
                bin_min = 0.0 - bin_add;
                bin_max = 1.0 + bin_add;
            }
            below = 0; 
            above = 0;
            count_max = 0;
        
            for (uint16_t i = 0; i < num_bins; ++i) count[i] = 0;

            for (uint64_t i = 0; i < num; ++i) {
                REAL_TYPE value = (*rng)();
                if (rep == 1) {
                    if (normal) {
                        // Box-Muller transform gives normal distribution from uniform distribution
                        // we use only one number
                        REAL_TYPE R     = (old > 0.0) ? std::sqrt(-2*std::log(old)) : 0.0;
                        REAL_TYPE theta = 2*Pi*value;
                        REAL_TYPE Z0 = R*std::cos(theta);
                        //REAL_TYPE Z1 = R*std::sin(theta);
                        old = value;
                        value = Z0;
                    }
                    else {
                        // inverse Box-Muller transform gives uniform distribution from normal distribution
                        // we use only one number
                        REAL_TYPE U1 = std::exp(-(old*old+value*value)/2.0);
                        //REAL_TYPE U2 = std::atan2(value,old)/(2*pi);
                        old = value;
                        value = U1;
                    }
                }
                if      (value > max) max = value;
                else if (value < min) min = value;
                sum     += value;
                sum_sqr += value*value;
                
                int64_t bin = static_cast<uint64_t>(std::floor((value - bin_min)*num_bins/(bin_max-bin_min)));
                if      (bin <  0                  ) ++below;
                else if (bin >= ((int64_t)num_bins)) ++above; 
                else {
                    ++count[bin];
                    if (count[bin] > count_max) count_max = count[bin];
                }
            }
            
            std::cout << "histogram between min = " << bin_min << " and max = " << bin_max << ": max. count " << count_max;
            if (below || above) {
                std::cout << ". warning: ";
                if (below)          std::cout << below << " counts below min!";  
                if (below && above) std::cout << " and ";
                if (above)          std::cout << above << " counts above max!";
            }
            std::cout << std::endl;
            
            // draw histogram
            for (uint8_t row = 0; row < draw_rows; ++row) {
                double f = ((double)(draw_rows-row-1))/(draw_rows-1);
                std::cout << std::setw(6) << static_cast<uint64_t>(std::floor(f*count_max)) << " ";
                for (uint16_t i = 0; i < num_bins; ++i) {
                    double g = ((double)count[i])/count_max;
                    std::cout << (((g > 0) && (g >= f)) ? '+' : '.');
                }
                std::cout << std::endl;
            }
            // axis coordinates (note: placement is not 100% correct)
            std::cout << std::setw(7) << " ";
            REAL_TYPE dx = (bin_max-bin_min)/num_bins;
            for (uint16_t i = 0; i < num_bins; ++i) {
                REAL_TYPE x = bin_min + i*dx;                           // left coordinate of bin
                int16_t pr = static_cast<int16_t>(std::ceil(x));        // next higher integer from left
                int16_t pl = static_cast<int16_t>(std::floor(x+dx));    // next lower integer from right
                if (pl == pr) {
                    std::cout << (pl % 10);
                    if (pl < 0) ++i;                                    // negative sign needs 2 places
                }
                else {
                    std::cout << ' ';
                }
            }
            std::cout << std::endl;
            //std::cout << info << " generated." << std::endl;

            //for (uint16_t i = 0; i < num_bins; ++i) std::cout << std::setw(3) << i << ": " << count[i] << std::endl;
            
            // calculate mean and standarddeviation and compare with expected values
            REAL_TYPE mean    = sum / num;
            REAL_TYPE stddev  = std::sqrt((sum_sqr - num*mean*mean)/(num-1.0));
            REAL_TYPE std_err = std::sqrt (sum_sqr + num*mean*mean)/(num-1.0) *2.0; // error propagation of stddev
            // TODO: multiplied std_err * 2.0 since sometimes was complaining!
            if (normal) {
                if ((abs(mean) >= std_err) || (abs(stddev - 1.0) >= std_err)) error = -190;
            }
            else {
                if ((abs(mean - 0.5) >= std_err) || (abs(stddev - STDDEV_UNIFORM) >= std_err)) error = -191;
            }

            // output result            
            if (error) SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            else       SET_CONSOLE_COLOR_INFO (STDOUT_HANDLE);
            std::cout << std::setw(7) << " ";
            if (normal) {
                std::cout << "min = " << min << ", max = " << max;
            }
            else {
                std::cout <<   "min = " << std::scientific << std::setprecision(6) << min     << std::fixed << std::setprecision(OUTPUT_PRECISION) 
                          << ", max = " << max << ", 1.0-max = " 
                                        << std::scientific << std::setprecision(6) << 1.0-max << std::fixed << std::setprecision(OUTPUT_PRECISION);
            }
            std::cout << ", mean +/- stddev = " << mean << " +/- " << stddev << " (+/-" 
                      << std::scientific << std::setprecision(1) << std_err << std::fixed << std::setprecision(OUTPUT_PRECISION) << ")";
            if (error) std::cout << " (error!)" << std::endl;  
            else       std::cout << " (ok)" << std::endl;  
            std::cout << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);

            if (rep == 0) {
                // 2nd round: take histogram again but from mapped distribution from uniform to normal or vice versa
                normal = !normal;
                info += " mapped to ";
                info += ((normal) ? "normal distribution" : "uniform distribution");
            }
        }
        
        delete [] count;
    }
    //std::cout << std::endl;
    
    return error;
}

// test conversion of to_double and from_double num repetitions
int test_to_double(random_generator<RNG_INT_TYPE> *rng, uint64_t num) {
    int error = 0;  

#if REAL_PRECISION == REAL_DOUBLE
    std::cout << "test_to_double performing " << num << " int64 <-> double conversions ... " << std::endl;
#else
    std::cout << "test_to_double performing " << num << " int32 <-> float conversions ... " << std::endl;
#endif
    std::cout << std::setfill('0') << std::fixed << std::setprecision(8);
    for (uint64_t k = 0; k < num; ++k) {
        RNG_INT_TYPE value = (*rng)();
        REAL_TYPE d = to_double(value);
        RNG_INT_TYPE i   = from_double(d);
        if (i != (value & RNG_MASK) ) {
            std::cout << "0x" << std::hex << std::setw(16) << value << std::dec << " <-> " << d << " error!"<< std::endl;
            error = -1101;
            break;
        }
        else {
            std::cout << "0x" << std::hex << std::setw(16) << value << std::dec << " <-> " << d << " ok" << std::endl;
        }
    }

    std::cout << std::endl;

    return error;
}

// test seed
// this generates num seed_seq objects and looks how the output pattern changes for each thread
// notes: 
// - seed_seq constructor takes only lower 32bits and returnes uint32_t values only!
// - each call so generate repeats the same same numbers! so repeated calls cannot be used for seeding threads.
int test_seed(random_generator<RNG_INT_TYPE> *rng, uint64_t num_tests, uint8_t num_threads) {
    int error = 0;
    std::cout << "test_seed performing " << num_tests << " tests ... " << std::endl;
    for (uint64_t i = 0; i < num_tests; ++i) {
        RNG_INT_TYPE value = (*rng)();
#if (REAL_PRECISION == REAL_DOUBLE)        
        std::seed_seq seq{value & 0xffffffff, value>>32};
#else
        std::seed_seq seq{value};
#endif
        std::vector<uint32_t> seed1(4);
#ifdef _DEBUG        
        // test repeated calls to generate: creates identical sequence of numbers!
        for (uint8_t k = 0; k < 2; ++k)
#endif
        {
            seq.generate(seed1.begin(), seed1.end());
            std::cout << std::setw(4) << i << ": " << std::hex << std::setw(16) << value;
            for (uint8_t j = 0; j < num_threads; ++j) {
                std::cout << " " << std::hex << std::setw(8) << seed1[j];
            }
            std::cout << std::dec << std::endl;
        }
    }
    
    std::cout << std::endl;
    
    return error;
}

// test uint128. returns 0 if ok, otherwise error code.
// notes: 
// - this code uses only double/uint64_t as data type.
// - TODO we must ensure this function is not optimized, but I am not sure if this is done properly?
#if defined(_WIN32) || defined(_WIN64) // windows
#pragma optimize("", off)
int test_uint128(uint64_t num, random_generator<uint32_t> & seed) {
#else
int __attribute__((optimize("O0")))
//#pragma GCC push_options
//#pragma GCC optimize ("-O0")
//#pragma GCC pop_options
test_uint128(uint64_t num, random_generator<uint32_t> & seed) {
#endif
    int error = 0;
#ifdef CLASS_UINT128    
    std::string info("test_uint128 (sim uint128)");
#else
    std::string info("test_uint128");
#endif

    // 0. check Lehmer128 multiplier (see random.hpp)
    if ((( Lehmer128_mult        & 0xffffffffffffffff) != uint64_t(0x2e714eb2b37916a5)) ||
        (((Lehmer128_mult >> 64) & 0xffffffffffffffff) != uint64_t(0x12e15e35b500f16e))) {
        error = -1150;
    }

    // test (a+b)*(c+d) == a*c + b*c + a*d + b*d
    random_generator<uint64_t>* rng = nullptr;
    error = get_RNG<double, uint64_t>(RND_MERSENNE_TWISTER, seed, NULL, NULL, &rng);
    if (!error) {
        if (rng == NULL) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << info << " error: could not get MersenneTwister!" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            error = -1151;
        }
        else {
            std::cout << info << ": " << num << " sum/mult test ... ";
            uint64_t value[8];
            uint64_t bits; // bits keeps the highest bits of each step and should ensure that compiler does not optimize steps away.  
            for (uint64_t k = 0; k < num; ++k) {
                (*rng)(value, 8);
#ifdef CLASS_UINT128
                uint128 a = uint128(value[0], value[1]);
                uint128 b = uint128(value[2], value[3]);
                uint128 c = uint128(value[4], value[5]);
                uint128 d = uint128(value[6], value[7]);
#else
                uint128 a = uint128(value[0]) | (uint128(value[1]) << 64);
                uint128 b = uint128(value[2]) | (uint128(value[3]) << 64);
                uint128 c = uint128(value[4]) | (uint128(value[5]) << 64);
                uint128 d = uint128(value[6]) | (uint128(value[7]) << 64);
#endif
                bits = ((a >> 127) & 1) | ((b >> 127) & 1) | ((c >> 127) & 1) | ((d >> 127) & 1);
                uint128 m = a + b;
                uint128 n = c + d;
                bits |= ((m >> 126) & 2) | ((n >> 126) & 2);
                uint128 x = m * n;
                bits |= ((x >> 125) & 4);
                uint128 o = a * c;
                uint128 p = b * c;
                uint128 q = a * d;
                uint128 r = b * d;
                bits |= ((o >> 124) & 8) | ((p >> 124) & 8) | ((q >> 124) & 8) | ((r >> 124) & 8);
                uint128 y = o + p + q + r;
                bits |= ((y >> 123) & 16);
                uint128 z = x - y;
                if (z != uint128(0)) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << "error:" << std::endl;
#ifdef CLASS_UINT128
                    std::cout << std::setw(6) << std::setfill('0') << k << " : ( " << std::hex << std::setw(16) <<
                                         a.hi << "_" << a.lo << "  +  " << b.hi << "_" << b.lo << " ) *" << std::endl <<
                        "         ( " << c.hi << "_" << c.lo << "  +  " << d.hi << "_" << d.lo << " ) =" << std::endl << 
                        "           " << y.hi << "_" << y.lo << " !=  " << x.hi << "_" << x.lo << " !"   << std::endl <<
                        "difference " << z.hi << "_" << z.lo << " !"                         << std::dec << std::endl;
#else
                    std::cout << std::setw(6) << std::setfill('0') << k << " : ( " << std::hex << std::setw(16) <<
                                         ((uint64_t)(a >> 64)) << "_" << ((uint64_t)(a & 0xffffffffffffffff)) << "  +  " << 
                                         ((uint64_t)(b >> 64)) << "_" << ((uint64_t)(b & 0xffffffffffffffff)) << " ) *"  << std::endl <<
                        "         ( " << ((uint64_t)(c >> 64)) << "_" << ((uint64_t)(c & 0xffffffffffffffff)) << "  +  " << 
                                         ((uint64_t)(d >> 64)) << "_" << ((uint64_t)(d & 0xffffffffffffffff)) << " ) ="  << std::endl <<
                        "           " << ((uint64_t)(y >> 64)) << "_" << ((uint64_t)(y & 0xffffffffffffffff)) << " !=  " << 
                                         ((uint64_t)(x >> 64)) << "_" << ((uint64_t)(x & 0xffffffffffffffff)) << " !"    << std::endl <<
                        "difference " << ((uint64_t)(z >> 64)) << "_" << ((uint64_t)(z & 0xffffffffffffffff)) << " !"    << std::dec  <<
                                                                                                                            std::endl;
#endif                        
                    std::cout << "bits: 0x" << std::hex << bits << std::dec << std::endl;
                    error = -1152;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    break;
                }
            }            
            if (!error) {
                SET_CONSOLE_COLOR_INFO(STDOUT_HANDLE);
                std::cout << "ok" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            }
            
            delete rng;
            rng = nullptr;
        }
    }
    
    std::cout << std::endl;

    return error;
}
#if defined(_WIN32) || defined(_WIN64) // windows
#pragma optimize("", on)
#endif

std::stringstream get_start_time(std::chrono::time_point<std::chrono::system_clock> &start) {
    start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    // works but output is not nice: std::ctime(&start_time)
    std::stringstream s_start;
    s_start << std::put_time(std::localtime(&start_time), FMT_DATE);
    return s_start;
}

std::stringstream get_end_time(std::chrono::time_point<std::chrono::system_clock> start, std::chrono::duration<REAL_TYPE> &duration) {
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    duration = end-start;
    std::stringstream s_end;
    s_end << std::put_time(std::localtime(&end_time), FMT_DATE);
    return s_end;
}

// command line argument options
static const std::string arg_gen     = "-g";
static const std::string arg_seed    = "-s";
static const std::string arg_param   = "-p";
static const std::string arg_test_tI = "-tI";
static const std::string arg_test_ti = "-ti";
static const std::string arg_test_tD = "-tD";
static const std::string arg_test_td = "-td";
static const std::string arg_test_tN = "-tN";
static const std::string arg_test_tn = "-tn";

int main(int argc, char* argv[])
{
    int error = 0;
    uint32_t t_total = GetTickCount();
    std::chrono::time_point<std::chrono::system_clock> tc_total, tc_start;
    std::stringstream ss_time = get_start_time(tc_total);
    std::chrono::duration<REAL_TYPE> tc_duration;
    struct thread_data *pthreads = NULL;
    class queue_entry *entry     = NULL;
    class queue *to_helper       = NULL;
    class queue *from_helper     = NULL;
    uint64_t count[3]            = {0,0,0};
    uint64_t num                 = 0;
    random_generator<uint32_t> *seed_generator = nullptr;   // RNG used for seed generator of threads
    std::vector     <uint32_t>  seed_values;                // seed values
    REAL_TYPE k_gamma            = 0.0;
    std::string filename; // parameter file name
    uint8_t rnd_index = RND_DEFAULT;
    bool seed_set  = false;
    bool gen_set   = false;
    bool param_set = false;
    uint8_t do_RNG_test = RNG_TEST_NONE;
    
    // mode_info printed each repetition and saved into info field of result file.
    std::string mode_info = VERSION_INFO;
#if REAL_PRECISION == REAL_SINGLE
    mode_info += ", single precision";
    SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
    std::cerr << std::endl << "* * * A T T E N T I O N ! * * *"                     << std::endl;
    std::cerr << argv[0] << " is compiled with 32bit single precision!"             << std::endl;
    std::cerr << "it is highly recommended to compile with 64bit double precision!" << std::endl;
    std::cerr << "use at your own risk and only for testing!"                       << std::endl;
    std::cerr << "to change this see MolecularConversion.h"                         << std::endl;
    std::cerr << "change definition of REAL_PRECISION to REAL_DOUBLE and recompile" << std::endl;
    std::cerr << std::endl;
    RESET_CONSOLE_COLOR(STDERR_HANDLE);
#endif
#ifdef USE_OLD_MOLSEARCH
    mode_info += ", v1.2 mol.search";
#endif
#ifdef CLASS_UINT128
    mode_info += ", sim uint128";
#else
#endif
#ifdef _DEBUG
    mode_info += ", debug";
#endif    
     
    // parse command line arguments, 0 = name of program
    for (uint8_t i = 1; (i < argc) && (!error); ++i) {
        //std::cerr << i << " " << argv[i] << std::endl; 
        if (arg_param.compare(argv[i]) == 0) {
            // read parameter file given in next argument
            if (param_set) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: parameter file given twice!?" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -110;
            } 
            else if (do_RNG_test) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: parameter file cannot be used with RNG test!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -111;
            }
            else if ((i+1) < argc ) {
                i += 1;
                filename = argv[i];
                param_set = true;
            }
            else {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: expected filename missing!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -112;
            }
        }
        else if (arg_gen.compare(argv[i]) == 0) {
            // select generator in next argument
            if (gen_set) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: generator set twice!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -120;
            }
            else if ((i+1) < argc ) {
                i += 1;
                uint8_t j = 0;
                for (; j < NUM_RND; ++j) {
                    if (argv[i] == rnd_all[j]) {
                        rnd_index = j;
                        rnd_use = rnd_all[rnd_index];
                        break; 
                    }
                }
                if ( j >= NUM_RND ) {
                    SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                    std::cerr << "error: generator name '" << argv[i] << "' invalid!" << std::endl;
                    RESET_CONSOLE_COLOR(STDERR_HANDLE);
                    std::cerr << "possible generators are:" << std::endl;
                    for (j = 0; j < NUM_RND; ++j) {
                        std::cerr << rnd_all[j] << std::endl;
                    }
#if (REAL_PRECISION == REAL_DOUBLE)
                    std::cerr << std::endl;
                    std::cerr << "note: " << argv[0] << " is compiled with 64bit double precision which disables 32bit generators!" << std::endl;
                    std::cerr << "      to enable 32bit generators compile with REAL_PRECISION = REAL_SINGLE but which is not recommended!"      << std::endl;
                    std::cerr << "      this is only for testing purpose! use at your own risk!"                                                 << std::endl;
#endif
                    error = -121;
                }
                else {
                    std::cerr << "using random number generator '" << rnd_use << "'" << std::endl;
                    gen_set = true;
                }
            }
            else {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: expected generator name missing!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -122;
            }
        }
        else if (arg_seed.compare(argv[i]) == 0) {
            // set new seed values
            if (seed_set) {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: seed set twice!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -130;
            }
            else if ((i+1) < argc ) {
                i += 1;
                seed_use = argv[i];
                seed_set = true;
            }
            else {
                SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
                std::cerr << "error: expected seed value missing!" << std::endl;
                RESET_CONSOLE_COLOR(STDERR_HANDLE);
                error = -131;
            }
        }
        else if (arg_test_tI.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_INT64;
        }
        else if (arg_test_ti.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_INT32;
        }
        else if (arg_test_tD.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_DOUBLE;
        }
        else if (arg_test_td.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_FLOAT;
        }
        else if (arg_test_tN.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_NORMAL_DOUBLE;
        }
        else if (arg_test_tn.compare(argv[i]) == 0) {
            do_RNG_test = RNG_TEST_NORMAL_FLOAT;
        }
        else {
           // unknown option
           error = -150;
           SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
           std::cerr << "unknown option " << argv[i] << std::endl << std::endl;
           RESET_CONSOLE_COLOR(STDERR_HANDLE);
           std::cerr << "usage:" << std::endl;
           std::cerr << "-g <RNG>  = select random number generator <RNG name>"                                          << std::endl;
           std::cerr << "-s <seed> = give seed value(s) as \"{1,2,..}\" (minimum 1 value)"                               << std::endl;
           std::cerr << "-p <file> = read options from parameter file <file name>"                                       << std::endl;
           std::cerr << "-tI       = stream uniform distributed 64bit integer to stdout (PractRand)"                     << std::endl;
           std::cerr << "-tD       = stream uniform distributed double converted to 64bit integer to stdout (PractRand)" << std::endl;
           std::cerr << "-tN       = stream normal  distributed double converted to 64bit integer to stdout (PractRand)" << std::endl;
           std::cerr << "-ti       = stream uniform distributed 32bit integer to stdout (PractRand)"                     << std::endl;
           std::cerr << "-td       = stream uniform distributed float converted to 32bit integer to stdout (PractRand)"  << std::endl;
           std::cerr << "-tn       = stream normal  distributed float converted to 32bit integer to stdout (PractRand)"  << std::endl;
           break;
        }
    }
    
    if ( (!error) && (do_RNG_test == RNG_TEST_NONE) && (filename.length() > 0) ) {
        
        // read parameters if parameter file is given as argument '-p [filename]'
        std::cout << "reading '" << filename << "' ..." << std::endl;
        error = ReadParams(filename.c_str(), argv);
        
        if ( (!error) && ((num_threads < 1) || (num_threads > THREADS_MAX_NUM)) ) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << "error: invalid number of threads " << +num_threads << " must be within [1," << THREADS_MAX_NUM << "]" << std::endl;
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            error = -160; 
        }
        else if ( num_bins < 1) {
            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
            std::cout << "error: number of bins must be >= 1!";
            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
            error = -161;
        }
        else if ( ( ::atom_num < 1) || ( ::atom_num > 2) ) {
            std::cout << "number of species must be 1 or 2!";
            error = -162;
        }
        else if ( ::atom_num == 1 ) {
            // check if any species 1 parameter is set
            for (uint8_t j = 0; j < NUM_VAR; ++j) {
                if ( (var_params[j].flag & FLAG_SP1) && (var_params[j].flag & FLAG_SET) ) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << "error: parameter '" << var_params[j].name << "' cannot be set with single species!" << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -163;
                    break;
                } 
            }
            if (atom_N[0] <= 0) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: N_0 must be > 0!" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -164;
            }
            else if (!error) {
                // reset unused species to be sure its not used somwhere by mistake
                ::atom_name[1] = atom_stat_use[1] = "";
                atom_mass[1] = 0.0; //atom_mass[0]; // needed for Ekin
                atom_N[1] = 0;
                atom_T[1] = 0.0;
                atom_wx[1] = 0.0; //atom_wx[0];
                atom_wy[1] = 0.0; //atom_wy[0];
                atom_wz[1] = 0.0; //atom_wz[0];
                CalculationSize[1] = 0;
            }
        }
        else {
            if (atom_N[0] <= 0) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: N 0 must be > 0!" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -170;
            }
            else if (atom_N[1] <= 0) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: N 1 must be > 0!" << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -171;
            }
        }
        
        if(!error) {
            for (uint8_t j = 0; j < MAX_NUM_SPECIES; ++j) {
                if ((::atom_num == 1) && (j == 1)) continue; // skip second species if not selected
                uint8_t i = 0;
                for (; i < NUM_STAT; ++i) {
                    if (atom_stat_use[j] == stat_all[i]) {
                        atom_stat_index[j] = i;
                        break; 
                    }
                }
                if (i >= NUM_STAT) {
                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                    std::cout << "error: statistics for " << ::atom_name[j] << " '" << atom_stat_use[j] << "' invalid! choose one of the following: " << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    error = -180;
                    for (i = 0; i < NUM_STAT; ++i) {
                        std::cout << stat_all[i] << std::endl;
                    }
                    std::cout << std::endl;
                    break;
                }
            }
        }

        if (!error) {
            uint8_t i = 0;
            for (; i < NUM_STAT; ++i) {
                if (mol_stat_use == stat_all[i]) {
                    mol_stat_index = i;
                    break;
                }
            }
            if (i >= NUM_STAT) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: statistics for molecule '" << mol_stat_use << "' invalid! choose one of the following: " << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -181;
                for (i = 0; i < NUM_STAT; ++i) {
                    std::cout << stat_all[i] << std::endl;
                }
                std::cout << std::endl;
            }
            else if (::atom_num == 1) { 
                // single species
                if ( mol_stat_index == STAT_FD ) {
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: '" << stat_all[mol_stat_index] << "' molecule (dimer) cannot be generated from sincle-species atom! '" << ::atom_name[0] << "' statistics is '" << atom_stat_use[0] << "'." << std::endl << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                }
                else if ( (mol_stat_index == STAT_BE) && (atom_stat_index[0] == STAT_MB) ) {
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: '" << stat_all[mol_stat_index] << "' molecule (dimer) expects both atoms '" << stat_all[STAT_BE] << "' or '" << stat_all[STAT_FD] << "', but selection for '" << ::atom_name[0] << "' is '" << atom_stat_use[0] << "' (single species)!" << std::endl << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                }
            }
            else {
                // double species
                if ( (mol_stat_index == STAT_FD) && ((atom_stat_index[0] == atom_stat_index[1]) || (atom_stat_index[0] == STAT_MB)) ) {
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: '" << stat_all[mol_stat_index] << "' molecule (dimer) expects one atom '" << stat_all[STAT_BE] << "' and other '" << stat_all[STAT_FD] << "', but selection for '" << ::atom_name[0] << "' is '" << atom_stat_use[0] << "' and for '" << ::atom_name[1] << "' is '" << atom_stat_use[1] << "'!" << std::endl << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                }
                else if ( (mol_stat_index == STAT_BE) && ((atom_stat_index[0] != atom_stat_index[1]) || (atom_stat_index[0] == STAT_MB) || (atom_stat_index[1] == STAT_MB)) ) {
                    SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                    std::cout << "warning: '" << stat_all[mol_stat_index] << "' molecule (dimer) expects both atoms '" << stat_all[STAT_BE] << "' or '" << stat_all[STAT_FD] << "', but selection for '" << ::atom_name[0] << "' is '" << atom_stat_use[0] << "' and for '" << ::atom_name[1] << "' is '" << atom_stat_use[1] << "'!" << std::endl << std::endl;
                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                }
            }
        }
        
        if (!error) {
            uint8_t i = 0;
            for (; i < NUM_RND; ++i) {
                if (rnd_use == rnd_all[i]) {
                    rnd_index = i;
                    break; 
                }
            }
            if (i >= NUM_RND) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: random number generator '" << rnd_use << "' invalid! choose one of the following: " << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -182;
                for (i = 0; i < NUM_RND; ++i) {
                    std::cout << rnd_all[i] << std::endl;
                }
                std::cout << std::endl;
            }
        }
        
        if (!error) {
            uint8_t i = 0;
            for (; i < NUM_DSTRB; ++i) {
                if (dstrb_use == dstrb_all[i]) {
                    dstrb_index = i;
                    break; 
                }
            }
            if (i >= NUM_DSTRB) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: distribution generator '" << dstrb_use << "' invalid! choose one of the following: " << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -183;
                for (i = 0; i < NUM_DSTRB; ++i) {
                    std::cout << dstrb_all[i] << std::endl;
                }
                std::cout << std::endl;
            }
        }

        if (!error) {
            uint8_t i = 0;
            for (; i < NUM_DISTANCE; ++i) {
                if (mol_distance_use == distance_all[i]) {
                    mol_distance_index = i;
                    break; 
                }
            }
            if (i >= NUM_DISTANCE) {
                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                std::cout << "error: distance '" << mol_distance_use << "' invalid! choose one of the following: " << std::endl;
                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                error = -184;
                for (i = 0; i < NUM_DISTANCE; ++i) {
                    std::cout << distance_all[i] << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }
    
    if (!error) {
        // get seed generator from ::seed_use or random device.
        error = get_seed(seed_generator, seed_values);
        if (!error) {
            if (do_RNG_test != RNG_TEST_NONE) {
                // stream random numbers on stdout for PractRand or other testing suites.
                // notes: 
                // - RNG_stream returns only on error!
                // - ensure that nothing is written to stdout before this test starts! std::cerr is ok.
                // - on windows console must be set into binary mode!

                error = RNG_stream(rnd_index, *seed_generator, do_RNG_test);
            }
            else {
                {
                    // perform random number self-tests
                    uint64_t num = 8192;
                    random_generator<RNG_INT_TYPE> *rng_int     = nullptr;
                    random_generator<REAL_TYPE   > *rng_uniform = nullptr;
                    random_generator<REAL_TYPE   > *rng_normal  = nullptr;
                    error = get_RNG<REAL_TYPE, RNG_INT_TYPE>(rnd_index, *seed_generator, &rng_uniform, &rng_normal, &rng_int);
                    if (error) {
                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                        std::cout << "error: could not get random number generator " << rnd_use << " (" << +rnd_index << ")!" << std::endl;
                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                    }
                    else {
                        error = test_RNG<RNG_INT_TYPE>(rng_int, num, 5, std::string("test ") + rnd_use);
                        if (!error) error = test_distribution(rng_uniform, num, false, 0.1, 4.1, 80, 10); 
                        if (!error) error = test_distribution(rng_normal , num, true , 0.1, 4.1, 80, 10); 
                        if (!error) error = test_to_double(rng_int, 5); 
                        if (!error) error = test_seed(rng_int, 5, 4);
                        if (!error) error = test_uint128(num, *seed_generator);
                    }
                    if (rng_int    ) delete rng_int;
                    if (rng_uniform) delete rng_uniform;
                    if (rng_normal ) delete rng_normal;
                }
    
                if(!error) {
                    // reset seed generator to original seed_values.
                    // this ensures that Monte Carlo is done with reproduceable seed values regardless of previous tests.
                    error = get_seed(seed_generator, seed_values);
                }
                
                if (!error) {

                    // set default floating point output style
                    std::cout << std::fixed << std::setprecision(OUTPUT_PRECISION); 
                    
                    if (!error) {
                        // init global mutex and barrier
                        if (!MUTEX_INIT_AND_OK(::global_mutex)) {
                            error = -200;
                        }
                        if (!BARRIER_INIT_AND_OK(::global_barrier, num_threads)) {
                            error = -201;
                        }
                    }
                    
                    if (!error) {
                        // save molecule name
                        mol_name = (::atom_num == 1) ? ::atom_name[0]+::atom_name[0] : ::atom_name[0]+::atom_name[1];

                        // convert parameters to appropriate units
                        //Mass[0] *= mn;
                        //Mass[1] *= mn;
                        atom_wx[0] *= 2*Pi;
                        atom_wx[1] *= 2*Pi;
                        atom_wy[0] *= 2*Pi;
                        atom_wy[1] *= 2*Pi;
                        atom_wz[0] *= 2*Pi;
                        atom_wz[1] *= 2*Pi;

                        if (mol_distance_index == DISTANCE_X) {
                            // squared real-space distance in um
                            k_gamma = mol_gamma*mol_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_V) {
                            // squared velocity distance in mm/s
                            k_gamma = mol_gamma*mol_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_PV) {
                            // squared velocity space distance in units of Planck constant divided by reduced mass
                            // this is the same as center-of-mass momentum distance
                            // notes: 
                            // - I get without factor 2.0 in calculation of relative momentum in center of mass frame.
                            //   but this is per particle, and maybe factor 2 is for both particles together?
                            //   comparison with Hodby figure 3 and Yamakoshi figure 2 shows that factor 2.0 should be there.
                            // - k_gamma = (mol_gamma*h/(ReducedMass*mn*1e-9))^2
                            REAL_TYPE ReducedMass = (::atom_num == 1) ? atom_mass[0] : 2.0*atom_mass[0]*atom_mass[1]/(atom_mass[0]+atom_mass[1]); // amu
                            k_gamma = mol_gamma*k_h/ReducedMass; // k_h = h/(mn*1e-9)
                            k_gamma = k_gamma*k_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_PP) {
                            // squared momentum-space distance in units of Planck constant
                            k_gamma = mol_gamma*k_h; // k_h = h/(mn*1e-9)
                            k_gamma = k_gamma*k_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_CROSS) {
                            // squared cross momentum-space distance |delta x  x  delta p| in units of reduced Planck constant
                            k_gamma = mol_gamma*k_h/(2*Pi); // k_h = h/(mn*1e-9)
                            k_gamma = k_gamma*k_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_MEAN) {
                            // cubic mean momentum-space distance: <|delta x| * |delta p|> in units of Planck constant
                            k_gamma = mol_gamma*k_h; // k_h = h/(mn*1e-9)
                            k_gamma = k_gamma*k_gamma*k_gamma;
                        }
                        else if (mol_distance_index == DISTANCE_MAX) {
                            // squared maximum momentum-space distance: max|delta xi * delta pi| in units of Planck constant
                            k_gamma = mol_gamma*k_h; // k_h = h/(mn*1e-9)
                            k_gamma = k_gamma*k_gamma;
                        }
                        else {
                            std::cout << "error distance measure ' " << distance_all[mol_distance_index] << " ' (" << mol_distance_index << ") not implemented!" << std::endl;
                            return -1;
                        }
                                
                        // diplay conversion constants to check order of magnitude
                        std::cout << "conversion constants:" << std::endl;
                        std::cout << std::scientific << std::setprecision(6); 
                        std::cout << "k_Epot         = " << k_Epot << std::endl;
                        std::cout << "k_Ekin         = " << k_Ekin << std::endl;
                        std::cout << "k_nK           = " << k_nK   << std::endl;
                        std::cout << "k_h            = " << k_h    << std::endl << std::endl;   

                        std::cout << "maximum " << mol_distance_use << " distance for molecule formation: (index " << +mol_distance_index << ")" << std::endl;     
                        std::cout << "gamma          = " << mol_gamma << std::endl;
                        std::cout << "gamma scaled^" << ((mol_distance_index==DISTANCE_MEAN)?"3":"2") <<" = " << k_gamma << std::endl << std::endl;
                        std::cout << std::fixed << std::setprecision(OUTPUT_PRECISION);
                        
                        // thread data and queues
                        pthreads    = new struct thread_data[num_threads];
                        to_helper   = new class queue;
                        from_helper = new class queue;

#ifdef _DEBUG        
                        if (to_helper->get_status() != QUEUE_STATUS_OK) {
                            std::cout << "queue 'to helper' init error!" << std::endl << std::endl;
                            //std::cout << to_helper->get_status() << std::endl;
                            if (!error) error = -210;
                        }
                        if (from_helper->get_status() != QUEUE_STATUS_OK) {
                            std::cout << "queue 'from helper' init error!" << std::endl << std::endl;
                            //std::cout << from_helper->get_status() << std::endl;
                            if (!error) error = -211;
                        }
#endif

                        // molecule constants
                        REAL_TYPE mol_mass = (::atom_num == 1) ? atom_mass[0]*2.0 : atom_mass[0] + atom_mass[1];
                        // molecule f_trap scaled from omega = std::sqrt(U/m)
                        REAL_TYPE mol_wx   = (::atom_num == 1) ? atom_wx[0] : std::sqrt((atom_wx[0]*atom_wx[0]*atom_mass[0]+atom_wx[1]*atom_wx[1]*atom_mass[1])/mol_mass);
                        REAL_TYPE mol_wy   = (::atom_num == 1) ? atom_wy[0] : std::sqrt((atom_wy[0]*atom_wy[0]*atom_mass[0]+atom_wy[1]*atom_wy[1]*atom_mass[1])/mol_mass);
                        // TODO: axial direction is often given by B-field curvature and not by optical confinenemt which we assume here.
                        //       this entire calculation can be done later (outside of this code) and is not really needed here.
                        REAL_TYPE mol_wz   = (::atom_num == 1) ? atom_wz[0] : std::sqrt((atom_wz[0]*atom_wz[0]*atom_mass[0]+atom_wz[1]*atom_wz[1]*atom_mass[1])/mol_mass);
                        REAL_TYPE mol_who  = std::pow(mol_wx*mol_wy*mol_wz,1.0/3.0);

                        if (!error) {                   
                            // init thread data
                            for (uint8_t i = 0; i < num_threads; ++i) { 
                                if (!MUTEX_INIT_AND_OK(pthreads[i].mutex)) {
                                    error = -220;
                                    break;
                                }
                                pthreads[i].handle      = INVALID_THREAD;
                                pthreads[i].id          = i;
                                pthreads[i].num_threads = num_threads;
                                pthreads[i].to_helper   = to_helper;
                                pthreads[i].from_helper = from_helper;
                                pthreads[i].next        = (i < (num_threads-1)) ? &pthreads[i+1] : &pthreads[0];
#if defined(_WIN32) || defined (_WIN64)
                                pthreads[i].win_id      = 0;
#endif
                                for (uint8_t j = 0; j < NUM_CELLS_LISTS; ++j) {
                                    pthreads[i].num   [j] = 0;
                                    pthreads[i].mass  [j] = 0.0;
                                    pthreads[i].x     [j] = nullptr;
                                    pthreads[i].E     [j] = nullptr;
                                    pthreads[i].status[j] = nullptr;
                                } 
                                
                                // get random number generators seed individually for each thread
                                // note: on error do not break loop but wait until all threads are started.
                                if (get_RNG<REAL_TYPE, RNG_INT_TYPE>(
                                        rnd_index, 
                                        *seed_generator, 
                                        &pthreads[i].gen_uniform, 
                                        &pthreads[i].gen_normal, 
                                        NULL)) {
                                    if (!error) error = -221;
                                }
                                // insert first message into queue. 
                                // note: we insert one message for each thread but any running thread can pick up message and respond.
                                //       but since each thread waits for barrier after responding, so they should respond 1x each of them.
#ifdef _DEBUG
                                if (to_helper->put(new class queue_entry(THREAD_START, NULL)) != QUEUE_STATUS_OK) {
                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                    std::cout << std::endl << "error: queue 'to_helper' error 'put START'!" << std::endl << std::endl;
                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                    if (!error) error = -222;
                                }
#else
                                to_helper->put(new class queue_entry(THREAD_START, NULL));
#endif
                            }
                        }

                        if (!error) {
                            // create all cells on main thread and alternate owner thread over pdata list.
                            // TODO: cells could be created in parallel by all threads.
                            // 1. first and second species
                            uint8_t cells_count[3] = {0,0,0}; // number of cells per species0,species1,molecule
                            uint8_t cell_index = 0; // cell index
                            size_t max_name = 0;
                            REAL_TYPE *max = NULL, *R = NULL;
                            for (uint8_t species = 0; species < ::atom_num; ++species) {
                                if (atom_stat_index[species] == STAT_BE) {
                                    // Bosons: create separate species cells for non-condensed and condensed parts
                                    // get critical temperature and chemical potential
                                    atom_mu[species] = get_mu(species);
                                    // get calculation volume
                                    get_size(species, STAT_BE_NC, &max, &R);
                                    // name
                                    std::string name = ::atom_name[species] + "(thermal)";
                                    if (name.length() > max_name) max_name = name.length();
                                    // create cells
                                    ::cells[cell_index] = new cells_list(pthreads,
                                                CELLS_TYPE_ATOMS,
#if NUM_CELLS > 1
                                                NUM_CELLS, 
#else
                                                num_threads,
#endif
                                                atom_N[species]-atom_Nc[species],
                                                species, SPECIES_NONE,
                                                name, 
                                                STAT_BE_NC, 
                                                //distance_default[STAT_BE_NC], 
                                                mol_distance_index,
                                                //atom_gamma[species], 
                                                k_gamma,
                                                //false,
                                                max,
                                                R);
                                    ++cells_count[species];
                                    ++cell_index;
                                    if (atom_Nc[species] > 0) {
                                        // for BEC with rejection sampling need a much narrower range otherwise too slow.
                                        // therefore, we recalculate size for condensed part.
                                        // TODO: this calculation might be useful also for other statistics & Metropolis.
                                        get_size(species, STAT_BE_C, &max, &R);
                                        // name
                                        std::string name = ::atom_name[species] + "(condensed)";
                                        if (name.length() > max_name) max_name = name.length();
                                        ::cells[cell_index] = new cells_list(pthreads,
                                                CELLS_TYPE_ATOMS, 
#if NUM_CELLS > 1
                                                NUM_CELLS,
#else
                                                num_threads,
#endif
                                                atom_Nc[species],
                                                species, SPECIES_NONE,
                                                name, 
                                                STAT_BE_C, 
                                                //distance_default[STAT_BE_C], 
                                                mol_distance_index,
                                                //atom_gamma[species], 
                                                k_gamma,
                                                //false,
                                                max,
                                                R);
                                        ++cells_count[species];
                                        ++cell_index;
                                    }
                                }
                                else {
                                    // Fermions or thermal gas: single list of cells 
                                    // get critical temperature and chemical potential
                                    atom_mu[species] = get_mu(species);
                                    // get calculation volume
                                    get_size(species, atom_stat_index[species], &max, &R);
                                    // name
                                    std::string name = ::atom_name[species];
                                    if (name.length() > max_name) max_name = name.length();
                                    // create cells
                                    ::cells[cell_index] = new cells_list(pthreads,
                                                CELLS_TYPE_ATOMS, 
#if NUM_CELLS > 1
                                                NUM_CELLS,
#else
                                                num_threads,
#endif
                                                atom_N[species],
                                                species, SPECIES_NONE,
                                                name, 
                                                atom_stat_index[species], 
                                                //distance_default[atom_stat_index[species]], 
                                                mol_distance_index,
                                                //atom_gamma[species], 
                                                k_gamma,
                                                //false,
                                                max,
                                                R);
                                    ++cells_count[species];
                                    ++cell_index;
                                }
                                    
                                // output species specific information
                                std::cout << "N        " << ::atom_name[species] << " = " << atom_N[species]                     << std::endl;
                                std::cout << "T        " << ::atom_name[species] << " = " << atom_T[species]            << " nK" << std::endl;
                                if (atom_stat_index[species] == STAT_FD) {
                                    std::cout << "TF       " << ::atom_name[species] << " = " << atom_Tcrit[species]   << " nK" << std::endl;
                                    std::cout << "T/TF     " << ::atom_name[species] << " = " << atom_T[species]/atom_Tcrit[species] << std::endl;
                                    std::cout << "Rx       " << ::atom_name[species] << " = " << R[0] << " um" << std::endl;
                                    std::cout << "Ry       " << ::atom_name[species] << " = " << R[1] << " um" << std::endl;
                                    std::cout << "Rz       " << ::atom_name[species] << " = " << R[2] << " um" << std::endl;
                                }
                                else if ((atom_stat_index[species] == STAT_BE   ) || 
                                            (atom_stat_index[species] == STAT_BE_NC)) {
                                    std::cout << "Tc       " << ::atom_name[species] << " = " << atom_Tcrit[species]   << " nK" << std::endl;
                                    std::cout << "T/Tc     " << ::atom_name[species] << " = " << atom_T[species]/atom_Tcrit[species] << std::endl;
                                    std::cout << "Nc/N     " << ::atom_name[species] << " = " << atom_Nc[species] << " / " << atom_N[species] << " = " << (((REAL_TYPE)atom_Nc[species])/atom_N[species]) << std::endl;
                                    std::cout << "Rx       " << ::atom_name[species] << " = " << R[0] << " um" << std::endl;
                                    std::cout << "Ry       " << ::atom_name[species] << " = " << R[1] << " um" << std::endl;
                                    std::cout << "Rz       " << ::atom_name[species] << " = " << R[2] << " um" << std::endl;
                                }
                                else if (atom_stat_index[species] == STAT_BE_C) {
                                    std::cout << "Tc       " << ::atom_name[species] << " = " << atom_Tcrit[species]   << " nK" << std::endl;
                                    std::cout << "T/Tc     " << ::atom_name[species] << " = " << atom_T[species]/atom_Tcrit[species] << std::endl;
                                    std::cout << "Nc/N     " << ::atom_name[species] << " = " << atom_Nc[species] << " / " << atom_N[species] << " = " << (((REAL_TYPE)atom_Nc[species])/atom_N[species]) << std::endl;
                                    std::cout << "sigma_x  " << ::atom_name[species] << " = " << R[0] << " um" << std::endl;
                                    std::cout << "sigma_y  " << ::atom_name[species] << " = " << R[1] << " um" << std::endl;
                                    std::cout << "sigma_z  " << ::atom_name[species] << " = " << R[2] << " um" << std::endl;
                                }
                                else if (atom_stat_index[species] == STAT_MB ) {
                                    std::cout << "sigma_x  " << ::atom_name[species] << " = " << R[0] << " um" << std::endl;
                                    std::cout << "sigma_y  " << ::atom_name[species] << " = " << R[1] << " um" << std::endl;
                                    std::cout << "sigma_z  " << ::atom_name[species] << " = " << R[2] << " um" << std::endl;
                                }
                                std::cout << "Emax     " << ::atom_name[species] << " = " << max_E[species] << " nK" << std::endl<< std::endl;                    
                            }
                            // 2. molecule: create 1d cells per thread. 
                            // additional cells needed for all combinations of normal & condensed parts.
                            uint8_t sp0 = 0;  
                            uint8_t sp1 = (::atom_num == 1) ? 0 : 1;
                            for (uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                                if (::cells[i]) {
                                    if ((::cells[i]->get_type() == CELLS_TYPE_ATOMS) && (::cells[i]->get_index(0) == sp0)) { 
                                        for (uint8_t j = i; j < NUM_CELLS_LISTS; ++j) {
                                            if (::cells[j]) {
                                                if ((::cells[j]->get_type() == CELLS_TYPE_ATOMS) && (::cells[j]->get_index(0) == sp1)) { 
                                                    // name
                                                    std::string name = ::cells[i]->get_species_name()+::cells[j]->get_species_name();
                                                    if (name.length() > max_name) max_name = name.length();
                                                    ::cells[cell_index] = new cells_list(pthreads, 
                                                        CELLS_TYPE_MOLECULES,
                                                        (uint64_t)num_threads,
                                                        0, // initial num_atoms = 0 
                                                        i, j, // molecule between ::cells indices i and j 
                                                        name,
                                                        mol_stat_index, 
                                                        mol_distance_index, 
                                                        k_gamma, 
                                                        //false,
                                                        NULL,
                                                        NULL);
                                                    ++cells_count[2];
                                                    ++cell_index;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (cell_index > NUM_CELLS_LISTS) {
                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                std::cout << "error: number of cells list " << +cell_index << " exceeded maximum " << NUM_CELLS_LISTS << "! expect memory overflow and SEGFAULT!" << std::endl;
                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                error = -223;
                            }

                            if (!error) {
                                std::cout << "using " << +cell_index << " lists of cells:" << std::endl; 
                                for (uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                                    if (::cells[i]) {
                                        std::cout << "#" << std::setw(2) << std::setfill(' ') << std::right << +i << " " <<
                                        std::setw(max_name) << std::setfill(' ') << std::left << cells[i]->get_species_name() <<
                                        ": created " << std::setw(1) << +cells[i]->get_cells_per_dim() << "^" << +cells[i]->get_dim() << " = " << std::setw(7) << cells[i]->get_num() << std::setw(0) << +cells[i]->get_dim() << "d cells using " <<
                                        std::setw(distance_all[DISTANCE_MAX].length()+2) << std::setfill(' ') << std::left <<
                                        ("'"+distance_all[cells[i]->get_distance_index()]+"'") <<
                                        //(cells[i]->get_duplicate() ? " (duplicate)" : "") << 
                                        std::endl;
                                    }
                                }
                            }
                            /*
                            if (distance_default[atom_stat_index[0]] != mol_distance_index)
                                ::cells[CELLS_MOL_SP0] = new cells_list(pthreads, NUM_CELLS, 0, ::atom_name[0], mol_distance_index, k_gamma, true);
                            if (distance_default[atom_stat_index[1]] != mol_distance_index)
                                ::cells[CELLS_MOL_SP1] = new cells_list(pthreads, NUM_CELLS, 1, ::atom_name[1], mol_distance_index, k_gamma, true);
                            */
                            
                            if (!error) {
                                // start threads
                                // we wait until all of them are started and reply to THREAD_START message
                                std::cout << std::endl <<"calculation with " << +num_threads << " threads" << std::endl;
                                for (uint8_t i = 0; i < num_threads; ++i) { 
                                    int tmp = thread_start(helper_thread_func, &pthreads[i].handle, (void*)&pthreads[i]);
                                    if (tmp) {
                                        std::cout << "thread " << +entry->id << "/" << +num_threads << " creation error " << tmp << std::endl;
                                        if (!error) error = tmp;
                                    }
                                    else {
                                        entry = from_helper->get(1, 2*THREAD_TIMEOUT);
                                        if (entry) {
                                            if (entry->result) {
                                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                std::cout << "error: thread " << +entry->id << "/" << +num_threads << " created but responds is " << entry->result << " !?" << std::endl;
                                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                if (!error) error = -230;
                                            }
                                            else {
                                                std::cout << "thread " << +entry->id << "/" << +num_threads << " started ok" << std::endl;
                                            }
                                            delete entry;
                                        }
#ifdef _DEBUG
                                        else if (from_helper->get_status() != QUEUE_STATUS_OK) { // error
                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                            std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get START'!" << std::endl << std::endl;
                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                            if (!error) error = -231;
                                        }
#endif
                                        else { // timeout
                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                            std::cout << "error: thread " << i << "/" << num_threads << " created but no responds!?" << std::endl;
                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                            if (!error) error = -233;
                                        }
                                    }
                                }
                                //std::cout << +num_threads << " threads started" << std::endl; 
                            }
                             
                            if (!error) {
                    
                                // open results file                   
                                std::ofstream out; 
                                out.open(file_result, std::ios::app);
                                if (out.fail()) {
                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                    std::cout << std::endl << "error: could not open result file '" << file_result << "'!" << std::endl;
                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                    show_cwd(argv);
                                    error = -151;
                                }
                                else {
                                    REAL_TYPE sum = 0.0, sum_sqr = 0.0; // molecule conversion efficiency sum for mean and standard deviation
                                    REAL_TYPE t_sum = 0.0, t_sum_sqr = 0;
                                    uint32_t n_rep;
                                    for (n_rep = 1; (!error) && (n_rep <= repetitions); n_rep++) {
                                        uint32_t t_loop = GetTickCount(), t_start;
                                        ss_time = get_start_time(tc_start);
                                        
                                        std::cout << std::endl << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
                                        std::cout << ss_time.str() << std::endl;
                                        std::cout << "repetition " << n_rep << " / " << repetitions << " on " << +num_threads << " threads ... (" << mode_info << ")" << std::endl << std::endl;

                                        out << std::fixed << std::setprecision(FILE_PRECISION);
                                        out     << "////////////////////////////////////////////////////////////////////////////////"  << std::endl;
                                        out     << ss_time.str()                                                                       << std::endl;
                                        out     << "repetition                 = " << n_rep << " / " << repetitions                    << std::endl;
                                        out     << "number species             = " << +::atom_num                                      << std::endl;
                                        out     << "number threads             = " << +num_threads                                     << std::endl;
                                        out     << "number bins                = " << num_bins                                         << std::endl;
                                        out     << "disribution generator      = " << dstrb_use                                        << std::endl;
                                        out     << "random number generator    = " << rnd_use                                          << std::endl;
                                        if (seed_use.size() == 0) {
                                            out << "seed (hardware)            = {";
                                            for (uint32_t i = 0; i < seed_values.size(); ++i) {
                                                if (i == 0) out << "0x"  << std::hex << seed_values[i] << std::dec;
                                                else        out << ",0x" << std::hex << seed_values[i] << std::dec;    
                                            }
                                            out << "}" << std::endl;
                                        }
                                        else {
                                        out     << "seed (manual)              = " << seed_use                                         << std::endl;
                                        }
                                        out     << "size 0                     = " << CalculationSize[0] << " * " << EMAX_SCALING      << std::endl;
                                        if (::atom_num > 1) {
                                            out << "size 1                     = " << CalculationSize[1] << " * " << EMAX_SCALING      << std::endl;
                                        }
                                        out     << "step size                  = " << InitialStepSize                                  << std::endl;
                                        out     << "error allowed              = " << std::scientific << ErrorAllowed << std::fixed    << std::endl;
                                        out     << "find nearest               = " << +find_nearest                                    << std::endl;
                                        out     << "distance measure           = " << mol_distance_use                                 << std::endl;
                                        out     << "gamma                      = " << std::scientific << mol_gamma << std::fixed       << std::endl;
                                        out     << "gamma scaled               = " << std::scientific << k_gamma   << std::fixed       << std::endl;
                                        out     << "info                       = " << mode_info                                        << std::endl;
                                        out     << std::endl;

                                        // create atom distribution on num_threads threads
                                        std::cout << "random number generator = '" << rnd_use  << "'" << std::endl;
                                        std::cout << "distribution generator  = '" << dstrb_use << "'" << std::endl; 
                                        
                                        REAL_TYPE atom_t[MAX_NUM_SPECIES] = { 0.0, 0.0 };
                                        for (uint8_t species = 0; (!error) && (species < ::atom_num); ++species) {
                                            for (cell_index = 0; cell_index < NUM_CELLS_LISTS; ++cell_index) {
                                                if (!::cells[cell_index]) continue; 
                                                if ((::cells[cell_index]->get_type() != CELLS_TYPE_ATOMS) || (::cells[cell_index]->get_index(0) != species)) continue;
                                                t_start = GetTickCount();
                                                uint64_t num_atoms = ::cells[cell_index]->get_num_atoms(); // total number of atoms to create.
                                                uint64_t atom_index = 0, atoms_delta = (uint64_t)(num_atoms/num_threads);
                                                std::cout << std::endl << "creating " << num_atoms << " " << ::cells[cell_index]->get_species_name() << " atoms on " << +num_threads << " threads ..." << std::endl;
                                                for (uint8_t i = 0; i < num_threads; ++i, atom_index += atoms_delta) { 
                                                    // number of atoms to create per species: int(N[i]/num_threads). for last we take remaining.
                                                    struct cmd_data *data  = new struct cmd_data;
                                                    data->in_index [0]     = species;
                                                    data->in_index [1]     = ::cells[cell_index]->get_statistics();
                                                    data->out_index[0]     = cell_index; // save atoms into ::cells[cell_index]
                                                    data->out_index[1]     = 0; 
                                                    data->count    [0]     = atom_index; // starting atom index to assign
                                                    data->count    [1]     = 0;
                                                    data->count    [2]     = 0;
                                                    data->total            = ( i < (num_threads-1) ) ? atoms_delta : num_atoms - atoms_delta*i;
                                                    for (uint8_t j = 0; j < DIM_ENERGY; ++j) data->E[j] = 0.0;                        
#ifdef _DEBUG                
                                                    if (to_helper->put(new class queue_entry(dstrb_index == DSTRB_BASIC ? THREAD_CREATE_BASIC : THREAD_CREATE_METROPOLIS, data)) != QUEUE_STATUS_OK) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: queue 'to_helper' error 'put create atoms'!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -240;
                                                    }
#else
                                                    to_helper->put(new class queue_entry(dstrb_index == DSTRB_BASIC ? THREAD_CREATE_BASIC : THREAD_CREATE_METROPOLIS, data));
#endif
                                                }
                                                // wait until atoms are created
                                                //std::atomic_thread_fence(std::memory_order_release)
                                                count[0] = count[1] = count[2] = num = 0;
                                                REAL_TYPE Emax_atoms = 0.0;
                                                for (uint8_t i = 0; i < num_threads; ++i) { 
                                                    entry = from_helper->get(1, WAIT_INFINITE);
                                                    if (entry) {
                                                        if (entry->result) {
                                                            if (!error) error = entry->result;
                                                            std::cout << "thread " << +entry->id << " creating atoms error " << entry->result << std::endl;
                                                        }
                                                        else if (entry->data) {
                                                            // ok: count number of atoms
                                                            struct cmd_data *data = reinterpret_cast<struct cmd_data*>(entry->data);
                                                            num           += data->total;     // number of created atoms
                                                            count[0]      += data->count[0];  // loop counter
                                                            count[1]      += data->count[1];  // E < Emin counter for stat_BE
                                                            count[2]      += data->count[2];  // E > Emax counter
                                                            if (data->E[E_MAX] > Emax_atoms) Emax_atoms = data->E[E_MAX]; // maximum energy of created atoms
                                                            delete data;
                                                        }
                                                        else {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " no atom data from thread?!" << std::endl << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -241;
                                                        }
                                                        delete entry;
                                                    }
#ifdef _DEBUG
                                                    else if (from_helper->get_status() != QUEUE_STATUS_OK) { // error
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get create FD'!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -242;
                                                    }
                                                    else { // timeout unexpected!
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get create FD' timeout!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -243;
                                                    }
#else
                                                    else { // there should be no timeout!
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " timeout get atom data?!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -244;
                                                    }
#endif                            
                                                }
                                                if (error) break;

                                                // check total number of created atoms 
                                                if ( num != num_atoms ) {
                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                    std::cout << "error: creating " << num_atoms << " " << ::cells[cell_index]->get_species_name() << " atoms, but " << num << " counted!" << std::endl;
                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    error = -245;
                                                    break;
                                                }

                                                // check atoms in cell
                                                num = 0;
                                                class cell_data *cell = cells[cell_index]->get_first();
                                                while(cell) {
#ifdef _DEBUG
                                                    // ensure all atoms are in proper cells, i.e. each cell index == atom hash key
                                                    // here we also check cell->get_num() vs. counted number of atoms.
                                                    uint64_t index = cell->get_hash_key();
                                                    uint64_t tmp = 0;
                                                    class species_entry *next = cell->species.get_first();
                                                    if (dstrb_index == DSTRB_BASIC) {
                                                        // for 'basic' distribution generator only first cell is populated, i.e. index != hash key
                                                        while(next) {
                                                            ++tmp;
                                                            next = next->get_next();
                                                        }
                                                    }
                                                    else {
                                                        while(next) {
                                                            if (next->get_hash_key() != index) {
                                                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                                std::cout << std::endl << "error: " << ::cells[cell_index]->get_species_name() << ": cell index does not match atom hash key: " << index << " != " << next->get_hash_key() << " !" << std::endl << std::endl;
                                                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                                error = -246;
                                                                break;
                                                            }
                                                            ++tmp;
                                                            next = next->get_next();
                                                        }
                                                        if (error) break;
                                                    }
                                                    if (tmp != cell->species.get_num()) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: " << ::cells[cell_index]->get_species_name() << ": cell with wrong atom number: " << tmp << " != " << cell->species.get_num() << " !" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -247;
                                                        break;
                                                    }
                                                    num += tmp;
#else
                                                    // count only atom number per species
                                                    num += cell->species.get_num();
#endif
                                                    cell = cell->get_next();
                                                }
                                                if (error) break;

                                                if (num != num_atoms) {
                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                    std::cout << std::endl << "error: " << ::cells[cell_index]->get_species_name() << ": unexpected atom number in cell " << num << " != " << num_atoms << " !" << std::endl << std::endl;
                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    error = -248;
                                                    break;
                                                }
                                                
                                                if (count[2] || (Emax_atoms >= max_E[species])) {
                                                    if ((::dstrb_index == DSTRB_METROPOLIS) && (count[2] < EMAX_ERROR_NUM_ATOMS)) {
                                                        // for Metropolis neigher E nor p need to be normalized. 
                                                        // we accept E > Emax for fewer than EMAX_ERROR_NUM_ATOMS.
                                                        SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                                                        std::cout << std::endl << "warning: ";
                                                    }
                                                    else {
                                                        // for rejection-sampling E > Emax is not acceptable.
                                                        // we keep this strict since rejection-sampling is intended for testing and not for production.
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: "; 
                                                        error = -249;
                                                    }
                                                    std::cout << "E = " << Emax_atoms << " nK > Emax = " << max_E[species] << " nK for " << count[2] << " atoms! increase calc_size[" << +species << "] = " << CalculationSize[species] << " for " << ::cells[cell_index]->get_species_name() << " (" << +species << ") !" << std::endl << std::endl;
                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    if (error) break;
                                                }
                                                else {
                                                    // suggest possible reduction of max_E when actual max. E < EMIN_SCALING*max_E
                                                    // we suggest scaling in powers of 2 which is a bit coarse.
                                                    // TODO: for BEC this might suggest different scaling for condensed and thermal parts.
                                                    if (Emax_atoms < (EMIN_SCALING*max_E[species])) {
                                                        REAL_TYPE sc = std::pow(2.0, std::ceil(std::log2(CalculationSize[species]*Emax_atoms/max_E[species])));
                                                        SET_CONSOLE_COLOR_NOTE(STDOUT_HANDLE);
                                                        std::cout << "note: " << ::cells[cell_index]->get_species_name() << " max. E = " << Emax_atoms << " nK < " << EMIN_SCALING << " * Emax = " << max_E[species] << " nK, ratio " << Emax_atoms/max_E[species] << ", calc_size[" << +species << "] = " << CalculationSize[species] << " suggested reduction to 1/" << 1.0/sc << " = " << sc << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    }
                                                    else {
                                                        std::cout << "info: " << ::cells[cell_index]->get_species_name() << " max. E = " << Emax_atoms << " nK < Emax = " << max_E[species] << " nK (ratio " << Emax_atoms/max_E[species] << "), calc_size[" << +species << "] = " << CalculationSize[species] << " (ok) " << std::endl;
                                                    }
                                                }
                                                
                                                if (count[1]) {
                                                    std::cout << "info: E < Emin counter = " << count[1] << " (ok)" << std::endl;
                                                } 

                                                if (file_histogram.length() > 0) { 
                                                    // save histograms on primary thread
                                                    std::string species_name = ::cells[cell_index]->get_species_name();
                                                    std::string file_name = expand_filename(file_histogram, species_name, n_rep);
                                                    Histogram(file_name, species_name, stat_all[::cells[cell_index]->get_statistics()], 
                                                        n_rep, //rep, 
                                                        num_bins, //bins, 
                                                        cell_index,
                                                        atom_T[species], atom_mu[species],
                                                        atom_wx[species], atom_wy[species], atom_wz[species],
                                                        atom_mass[species],
#if (HIST_GET_MAX)                      
                                                        false
#else
                                                        true //atom_max[species]
#endif
                                                        );
                                                }

                                                // atom creation time
                                                REAL_TYPE dt = ((REAL_TYPE)get_ticks_delta(t_start))/1e3;
                                                atom_t[species] += dt;
                                                std::cout << "created " << num_atoms << " " << ::cells[cell_index]->get_species_name() << " atoms with 1:" << (uint64_t)std::ceil(((REAL_TYPE)count[0])/num_atoms) << " efficiency (" << dt << " s)" << std::endl;

                                                out     << "atomic species             = " << ::cells[cell_index]->get_species_name()           << std::endl;
                                                out     << "statistics                 = " << stat_all[::cells[cell_index]->get_statistics()]   << std::endl;
                                                out     << "calculation time           = " << dt                           << " s"     << std::endl;
                                                out     << "N                          = " << ::cells[cell_index]->get_num_atoms()              << std::endl;
                                                out     << "T                          = " << atom_T[species]              << " nK"    << std::endl;
                                                if (atom_stat_index[species] == STAT_FD) {
                                                    out << "Rx                         = " << ::cells[cell_index]->get_R(0)        << " um"    << std::endl;
                                                    out << "Ry                         = " << ::cells[cell_index]->get_R(1)        << " um"    << std::endl;
                                                    out << "Rz                         = " << ::cells[cell_index]->get_R(2)        << " um"    << std::endl;
                                                }
                                                else if ((atom_stat_index[species] == STAT_BE   ) || 
                                                         (atom_stat_index[species] == STAT_BE_NC)) {
                                                    out << "Rx                         = " << ::cells[cell_index]->get_R(0)        << " um"    << std::endl;
                                                    out << "Ry                         = " << ::cells[cell_index]->get_R(1)        << " um"    << std::endl;
                                                    out << "Rz                         = " << ::cells[cell_index]->get_R(2)        << " um"    << std::endl;
                                                }
                                                else if (atom_stat_index[species] == STAT_BE_C) {
                                                    out << "sigma_x                    = " << ::cells[cell_index]->get_R(0)        << " um"    << std::endl;
                                                    out << "sigma_y                    = " << ::cells[cell_index]->get_R(1)        << " um"    << std::endl;
                                                    out << "sigma_z                    = " << ::cells[cell_index]->get_R(2)        << " um"    << std::endl;
                                                }

                                                else if (atom_stat_index[species] == STAT_MB) {
                                                    out << "sigma_x                    = " << ::cells[cell_index]->get_R(0)        << " um"    << std::endl;
                                                    out << "sigma_y                    = " << ::cells[cell_index]->get_R(1)        << " um"    << std::endl;
                                                    out << "sigma_z                    = " << ::cells[cell_index]->get_R(2)        << " um"    << std::endl;
                                                }
                                                for (uint8_t k = 0; k < DIM_MAXVAL; ++k) {
                                                    out << "max " << std::setw(22) << std::left << cname[k] << std::setw(0) << " = " << ::cells[cell_index]->get_max(k) << " " << cunit[k] << std::endl;
                                                }
                                                out     << "counters                   = [" << count[0] << "," << count[1] << "," << count[2] << "]" << std::endl;
                                                out     << std::endl;

                                            } // next cell
                                            if (error) break;

                                            // summary of species
                                            //if (cells_count[species] > 1) {
                                            std::cout << std::endl << "created " << atom_N[species] << " " << ::atom_name[species] << " atoms total (" << atom_t[species] << " s)" << std::endl;
                                            out     << "atomic species             = " << ::atom_name    [species] << " [summary]"   << std::endl;
                                            out     << "statistics                 = " << atom_stat_use[species]                   << std::endl;
                                            out     << "calculation time           = " << atom_t       [species]       << " s"     << std::endl;
                                            out     << "mass                       = " << atom_mass    [species]       << " amu"   << std::endl;
                                            out     << "N                          = " << atom_N       [species]                   << std::endl;
                                            out     << "T                          = " << atom_T       [species]       << " nK"    << std::endl;
                                            if (atom_stat_index[species] == STAT_FD) {
                                                out << "TF                         = " << atom_Tcrit[species]          << " nK"    << std::endl;
                                                out << "T/TF                       = " << atom_T[species]/atom_Tcrit[species]      << std::endl;
                                            }
                                            else if ((atom_stat_index[species] == STAT_BE   ) || 
                                                     (atom_stat_index[species] == STAT_BE_NC) || 
                                                     (atom_stat_index[species] == STAT_BE_C )) {
                                                out << "Tc                         = " << atom_Tcrit[species]          << " nK"    << std::endl;
                                                out << "T/Tc                       = " << atom_T[species]/atom_Tcrit[species]      << std::endl;
                                                out << "Nc/N                       = " << atom_Nc[species] << " / " << atom_N[species] << " = " << (((REAL_TYPE)atom_Nc[species])/atom_N[species]) << std::endl;
                                            }
                                            out     << "fx                         = " << atom_wx[species]/(2*Pi)      << " Hz"    << std::endl;
                                            out     << "fy                         = " << atom_wy[species]/(2*Pi)      << " Hz"    << std::endl;
                                            out     << "fz                         = " << atom_wz[species]/(2*Pi)      << " Hz"    << std::endl;
                                            if (species == 1) {
                                                out << "displacement x             = " << atom_cx                      << " um"    << std::endl;
                                                out << "displacement y             = " << atom_cy                      << " um"    << std::endl;
                                                out << "displacement z             = " << atom_cz                      << " um"    << std::endl;
                                            } 
                                            out     << "chemical potential         = " << atom_mu[species]             << " nK"    << std::endl;
                                            out     << std::endl;

                                        } // next species
                                        
                                        if (!error) {
                                                        
                                            // create molecules on num_threads threads
                                            // we loop through cells and create molecules for all combinations of species
                                            REAL_TYPE ReductionOfEnergy = 0, PotentialEnergyOfMolecules = 0, KineticEnergyOfMolecules = 0;
                                            uint64_t mol_N = 0;
                                            REAL_TYPE mol_t = 0.0, mol_E = 0.0;
                                            uint8_t checked_counter = 0;
                                            uint8_t atom_index[2], mol_index; // atom cell index for molecule mol_index
                                            for (mol_index = 0; mol_index < NUM_CELLS_LISTS; ++mol_index) {
                                                if (!::cells[mol_index]) continue;
                                                if (::cells[mol_index]->get_type() != CELLS_TYPE_MOLECULES) continue;
                                                t_start = GetTickCount();
                                                count[0] = count[1] = count[2] = num = 0;
                                                atom_index[0] = ::cells[mol_index]->get_index(0);
                                                atom_index[1] = ::cells[mol_index]->get_index(1);
                                                std::string info("");
                                                if (atom_index[0] == atom_index[1]) info = " (single species)";
                                                std::cout << std::endl << "creating " << ::cells[mol_index]->get_species_name() << " molecules" << info << " on " << +num_threads << " threads ... (" << (find_nearest ? "nearest" : "stop at first") << ")" << std::endl;
                                                for (uint8_t i = 0; i < num_threads; ++i) { 
                                                    struct cmd_data *data = new struct cmd_data;
                                                    data->in_index [0]   = atom_index[0]; // cells index i
                                                    data->in_index [1]   = atom_index[1]; // cells index j
                                                    data->out_index[0]   = mol_index; // save molecules back into same ::cells[mol_index]
                                                    data->out_index[1]   = checked_counter; // unique counter for each molecule search
                                                    data->count  [0]     = 0;
                                                    data->count  [1]     = 0;
                                                    data->count  [2]     = 0;
                                                    data->total          = 0;
#ifdef _DEBUG                
                                                    if (to_helper->put(new class queue_entry(THREAD_CREATE_MOL, data)) != QUEUE_STATUS_OK) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: queue 'to_helper' error 'put create MOL'!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -250;
                                                    }
#else
                                                    to_helper->put(new class queue_entry(THREAD_CREATE_MOL, data));
#endif
                                                }
                                                // wait until all threads are finished
                                                for (uint8_t i = 0; i < num_threads; ++i) { 
                                                    entry = from_helper->get(1, WAIT_INFINITE);
                                                    if (entry) {
                                                        if (entry->result) {
                                                            if (!error) error = entry->result;
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: thread " << +entry->id << " create molecule error " << entry->result << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        }
                                                        else if (entry->data) {
                                                            // count molecules and energies
                                                            struct cmd_data *data = reinterpret_cast<struct cmd_data*>(entry->data);
                                                            count[0]                   += data->count[0];
                                                            count[1]                   += data->count[1];
                                                            count[2]                   += data->count[2];
                                                            num                        += data->total;
                                                            ReductionOfEnergy          += data->E[E_RED];
                                                            PotentialEnergyOfMolecules += data->E[E_POT];
                                                            KineticEnergyOfMolecules   += data->E[E_KIN];
                                                            delete data;
                                                        }
                                                        else {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " no molecule data from thread?!" << std::endl << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -251;
                                                        }
                                                        delete entry;
                                                    }
#ifdef _DEBUG
                                                    else if (from_helper->get_status() != QUEUE_STATUS_OK) { // error
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get create MOL'!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -253;
                                                    }
                                                    else { // timeout unexpected!
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get create MOL' timeout!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        if (!error) error = -254;
                                                    }
#else
                                                    else { // there should be no timeout!
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " timeout molecule result?!" << std::endl << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -252;
                                                    }
#endif
                                                }
                                                if (error) break;

                                                // count total number of molecules
                                                mol_N += num;

                                                // get number of molecules
                                                uint64_t tmp = 0;
#ifdef USE_OLD_MOLSEARCH
                                                // old code saves molecules in each pthreads[i]->x[mol_index] entry
                                                for (uint8_t i = 0; i < num_threads; ++i) { 
                                                    tmp += pthreads[i].num[mol_index];
                                                }
#else
                                                // new code saves ::cells[mol_index] with created molecules of all threads
                                                class cell_data *cell = cells[mol_index]->get_first();
                                                while (cell) {
                                                    tmp += cell->species.get_num();
                                                    cell = cell->get_next();
                                                }
#endif                                                     
                                                if (tmp != num) {
                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                    std::cout << std::endl << "error: total number of molecules not matching! " << tmp << " != " << num << std::endl << std::endl;
                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    error = -255;
                                                    break;
                                                }
                                                
                                                if (count[0] || count[1] || count[2]) {
                                                    std::cout << "info: pairing failed counters = {" << count[0] << " " << count[1] << " " << count[2] << "} (ok)"<< std::endl;
                                                }

                                                // calculate energies per molecule
                                                if (num > 0) {
                                                    ReductionOfEnergy          /= num;
                                                    PotentialEnergyOfMolecules /= num;
                                                    KineticEnergyOfMolecules   /= num;
                                                }
                                                
                                                REAL_TYPE mol_Etot = KineticEnergyOfMolecules+PotentialEnergyOfMolecules; // nK
                                                mol_E += mol_Etot*num; // total energy in nK used to calculate mol_T in summary
                                                REAL_TYPE mol_T = mol_Etot/3.0; // note: use Ekin = Etot/2 and factor 3/2 from Ekin = 3/2*kB*T.
                                                REAL_TYPE mol_Tcrit, mol_lambda_dB, mol_n3d, mol_psd;
                                                REAL_TYPE mol_R[3] = {0.0, 0.0, 0.0};
                                                
                                                if ((mol_stat_index == STAT_BE) || (mol_stat_index == STAT_BE_NC)) {
                                                    mol_Tcrit=std::pow(num/zeta_3,1.0/3.0)*mol_who*k_nK; // k_nK = (hbar*1e9)/kB; omega -> nK

                                                    // calculattion of phase-space density
                                                    // TODO: this is the classical/Maxwell-Boltzmann limit! see Ketterle-Zwierlein "Making, probing and understanding ultracold Fermi gases"
                                                    //       for Bosons/Fermions we need to insert n3d_mol = +/-Li(3/2,+/-exp((mu-V)*beta))/lambda_dB^3
                                                    //       Li(3/2) can be evaluated numerically.
                                                    //       for T=Tc: n_crit = Li(3/2,1)/lambda_dB^3 = 2.612/lambda_dB^3
                                                    // TODO: I think factor 2 is not correct here since Ekin = Epot = 1/2*kB*T per degree of freedom? 
                                                    mol_lambda_dB = (num > 0) ? std::sqrt(2.0*Pi/(mol_mass*mol_T)*k_nK*k_nK/k_Epot) : 0.0; // in mu
                                                    mol_R[0]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wx*mol_wx*k_Epot)); // in mu; k_Epot = (mn/kB)*1e-3; 1e6*std::sqrt(1e-9) = std::sqrt(1e3);
                                                    mol_R[1]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wy*mol_wy*k_Epot));
                                                    mol_R[2]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wz*mol_wz*k_Epot));
                                                    mol_n3d       = (num >0) ? num/(std::pow(Pi,1.5)*mol_R[0]*mol_R[1]*mol_R[2]) : 0.0; // 1/mu^3
                                                    mol_psd       = (mol_lambda_dB*mol_lambda_dB*mol_lambda_dB)*mol_n3d; // in units of 1
                                                }
                                                else if (mol_stat_index == STAT_FD) {
                                                    mol_Tcrit     = std::pow(6*num,1.0/3.0)*mol_who*k_nK; // TFermi in nK
                                                    mol_lambda_dB = (num > 0) ? std::sqrt(2.0*Pi/(mol_mass*mol_T)*k_nK*k_nK/k_Epot) : 0.0; // in mu
                                                    mol_R[0]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wx*mol_wx*k_Epot)); // mu
                                                    mol_R[1]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wy*mol_wy*k_Epot));
                                                    mol_R[2]      = std::sqrt(2.0*mol_T/(mol_mass*mol_wz*mol_wz*k_Epot));
                                                    mol_n3d       = 0.0; //TODO: n0 = -Li(3/2, -exp(mu/T))/lambda_dB^3
                                                    mol_psd       = (mol_lambda_dB*mol_lambda_dB*mol_lambda_dB)*mol_n3d; // in units of 1
                                                }
                                                else if (mol_stat_index == STAT_MB) {
                                                    mol_Tcrit     = 0.0;
                                                    mol_lambda_dB = (num > 0) ? std::sqrt(2.0*Pi/(mol_mass*mol_T)*k_nK*k_nK/k_Epot) : 0.0; // in mu
                                                    mol_R[0]      = std::sqrt(mol_T/(mol_mass*mol_wx*mol_wx*k_Epot)); // Gauss sigma in mu
                                                    mol_R[1]      = std::sqrt(mol_T/(mol_mass*mol_wy*mol_wy*k_Epot));
                                                    mol_R[2]      = std::sqrt(mol_T/(mol_mass*mol_wz*mol_wz*k_Epot));
                                                    mol_n3d       = (num > 0) ? num/(std::pow(2.0*Pi,1.5)*mol_R[0]*mol_R[1]*mol_R[2]) : 0.0; // 1/mu^3
                                                    mol_psd       = (mol_lambda_dB*mol_lambda_dB*mol_lambda_dB)*mol_n3d; // in units of 1
                                                }
                                                else {
                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                    std::cout << "error: molecule statistics " << stat_all[mol_stat_index] << " not implemented!" << std::endl;
                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                    error = -257;
                                                    break;
                                                }

                                                if (file_histogram.length() > 0) {
                                                    // TODO: will not work for old code since does not save molecules into cells.
                                                    std::string file_name = expand_filename(file_histogram, ::cells[mol_index]->get_species_name(), n_rep);
                                                    Histogram(file_name, ::cells[mol_index]->get_species_name(), mol_stat_use, 
                                                            n_rep, //rep, 
                                                            num_bins, //bins, 
                                                            mol_index,
                                                            mol_T, mol_Tcrit,
                                                            mol_wx, mol_wy, mol_wz,
                                                            mol_mass,
                                                            false
                                                            );
                                                }

                                                // molecule creation time in seconds
                                                REAL_TYPE dt = ((REAL_TYPE)get_ticks_delta(t_start))/1e3;
                                                mol_t += dt;
                                                
                                                REAL_TYPE mol_eff;
                                                {
                                                    uint64_t N0 = ::cells[atom_index[0]]->get_num_atoms();
                                                    uint64_t N1 = ::cells[atom_index[1]]->get_num_atoms();
                                                    if (atom_index[0] == atom_index[1]) {
                                                        if (::atom_num != 1) {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: unexpected same cell list index " << atom_index[0] << " for multi-species molecule?" << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -258;
                                                            break;
                                                        }
                                                        mol_eff = mol_efficiency(num, &N0, 0)*100.0;
                                                    }
                                                    else {
                                                        mol_eff = mol_efficiency(num, &N0, N1)*100.0;
                                                    } 
                                                    if (cells_count[2] > 1) {
                                                        // output of conversion efficiency for each molecule type, otherwise just output total below.
                                                        std::cout << "created " << num << " " << ::cells[mol_index]->get_species_name() << " molecules" << info << 
                                                            ": conversion efficiency = " << num << " / " << N0 << " = " << mol_eff << 
                                                            " % (" << dt << " s)" << std::endl;
                                                    }
                                                }
                                                out     << "molecule                   = " << ::cells[mol_index]->get_species_name() + info    << std::endl;
                                                out     << "statistics                 = " << stat_all[::cells[mol_index]->get_statistics()]   << std::endl;
                                                out     << "calculation time           = " << dt                           << " s"     << std::endl;
                                                out     << "N_mol                      = " << num                                      << std::endl;
                                                out     << "conv. efficiency           = " << mol_eff                      << " %"     << std::endl;
                                                
                                                //out     << "ReductionOfEnergy          = " << ReductionOfEnergy          << " nK" << std::endl; // Andi: fluctuates a lot from shot-shot!
                                                out     << "Epot_mol                   = " << PotentialEnergyOfMolecules   << " nK"    << std::endl;
                                                out     << "Ekin_mol                   = " << KineticEnergyOfMolecules     << " nK"    << std::endl;
                                                out     << "Etot_mol                   = " << mol_Etot                     << " nK"    << std::endl;
                                                out     << "T_mol                      = " << mol_T                        << " nK"    << std::endl;
                                                if ((mol_stat_index == STAT_BE) || (mol_stat_index == STAT_BE_NC)) {
                                                    out << "Rx_mol                     = " << mol_R[0]                     << " um"    << std::endl;
                                                    out << "Ry_mol                     = " << mol_R[1]                     << " um"    << std::endl;
                                                    out << "Rz_mol                     = " << mol_R[2]                     << " um"    << std::endl;
                                                    out << "Tc_mol                     = " << mol_Tcrit                    << " nK"    << std::endl;
                                                    out << "T_mol/Tc_mol               = " << ((num > 0) ? mol_T/mol_Tcrit : 0.0)      << std::endl;
                                                    out << "lambda_dB                  = " << mol_lambda_dB                << " um"    << std::endl;
                                                    out << "n3d_mol                    = " << mol_n3d                      << " 1/um^3"<< std::endl;
                                                    out << "psd_mol                    = " << mol_psd                                  << std::endl;
                                                }
                                                else if (mol_stat_index == STAT_FD) {
                                                    out << "Rx_mol                     = " << mol_R[0]                     << " um"    << std::endl;
                                                    out << "Ry_mol                     = " << mol_R[1]                     << " um"    << std::endl;
                                                    out << "Rz_mol                     = " << mol_R[2]                     << " um"    << std::endl;
                                                    out << "TF_mol                     = " << mol_Tcrit                    << " nK"    << std::endl;
                                                    out << "T_mol/TF_mol               = " << ((num > 0) ? mol_T/mol_Tcrit : 0.0)      << std::endl;
                                                }
                                                else if (mol_stat_index == STAT_MB) {
                                                    out << "sigma_x_mol                = " << mol_R[0]                     << " um"    << std::endl;
                                                    out << "sigma_y_mol                = " << mol_R[1]                     << " um"    << std::endl;
                                                    out << "sigma_z_mol                = " << mol_R[2]                     << " um"    << std::endl;
                                                    out << "n3d_mol                    = " << mol_n3d                      << " 1/um^3"<< std::endl;
                                                }
                                                out     << "counters                   = [" << count[0] << "," << count[1] << "," << count[2] << "]" << std::endl;
                                                out     << std::endl;

                                                // increase molecule search counter = number of molecule cells
                                                ++checked_counter;
                                            
                                            } // next cell (molecule)
                                        
                                            if (!error) {
                                                // check atom status consistency of created molecules
                                                // all first species atoms must be paired or unpaired, all second species must be free or paired
                                                // and #paired must be consistent with #molecules
                                                // TODO: in final code could be moved to _DEBUG but so far would keep it active always.
                                                uint64_t free        [NUM_ATOMS_LISTS] = {0},
                                                         paired      [NUM_ATOMS_LISTS] = {0}, 
                                                         unpaired    [NUM_ATOMS_LISTS] = {0},
                                                         atoms_paired[2              ] = {0},
                                                         atoms_tot   [2              ] = {0};
                                                uint8_t  cell_ok     [NUM_ATOMS_LISTS] = {0},
                                                         cells_ok    [3              ] = {0};
                                                for (cell_index = 0; cell_index < NUM_CELLS_LISTS; ++cell_index) {
                                                    if (!::cells[cell_index]) continue;
                                                    std::string name     = ::cells[cell_index]->get_species_name();
                                                    uint64_t num_entries = ::cells[cell_index]->get_num_atoms();
                                                    uint64_t count       = 0;
                                                    if (::cells[cell_index]->get_type() == CELLS_TYPE_ATOMS) {
                                                        // atom
                                                        // count free, paired and unpaired atoms
                                                        
                                                        if (cell_index >= NUM_ATOMS_LISTS) {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: " << name << " list index " << +cell_index << " >= " << NUM_ATOMS_LISTS << " unexpected!" << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -260;
                                                            break;
                                                        }
#if NUM_CELLS > 1
                                                        class cell_data *cell = ::cells[cell_index]->get_first();
                                                        for (; cell && (!error); cell = cell->get_next()) { 
                                                            class species_entry *entry = cell->species.get_first(); 
                                                            for (; entry != NULL; ++count, entry = entry->get_next()) {
                                                                uint8_t status = entry->species.status;
                                                                if      (IS_STATUS_FREE    (status)) ++free    [cell_index];
                                                                else if (IS_STATUS_PAIRED  (status)) ++paired  [cell_index];
                                                                else if (IS_STATUS_UNPAIRED(status)) ++unpaired[cell_index];
                                                            }
                                                        }
#else
                                                        for (uint8_t i = 0; i < num_threads; ++i) {
                                                            uint8_t *pstatus = pthreads[i].status[cell_index]; 
                                                            uint64_t num     = pthreads[i].num   [cell_index]; 
                                                            for (uint64_t atom = 0; atom < num; ++atom, ++pstatus) {
                                                                uint8_t status = (*pstatus);
                                                                if      (IS_STATUS_FREE    (status)) ++free    [cell_index];                                                
                                                                else if (IS_STATUS_PAIRED  (status)) ++paired  [cell_index];
                                                                else if (IS_STATUS_UNPAIRED(status)) ++unpaired[cell_index];
                                                             }
                                                             count += num;
                                                        }
#endif
                                                        if (count != num_entries) {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: " << name << " #atoms " << num_entries << " != " << count << " found!" << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -261;
                                                            break;
                                                        }
                                                        else if (count != (free[cell_index]+paired[cell_index]+unpaired[cell_index])) {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: " << name << " #atoms " << count << " != free + paired + unpaired = " << free[cell_index] << " + " << paired[cell_index] << " + " << unpaired[cell_index] << " = " << (free[cell_index]+paired[cell_index]+unpaired[cell_index]) << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -262;
                                                            break;
                                                        }
                                                    }
                                                    else if (::cells[cell_index]->get_type() == CELLS_TYPE_MOLECULES) {
                                                        // molecule
                                                        // sum up paired atom0 and atom1 individually
                                                        // note: atom cells are always before molecule cells, so we can just use atom counters.
                                                        atom_index[0] = ::cells[cell_index]->get_index(0);
                                                        atom_index[1] = ::cells[cell_index]->get_index(1);
                                                        for (uint8_t i = 0; i < ::atom_num; ++i) {
                                                            if (cell_ok[atom_index[i]] == 0) { // ensure we count each cell only once
                                                                atoms_tot[i] += free[atom_index[i]] + paired[atom_index[i]] + unpaired[atom_index[i]];
                                                                if (atom_index[0] == atom_index[1]) // single species: pairs = 2*mol_N
                                                                    atoms_paired[i] += paired[atom_index[i]]>>1;
                                                                else
                                                                    atoms_paired[i] += paired[atom_index[i]];
                                                                // check free/unpaired which must be 0 for atom 0/1
                                                                if ((i == 0) && free[atom_index[i]]) {
                                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                                    std::cout << "error: " << ::cells[atom_index[i]]->get_species_name() << " atom (0) has " << free[atom_index[i]] << " free atoms! this is a bug." << std::endl;
                                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                                    error = -263;
                                                                    break;
                                                                }
                                                                else if ((i == 1) && unpaired[atom_index[i]]) {
                                                                    SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                                    std::cout << "error: " << ::cells[atom_index[i]]->get_species_name() << " atom (1) has " << unpaired[atom_index[i]] << " unpaired atoms! this is a bug." << std::endl;
                                                                    RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                                    error = -264;
                                                                    break;
                                                                }
                                                                cell_ok[atom_index[i]] = 1;
                                                                ++cells_ok[i]; 
                                                            }
                                                        }
                                                        if(error) break;
                                                        
                                                        // count molecules
#if NUM_CELL > 1
                                                        class cell_data *cell = ::cells[cell_index]->get_first();
                                                        for (; cell && (!error); cell = cell->get_next()) { 
                                                            class species_entry *entry = cell->species.get_first(); 
                                                            for (; entry != NULL; ++count, entry = entry->get_next());
                                                        }
#else
                                                        for (uint8_t i = 0; i < num_threads; ++i) {
                                                            count += pthreads[i].num[cell_index];
                                                        }
#endif
#ifndef USE_OLD_MOLSEARCH
                                                        if (count != num_entries) {
                                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                            std::cout << "error: " << name << " #molecules " << num_entries << " != " << count << " found!" << std::endl;
                                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                            error = -265;
                                                            break;
                                                        }
#endif                                
                                                        ++cells_ok[2];
                                                    }
                                                } // next cell
                                                
                                                if (!error) {
                                                    if (cells_ok[0] != cells_count[0]) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[0] << " atom (0) cell count " << cells_ok[0] << " != " << cells_count[0] << " !" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -266;
                                                    }
                                                    else if (cells_ok[1] != cells_count[1]) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[1] << " atom (1) cell count " << cells_ok[1] << " != " << cells_count[1] << " !" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -267;
                                                    }
                                                    else if (cells_ok[2] != cells_count[2]) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[0] << " molecule cell count " << cells_ok[2] << " != " << cells_count[2] << " !" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -268;
                                                    }
                                                    else if (atoms_tot[0] != ::atom_N[0]) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[0] << " #atoms (0) " << atoms_tot[0] << " != " << ::atom_N[0] << " found!" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -269;
                                                    }
                                                    else if (atoms_paired[0] != mol_N) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[0] << " #paired atoms (0) " << atoms_paired[0] << " != molecule number = " << mol_N << " !" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -270;
                                                    }
                                                    else if ((::atom_num == 2) && (atoms_tot[1] != ::atom_N[1])) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[1] << " #atoms (1) " << atoms_tot[1] << " != " << ::atom_N[1] << " found!" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -271;
                                                    }
                                                    else if ((::atom_num == 2) && (atoms_paired[1] != mol_N)) {
                                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                        std::cout << "error: " << ::atom_name[1] << " #paired atoms (1) " << atoms_paired[1] << " != molecule number = " << mol_N << " !" << std::endl;
                                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                        error = -272;
                                                    }
                                                    else std::cout << ((::atom_num == 1) ? ::atom_name[0] : (::atom_name[0]+::atom_name[1])) << " check atoms status field vs. number of molecules ok" << std::endl;
                                                }
                                            } // end check atom status

                                            // checked_counter must be the same as the cells_count for molecules
                                            if ((!error) && (checked_counter != cells_count[2])) {
                                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                                std::cout << "error: inconsistent number of molecule cells " << +checked_counter << " != " << +cells_count[2] << "!" << std::endl;
                                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                error = -276;  
                                            }
                                            else if (!error) {

                                                // summary molecule
                                                // note: efficiency = number of molecules / maximum possible number of molecules.
                                                //       for one species: efficiency = mol_N/(N0/2), N1 = 0 
                                                //       for two species: efficiency = mol_N/(min(N0, N1))
                                                uint64_t N0 = atom_N[0];
                                                REAL_TYPE mol_T = (mol_N > 0) ? mol_E/(3.0*mol_N) : 0.0; // see calculation above
                                                REAL_TYPE mol_eff = mol_efficiency(mol_N, &N0, atom_N[1])*100.0;
                                                SET_CONSOLE_COLOR_INFO(STDOUT_HANDLE);
                                                std::cout << std::endl << "created " << mol_N << " " << mol_name << " molecules total: conversion efficiency = " << 
                                                        mol_N << " / " << N0 << " = " << mol_eff << " % (" << mol_t << " s)" << std::endl << std::endl;
                                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                                out << "molecule                   = " << mol_name << " [summary]"                 << std::endl;
                                                out << "statistics                 = " << mol_stat_use                             << std::endl;
                                                out << "calculation time           = " << mol_t                        << " s"     << std::endl;
                                                out << "mass                       = " << mol_mass                     << " amu"   << std::endl;
                                                out << "N_mol                      = " << mol_N                                    << std::endl;
                                                out << "T_mol                      = " << mol_T                        << " nK"    << std::endl;
                                                out << "fx_mol                     = " << mol_wx/(2*Pi)                << " Hz"    << std::endl;
                                                out << "fy_mol                     = " << mol_wy/(2*Pi)                << " Hz"    << std::endl;
                                                out << "fz_mol                     = " << mol_wz/(2*Pi)                << " Hz"    << std::endl;
                                                out << "conv. efficiency           = " << mol_eff                      << " %"     << std::endl;

                                                // molecule conversion efficiency sum for mean and standard deviation
                                                sum     += mol_eff;
                                                sum_sqr += mol_eff * mol_eff;

                                                // export atoms, molecules and neighbors
                                                if ((export_atom[0].length() > 0) || 
                                                    (export_atom[1].length() > 0) || 
                                                    (export_mol.length()     > 0) ||
                                                    (neighbors_file.length() > 0) ) {
                                                    for (cell_index = 0; cell_index < NUM_CELLS_LISTS; ++cell_index) {
                                                        if (!::cells[cell_index]) continue;
                                                        if (::cells[cell_index]->get_type() == CELLS_TYPE_ATOMS) {
                                                            // export atoms
                                                            uint8_t species = ::cells[cell_index]->get_index(0);
                                                            if (export_atom[species].length() > 0) {
                                                                std::string species_name = ::cells[cell_index]->get_species_name();
                                                                std::string file_name = expand_filename(export_atom[species], species_name, n_rep);
                                                                ExportAtoms(file_name, cell_index);
                                                            }
                                                        }
                                                        else if (::cells[cell_index]->get_type() == CELLS_TYPE_MOLECULES) {
                                                            // export molecules
                                                            if (export_mol.length() > 0) {
                                                                std::string file_name = expand_filename(export_mol, ::cells[cell_index]->get_species_name(), n_rep);
                                                                ExportAtoms(file_name, cell_index);
                                                            }
                                                            /* export neighbors
                                                            // note: this uses molecule species. 
                                                            //       function does not use cells!
                                                            if (neighbors_file.length() > 0) {
                                                                std::string name = ::cells[cell_index]->get_species_name();
                                                                std::string file_name = expand_filename(neighbors_file, name, n);
                                                                error = ExportNeighbors(
                                                                            file_name, 
                                                                            ::cells[cell_index]->get_index(0),
                                                                            ::cells[cell_index]->get_index(1), 
                                                                            neighbors_num,
                                                                            mol_distance_index);
                                                                if (error) break;
                                                                name += "_inverted";
                                                                file_name = expand_filename(neighbors_file, name, n);
                                                                error = ExportNeighbors(
                                                                            file_name, 
                                                                            ::cells[cell_index]->get_index(1),
                                                                            ::cells[cell_index]->get_index(0), 
                                                                            neighbors_num,
                                                                            mol_distance_index);
                                                                if (error) break;
                                                            }*/
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        // reset generated species for next repetition
                                        // update: this never deletes atoms
                                        for (uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                                            if (::cells[i]) {
#ifdef _DEBUG
                                                std::cout << "delete species list " << +i << std::endl;
#endif
                                                ::cells[i]->reset(false);
                                            } 
                                        }
                                        
                                        // delete generated atoms and molecules from all threads
                                        // note: on error this is not done!
                                        for (uint8_t i = 0; i < num_threads; ++i) { 
                                            for (uint8_t j = 0; j < NUM_CELLS_LISTS; ++j) {
                                                if (pthreads[i].x[j] != nullptr) {
#ifdef _DEBUG
                                                    std::cout << "thread " << +i << " delete species " << +j << " with " << pthreads[i].num[j] << " atoms" << std::endl;
#endif
                                                    for (uint8_t d = 0; d < DIM_PHASESPACE; ++d) {
                                                        delete [] pthreads[i].x[j][d];
                                                    }
                                                    delete [] pthreads[i].x[j];
                                                    pthreads[i].x[j] = nullptr;
                                                }
                                                if (pthreads[i].E[j] != nullptr) {
                                                    for (uint8_t d = 0; d < DIM_ENERGY; ++d) {
                                                        delete [] pthreads[i].E[j][d];
                                                    }
                                                    delete [] pthreads[i].E[j];
                                                    pthreads[i].E[j] = nullptr;
                                                }
                                                if (pthreads[i].status[j] != nullptr) {
                                                    delete [] pthreads[i].status[j];
                                                    pthreads[i].status[j] = nullptr;
                                                }
                                                pthreads[i].num [j] = 0;
                                                pthreads[i].mass[j] = 0.0;                                
                                            } 
                                        }

                                        REAL_TYPE loop_time = ((REAL_TYPE)get_ticks_delta(t_loop))/1e3;
                                        ss_time = get_end_time(tc_start, tc_duration);
                                        out << std::endl << "calculation time [total]   = " << loop_time << " s" << std::endl;
                                        out << ss_time.str() << " (" << tc_duration.count() << " seconds)"  << std::endl << std::endl;
                                        
                                        //std::cout << "////////////////////////////////////////////////////////////////////////////////" << std::endl;                        
                                        std::cout << "repetition " << n_rep << " / " << repetitions << " calculation time " << loop_time << " s" << std::endl;
                                        std::cout << ss_time.str() << " (" << tc_duration.count() << " seconds)" << std::endl;

                                        t_sum     += loop_time;
                                        t_sum_sqr += loop_time * loop_time;

                                    } // next repetition

                                    // repetitions finished (or error).                    
                                    REAL_TYPE rep_time = ((REAL_TYPE)get_ticks_delta(t_total))/1e3;
                                    ss_time = get_end_time(tc_total, tc_duration);
                                    
                                    if (error) {
                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                        std::cout << std::endl << n_rep << " repetitions finished in " << rep_time << " s with error " << error << std::endl;
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                        out       << "terminated with error      = " << error << std::endl;
                                    }
                                    else {
                                        // molecule conversion efficiency mean and standard deviation
                                        REAL_TYPE mean  = ((REAL_TYPE)sum) / repetitions;
                                        REAL_TYPE stdev = (repetitions == 1) ? 0.0 : std::sqrt((((REAL_TYPE)sum_sqr) - repetitions*mean*mean)/(repetitions-1));
                                        
                                        // calculation time mean and standard deviation
                                        REAL_TYPE t_mean  = ((REAL_TYPE)t_sum) / repetitions;
                                        REAL_TYPE t_stdev = (repetitions == 1) ? 0.0 : std::sqrt((((REAL_TYPE)t_sum_sqr) - repetitions*t_mean*t_mean)/(repetitions-1));

                                        SET_CONSOLE_COLOR_INFO(STDOUT_HANDLE);
                                        std::cout << ((t_stdev < 0.1) ? std::setprecision(3) : std::setprecision(1));
                                        std::cout << std::endl << repetitions << " repetitions finished in " << rep_time << " s (" << t_mean << " +/- " << t_stdev <<" s/repetition)" << std::endl;
                                        std::cout << "conversion efficiency = " << std::setprecision(3) << mean << " +/- " << stdev << "%" << std::setprecision(1) << std::endl << std::endl;
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);

                                        out       << "repetition       [summary] = " << repetitions << " / " << repetitions << std::endl;
                                        out       << "conv. efficiency [summary] = " << mean   << " +/- " << stdev   << " %" << std::endl;
                                        out       << "calculation time [summary] = " << t_mean << " +/- " << t_stdev << " s" << std::endl;
                                    }
                                    
                                    // total time before file is closed
                                    std::cout << ss_time.str() << " (" << tc_duration.count() << " seconds)" << std::endl;
                                    out       << ss_time.str() << " (" << tc_duration.count() << " seconds)" << std::endl << std::endl;
                                                        
                                    out.close();
                                } // output file
                            
                            }// threads started without error
                            
                            // shutdown threads
                            // we send num_threads shutdown commands and wait until all threads respond
                            std::cout << std::endl;
                            for (uint8_t i = 0; i < num_threads; ++i) { 
                                if (pthreads[i].handle != INVALID_THREAD) {
#ifdef _DEBUG            
                                    if (to_helper->put(new class queue_entry(THREAD_SHUTDOWN, NULL)) != QUEUE_STATUS_OK) {
                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                        std::cout << std::endl << "error: queue 'to_helper' error 'put'!" << std::endl << std::endl;
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                        if (!error) error = -280;
                                    }
#else
                                    to_helper->put(new class queue_entry(THREAD_SHUTDOWN, NULL));
#endif
                                    entry = from_helper->get(1, 2*THREAD_TIMEOUT);
                                    if (entry) {
                                        if (entry->result) {
                                            SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                            std::cout << "error: shutdown thread " << +entry->id << "/" << +num_threads << " error " << entry->result << std::endl; 
                                            RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                            if (!error) error = entry->result;
                                        }
                                        else {
                                            //std::cout << "shutdown thread " << +entry->id << "/" << +num_threads << " ..." << std::endl; 
                                        }
                                        delete entry;
                                    }
#ifdef _DEBUG
                                    else if (from_helper->get_status() != QUEUE_STATUS_OK) { // error
                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                        std::cout << std::endl << "error: thread " << +i << "/" << +num_threads << " queue 'from_helper' error 'get SHUTDOWN'!" << std::endl << std::endl;
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                        if (!error) error = -281;
                                    }
#endif
                                    else {
                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                        std::cout << "shutdown thread " << +i << "/" << +num_threads << " no responds!" << std::endl; 
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                        if (!error) error = -282;
                                    }
                                }
                            }

                            // now all threads should have terminated
                            for (uint8_t i = 0; i < num_threads; ++i) { 
                                if (pthreads[i].handle != INVALID_THREAD) {
                                    int tmp = thread_shutdown(pthreads[i].handle, 2*THREAD_TIMEOUT);
                                    if (tmp) {
                                        if (!error) error = tmp;
                                        SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                        std::cout << "error: shutdown thread " << +i << "/" << +num_threads << " timeout!" << std::endl; 
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                    }
                                    else {
                                        std::cout << "shutdown thread " << +i << "/" << +num_threads << " ok" << std::endl; 
                                    }
                                }
                            }
                            //std::cout << "shutdown " << +num_threads << " threads ok" << std::endl; 
                            
                            // release buffers for each helper threads
                            for (uint8_t i = 0; i < num_threads; ++i) {
                                /* atom 0, atom 1, molecule 
                                for (int j = 0; j < MAX_NUM_SPECIES; ++j) {
                                    if ( (::atom_num == 1) && (j == 1) ) continue;
                                    // dimension
                                    for (int k = 0; k < DIM_PHASESPACE; ++k) {
                                        delete [] pthreads[i].x[j][k];
                                    }
                                    delete [] pthreads[i].x[j];
                                    if (j < (MAX_NUM_SPECIES-1)) {
                                        delete [] pthreads[i].status[j];
                                    }
                                }*/
                                // delete random number generator
                                if ( pthreads[i].gen_uniform ) delete pthreads[i].gen_uniform;
                                if ( pthreads[i].gen_normal  ) delete pthreads[i].gen_normal;
                                
                                // ensure atoms are already deleted
                                // TODO: on error this is not the case!
                                for (uint8_t j = 0; j < NUM_ATOMS_LISTS; ++j) {
                                    if ((pthreads[i].x[j] != nullptr) || (pthreads[i].status[j] != nullptr) || (pthreads[i].num[j] != 0)) {
                                        SET_CONSOLE_COLOR_WARN(STDOUT_HANDLE);
                                        std::cout << "warning: thread " << +i << " species " << +j << " not deleted!" << std::endl;
                                        RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                    }
                                } 
                                
                                // destroy mutex
                                if (!MUTEX_DELETE_AND_OK(pthreads[i].mutex)) {
                                    if (!error) error = -290;
                                }
                            }

#ifdef _DEBUG
                            if (to_helper->get_status() != QUEUE_STATUS_OK) {
                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                std::cout << std::endl << "error: queue 'to helper' in error state!" << std::endl << std::endl;
                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                if (!error) error = -291;
                            }
                            if (from_helper->get_status() != QUEUE_STATUS_OK) {
                                SET_CONSOLE_COLOR_ERROR(STDOUT_HANDLE);
                                std::cout << std::endl << "error: queue 'from helper' in error state!" << std::endl << std::endl;
                                RESET_CONSOLE_COLOR(STDOUT_HANDLE);
                                if (!error) error = -292;
                            }
#endif
                            // delete threads data and queues            
                            delete [] pthreads;
                            delete to_helper;
                            delete from_helper;

                            // destroy global mutex and barrier
                            if (!MUTEX_DELETE_AND_OK(::global_mutex)) {
                                if (!error) error = -292;
                            }
                            if (!BARRIER_DELETE_AND_OK(::global_barrier)) {
                                if (!error) error = -293;
                            }

                            // delete all species and cells
                            for (uint8_t i = 0; i < NUM_CELLS_LISTS; ++i) {
                                if (cells[i]) {
#ifdef _DEBUG
                                    std::cout << "delete cell " << +i << std::endl;
#endif
                                    delete cells[i];
                                }
                            }
                            
                        }

                    }
                }
            }
        }
    }

    if (seed_generator) delete seed_generator;
        
    if (error) {
        SET_CONSOLE_COLOR_ERROR(STDERR_HANDLE);
        std::cerr << std::endl << "done with error " << error << " (total " << get_ticks_delta(t_total)/1e3 << " s)" << std::endl << std::endl;
        RESET_CONSOLE_COLOR(STDERR_HANDLE);
    }
    else {
        std::cerr << std::endl << "done ok (total " << get_ticks_delta(t_total)/1e3 << " s)" << std::endl << std::endl;
    }
    
    return error;
}



