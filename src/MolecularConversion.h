// cells and specied definition file
// last change 22/09/2026 by Andi

#ifndef _CELLS_AND_SPECIES
#define _CELLS_AND_SPECIES

// actual version string
#define VERSION_INFO    "v1.6"

// define for debugging (slows down code!)
//#define _DEBUG

// if defined use v1.2 molecule search for timing comparison
// note: this has not all distance measures implemented.
//#define USE_OLD_MOLSEARCH

// if defined use uint128.h to simulate uint128_t with my_uint128
//#define CLASS_UINT128

// for windows we must simulate uint128
#if defined(_WIN32) || defined(_WIN64) // windows
#ifndef CLASS_UINT128
#define CLASS_UINT128
#endif
#endif

// data type float or double used in all calculations. 
// REAL_DOUBLE = uses double precision, default, recommended.
// REAL_SINGLE = uses single precision, maybe faster and less entropy draining but not good for large atom numbers! use only for testing!
// TODO: maybe add option to generate double-precision random numbers (higher entropy) but use floats for molecule search (faster)
// note:  32bit random number generators are disabled when REAL_DOUBLE is used. these are anyway not recommended.
#define REAL_SINGLE             32
#define REAL_DOUBLE             64
#define REAL_PRECISION          REAL_DOUBLE
#if (REAL_PRECISION == REAL_DOUBLE)
#define REAL_TYPE               double
#define RNG_INT_TYPE            uint64_t
#define RNG_MASK                ((((uint64_t)1<<52)-1)<<12)
#elif (REAL_PRECISION == REAL_SINGLE)
#define REAL_TYPE               float
#define RNG_INT_TYPE            uint32_t
#define RNG_MASK                ((((uint64_t)1<<23)-1)<<9)
#endif

// vector brackets and separator
#define V_OPEN                  '{'
#define V_CLOSE                 '}'
#define V_SEP                   ','

// indicates that following number is hexadecimal
// used for seed values
#define NUM_HEX                 "0x"

// if ::seed_use is empty use std::random_device to generate SEED_VALUES_NUM random seed values for each run (default).
// note: seed_use can be give with the '-s' option and with parameter file "seed" value.
#define SEED_VALUES_NUM         2

// double-linked list
#include "lists.hpp"

// thread specific definitions
#include "threads.h"

// MSVC defines min / max macros which interfer with std::numeric_limits.
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include <limits>             // std::numeric_limits<REAL_TYPE>::max()

// if defined and with options -tI and -ti output first STREAM_DATA_OUT random numbers in parallel on stdout and stderr.
// this is useful to debug on Windows issues with pipe which is not in binary mode.
#define STREAM_DATA_OUT     10

// define if RNG test should use produced doubles for test. this is more strict. (default). if not defined uses integer.
#define TEST_DOUBLE_CONVERSION

// if defined show cell offsets (requires _DEBUG flag)
#define SHOW_CELL_OFFSETS

// number of real-space and momentum space dimentions of atoms
#define DIM_SPACE               3
#define DIM_MOMENTUM            3
#define DIM_CELLS               3 // must be 3, for phase-space we use independent cell_data for position and velocity.
#define DIM_PHASESPACE          (DIM_SPACE+DIM_MOMENTUM)
#define DIM_MAXVAL              (DIM_PHASESPACE+1) // maximum values for phase-space + Etot

// number of neighboring cells
#define CELL_NEIGHBOURS_OFFSET  1                       // nearest neighbors to consider, 1=next nearest, 2=next-next nearest,..
#define CELL_NEIGHBOURS_PER_DIM (CELL_NEIGHBOURS_OFFSET*2+1)
#define CELL_NEIGHBOURS_POW3    (CELL_NEIGHBOURS_PER_DIM*CELL_NEIGHBOURS_PER_DIM*CELL_NEIGHBOURS_PER_DIM)
#define CELL_NEIGHBOURS         (CELL_NEIGHBOURS_POW3-1)
#define CELL_NEIGHBOURS_BEFORE  (CELL_NEIGHBOURS>>1)    // first half are cells created before actual cell
#define CELL_OFFSET             1

// fixed number of cells
//#define NUM_CELLS               125000      // cells_per_dim = 50 in 3d and 7 in 6d (slow)
//#define NUM_CELLS               20000      // cells_per_dim = 27 in 3d and 5 in 6d
//#define NUM_CELLS               64           // enough cells such that each thread has at least 1 (2^6). 
#define NUM_CELLS               1               // >=64 such that each thread has at least 1 cell
                                                // 1 = one cell per thread
                                                // 0 = same as 1 but molecule search bypasses cells (as in old code)


// maximum number of species = 2 atomic species
// note: for Bosons we create separate cells for non-condensed and condensed parts, i.e. with different statistics.
#define MAX_NUM_SPECIES         2

// index into 'cells' for species 0 and 1 and molecules and additional cells.
// molecule might require 2 additional species-selective cells when species have different distance measure than molecule.
// all cells except CELL_MOLECULE are created when atoms are generated.
//#define CELLS_SPECIES_0         0
//#define CELLS_SPECIES_0_BEC     1
//#define CELLS_SPECIES_1         2
//#define CELLS_SPECIES_1_BEC     3
//#define CELLS_MOLECULE          4

// cells type
#define CELLS_TYPE_ATOMS        0
#define CELLS_TYPE_MOLECULES    1
#define CELLS_TYPE_NEAREST      2

// number of cells lists with different distance measure
#define MAX_CELLS_PER_MOLECULE  4   // worst case with 2x BEC gives 2^2=4 combinations: n0+n1,n0+BEC1,BEC0+n1,BEC0+BEC1
#define NUM_ATOMS_LISTS         (MAX_NUM_SPECIES*2)
#define NUM_CELLS_LISTS         (NUM_ATOMS_LISTS + MAX_CELLS_PER_MOLECULE) 

// maximum tolerated number of atoms with E > Emax. 
// for Metropolis gives just a warning for fewer atoms.
// this allows user to increase calculation size next time without loosing the result.
#define EMAX_ERROR_NUM_ATOMS    10

// overall scaling factor for maximum energy and calculation size. adjusted to have calc_size = 1.0.
#define EMAX_SCALING            4.0

// gives a warning when maximum Energy of atoms < EMIN_SCALING * maximum calculation size/energy. 
#define EMIN_SCALING            0.5

// fraction of the maximum possible distance used for min/max values of cells.
// this can be <<1 since it is very unlikley that such a distance is ever measured.
// cells_list::get_hash_key is adding outside distances to first or last cell.
//#define CELL_MAX_DISTANCE_FRACTION  0.1

// no species for 2nd index of CELLS_TYPE_ATOMS
#define SPECIES_NONE            0xff

// energy coordinates
#define DIM_ENERGY              3
#define E_POT                   0       // potential energy
#define E_KIN                   1       // kinetic energy
#define E_RED                   2       // energy gain
#define E_MAX                   2       // total maximum energy of atoms

// atom status
#ifdef __DO_NOT_USE__
#define STATUS_PAIRED_MASK                                  0x01
#define STATUS_UNPAIRED_MASK                                0xfe
#define GET_CHECKED_CODE(status)                            ((status)  & STATUS_UNPAIRED_MASK)
//#define _STATUS_FREE(status)                                (((status) & STATUS_PAIRED_MASK)==0)
//#define _STATUS_PAIRED(status)                              (((status) & STATUS_PAIRED_MASK)!=0)
//#define _STATUS_CHECKED(status, checked_code)               (GET_CHECKED_CODE(status)==(checked_code))
#define _STATUS_FREE(status)                                ((status & STATUS_PAIRED_MASK)==0)
#define _STATUS_PAIRED(status)                              ((status & STATUS_PAIRED_MASK)!=0)
// start molecule search: atom is free and not checked with current checked_code:
// TODO: something is not working when I use this? at the moment use it only for testing.
#ifdef _DEBUG
//#define STATUS_FREE(status, checked_code)                   (_STATUS_FREE  (status) && (!_STATUS_CHECKED(status,checked_code)))
#define STATUS_FREE(status, checked_code)                   ((status & STATUS_PAIRED_MASK)==0)
#endif
// after molecule search: status is paired or unpaired with current checked_code
//#define STATUS_PAIRED(status, checked_code)                 (_STATUS_PAIRED(status) && ( _STATUS_CHECKED(status,checked_code)))
//#define STATUS_UNPAIRED(status, checked_code)               (_STATUS_FREE  (status) && ( _STATUS_CHECKED(status,checked_code)))
#define STATUS_PAIRED(status, checked_code)                 (status==(checked_code | 1))
#define STATUS_UNPAIRED(status, checked_code)               ((status & STATUS_PAIRED_MASK)==0)
#define SET_STATUS_PAIRED(status, checked_code)             status=((checked_code) | 1)
//#define SET_STATUS_UNPAIRED(status, checked_code)           status=((checked_code) | 0)
//#define SET_STATUS_FREE(status)                             SET_STATUS_UNPAIRED(status, STATUS_UNPAIRED_MASK)
//#define SET_STATUS_PAIRED(status, checked_code)             status=1
#define SET_STATUS_UNPAIRED(status, checked_code)           status=((checked_code) | 0)
#define SET_STATUS_FREE(status)                             status=0
#define MAKE_CHECKED_CODE(counter)                          (((1+(counter))<<1)&(STATUS_UNPAIRED_MASK))

#else

#define STATUS_FREE                                         0
#define STATUS_PAIRED                                       1
#define STATUS_UNPAIRED                                     2

#define IS_STATUS_FREE(status)                              ((status)==STATUS_FREE)
#define IS_STATUS_PAIRED(status)                            ((status)==STATUS_PAIRED)
#define IS_STATUS_UNPAIRED(status)                          ((status)==STATUS_UNPAIRED)

#define SET_STATUS_FREE(status)                             (status=STATUS_FREE)
#define SET_STATUS_PAIRED(status)                           (status=STATUS_PAIRED)
#define SET_STATUS_UNPAIRED(status)                         (status=STATUS_UNPAIRED)

#endif

// raise error codes
#define ERROR_UNKNOWN_DIST      -101    // stat_index[species] or use_stat[species] is unknown

// convert coordinate x in phase-space dimention dim to cell index: (atom->x[dim] - xL[i]) / dx[i]
//#define GET_CELL_INDEX(x, dim)  (uint8_t)( ((x)[dim] + ::max_x[dim] + cell_data::dx[dim]/2.0) / cell_data::dx[dim] )
//#define GET_CELL_INDEX(x, dim)  (uint8_t)( ((x)[dim] + ::max_x[dim] + ::dx[dim]/2.0) / ::dx[dim] )

// conversion factor of distance to uin64_t
//#define DISTANCE_TO_UINT64(d)   (uint64_t)(d*1e12)  // distance resolution 1e-12, maximum 1.8e7

// date/time format
#define FMT_DATE        "%Y-%m-%d %X"

// generic probability function
typedef REAL_TYPE (*prob_func)(REAL_TYPE mass, REAL_TYPE T, REAL_TYPE mu, REAL_TYPE omega[DIM_SPACE], REAL_TYPE x[DIM_PHASESPACE], REAL_TYPE *Epot, REAL_TYPE *Ekin, REAL_TYPE Emin);

// forward declarations
struct thread_data;
class cell_data;

// atoms or molecule properties
class species_data {
public:
    // mass in amu or mn
    REAL_TYPE mass;
    // real-space and momentum space coordinates
    REAL_TYPE x[DIM_PHASESPACE];
    // total energies in nK. use index E_POT, E_KIN, E_RED.
    REAL_TYPE E[DIM_ENERGY];
    // atom cell index for position and momentum
    uint8_t index[DIM_PHASESPACE];
    // status
    uint8_t status;
    // molecule distance^2
    REAL_TYPE distance;
    // unique atom index. for atom only first index is used, for molecule gives both atom indices.  
    // note: get_next returns atoms with unsorted index. this would require one extra walk through list of atoms. 
    uint64_t atom_index[2];

    species_data();
    ~species_data();

};

// list of species. 
// use double_linked_list::add_sorted to create sorted list with increasing hash_key
class species_entry {
private:
    // double-linked list of previous and next atom
    friend double_linked_list<class species_entry, uint64_t>;
    friend double_linked_list<class species_entry, REAL_TYPE>;
    class species_entry *next;
    class species_entry *prev;
    // pointer to one species or a pair of species, first is never NULL.
    // TODO: second species is never used but might be used for FindNearest.
    //class species_data *species[2];
    // distance = hash_key to sort species into double-linked list in cells
    //uint64_t hash_key;
public:    
    class species_data species; // update: single species, no pointer, public
    uint64_t hash_key; // update: public

    //species_entry(species_data *species0, species_data *species1, uint64_t hash_key);
    species_entry();
    ~species_entry();

    // return species
    //class species_data * get_species(uint8_t index) { return species[index]; };
    // return hash_key
    uint64_t get_hash_key(void) { return hash_key; };
    // returns next entry
    class species_entry * get_next(void) { return next; };
    // returns previous entry
    class species_entry * get_prev(void) { return prev; };
    
    /* returns next species entry with distance^2 <= distance_squared starting from entry (inclusive)
    class species_entry * get_next_distance(class species_entry *first, REAL_TYPE *distance_squared, uint8_t distance_index);
    // returns next cell with distance^2 <= distance_squared starting from input cell (inclusive)
    class cell_data * get_next_distance(
        class cell_data *first, 
        REAL_TYPE      *distance_squared, 
        uint8_t          distance_index,
        const REAL_TYPE *dx2,
        REAL_TYPE      mass1
    );*/
};

// double-linked list of cells. 
// first has prev = NULL, last has next = NULL.
// dimensions 'dim' and 'cells_per_dim' are saved in cells_list
class cell_data {
    // double-linked list of previous and next cells
    friend class double_linked_list<class cell_data, uint64_t>;
    friend species_entry;
    class cell_data *next;
    class cell_data *prev;
    // each cell is associated with a thread. the mutex of the thread is used to lock cell.
    struct thread_data *owner_thread;
    // hash key is used to assign atoms to cells = linear cell index
    uint64_t hash_key;
    // 6 or 12 neigboring cells for 3d or 6d cells
    class cell_data ** nb;
//#ifdef _DEBUG   
// TODO: can these be removed?
    // 3 or 6 cell index
    uint8_t * index;
    // 3 or 6 center coordinates of cell
    REAL_TYPE * x;
//#endif
public:
    // double-linked list of species within cell
    class double_linked_list<class species_entry, uint64_t> species;
    
    cell_data(
        struct thread_data * owner, 
        uint8_t dim, 
        uint64_t hash_key,
//#ifdef _DEBUG
        uint8_t *index, 
        REAL_TYPE *x
//#endif
    );
    ~cell_data();

    // returns next cell
    class cell_data * get_next(void) { return next; }
    // returns previous cell
    class cell_data * get_prev(void) { return prev; }

    // returns owner thread and id
    struct thread_data * get_owner(void);
    uint8_t get_owner_id(void);

    // lock and unlock cell
    void lock(void);
    void unlock(void);

    // return hash_key
    uint64_t get_hash_key(void) { return hash_key; }

    // returns absolute index of this cell
    //inline uint64_t get_index(void);
    // returns absolute index of given dimensional index
    //static inline uint64_t get_index(uint8_t index[DIM_PHASESPACE], uint8_t offset);

    // returns cell of given relative DIM_CELL-dimensional index starting from this cell
    //class cell_data * find_cell(int16_t _index[DIM_CELLS]);
};

// double-linked list of cells with min/max/delta values and get_hash_key.
class cells_list : public double_linked_list<class cell_data, uint64_t> {
private:
    uint8_t type;           // cell type used for CreateAtoms/Molecules/... to identify which cells to use.
    uint8_t index[2];       // species or cell index depending on type
    uint8_t statistics, distance_index;
    REAL_TYPE gamma;        // scaled gamma^2 for distance threshold
    std::string species_name;
    uint64_t num_atoms;     // atoms to be created
    //bool duplicate;         // if True cell->species contains copy of species which must not be deleted. update: all are duplicate!
    uint8_t dim;            // dimensions 3 or 6
    uint8_t offset;         // offset of ::atom_max phase-space index. 0 except 3 for DISTANCE_V
    uint8_t cells_per_dim;  // cells per dimension. total number of cells = cells_per_dim^dim.
    REAL_TYPE *x_min, *dx;  // 3 or 6-dim first center coordinate and step per cell 
    REAL_TYPE *max, *R;     // calculation size and cloud radius
public:
    cells_list(
        struct thread_data *pthreads,   // each cell gets one owner thread from this list
        uint8_t             type,
        uint64_t            num_cells,  // rounded to next lower cells_per_dim^dim.
        uint64_t            num_atoms,  // total number of atoms in cell
        uint8_t             index0, 
        uint8_t             index1,
        std::string         species_name, 
        uint8_t             statistics,
        uint8_t             distance_index, 
        REAL_TYPE           gamma, 
        //bool                duplicate,
        REAL_TYPE           *max,
        REAL_TYPE           *R
        );
    ~cells_list();
    
    // delete species and if delete_all delete cells
    int reset(bool delete_all); 

    // get private properties
    uint8_t      get_type          (void)       { return type;           };
    uint8_t      get_dim           (void)       { return dim;            };
    uint8_t      get_cells_per_dim (void)       { return cells_per_dim;  };
    uint64_t     get_num_atoms     (void)       { return num_atoms;      };
    uint8_t      get_index         (uint8_t i)  { return index[i];       };
    std::string& get_species_name  (void)       { return species_name;   };
    uint8_t      get_statistics    (void)       { return statistics;     };
    uint8_t      get_distance_index(void)       { return distance_index; };
    REAL_TYPE    get_gamma         (void)       { return gamma;          };
    //bool         get_duplicate     (void)       { return duplicate;      };
    REAL_TYPE *  get_x_min         (void)       { return x_min;          };
    REAL_TYPE *  get_dx            (void)       { return dx;             };
    uint8_t      get_offset        (void)       { return offset;         };
    REAL_TYPE    get_max           (uint8_t i)  { return (max && (i < DIM_MAXVAL)) ? max[i] : 0.0; }; // 3x mu, 3x mm/s, 1x nK
    REAL_TYPE    get_R             (uint8_t i)  { return (R   && (i < DIM_SPACE )) ? R  [i] : 0.0; }; // 3x mu

    // get hash_key for given species = cell index where to insert species with add_sorted()
    uint64_t get_hash_key(class species_data *species);
};

// data given to each thread
// TODO: convert to class and save as double-linked list?
struct thread_data {
    struct thread_data *next;               // pointer to next thread_data structure. wraps around from last to first.
    MUTEX mutex;                            // used to lock cells
    THREAD_HANDLE handle;                   // thread handle
    uint8_t id;                             // index 0..num-1 of thread
    uint8_t num_threads;                    // number of threads
    class queue *to_helper;                 // queue primary -> helper
    class queue *from_helper;               // queue helper -> primary
    class random_generator<REAL_TYPE> *gen_uniform;// random number generator for uniform real values from 0..1
    class random_generator<REAL_TYPE> *gen_normal; // random number generator for normal distributed real values with mean 0 and sigma 1 
#if defined(_WIN32) || defined (_WIN64)
    DWORD   win_id;                         // windows needs also some id
#endif
    uint64_t            num    [NUM_CELLS_LISTS]; // number of atoms/molecules per list
    REAL_TYPE           mass   [NUM_CELLS_LISTS]; // mass of atom/molecule per list
    REAL_TYPE           **x    [NUM_CELLS_LISTS]; // DIM_PHASESPACE coordinates, x[cell index][dim][# atoms].
    REAL_TYPE           **E    [NUM_CELLS_LISTS]; // DIM_ENERGY energies, E[cell index][dim][# atoms].
    uint8_t             *status[NUM_CELLS_LISTS]; // lists of status-values for atoms only. status[cell index][# atom].
};

// data for helper thread for given command put into queue
// TODO: check what information could be better saved into shared cells_list? need only results.
struct cmd_data {
    uint8_t in_index [2];                   // atoms: [species, statistics], molecules: cell index of species 
    uint8_t out_index[2];                   // cell index and checked_code to use for result
    uint64_t total;                         // in/out: num atoms or molecules to be/actual created by thread.
    uint64_t count[3];                      // in: atom_index, out: efficiency, loop or error/debug counter.
    REAL_TYPE E[DIM_ENERGY];                // out: atoms: max energy in nK, molecules: total energies in nK. use index E_POT, E_KIN, E_RED, E_MAX.
};

// number of species
//extern uint8_t atom_num;

// global double-linked list of shared atoms and molecules between threads
//extern class double_linked_list<class species_data> species[MAX_NUM_SPECIES];

// global double-linked list of shared cells containing all atoms and optional different distance measures
//extern class double_linked_list<class cell_data, REAL_TYPE> cells[NUM_CELLS_LISTS];
extern class cells_list *cells[NUM_CELLS_LISTS];

// global double-linked list of shared nearest neighbors found by FindNearest
extern class double_linked_list<class species_entry, uint64_t> nearest;

// functions called by helper thread
extern int CreateAtoms           (struct thread_data *pdata, struct cmd_data *cdata);
extern int CreateAtoms_Metropolis(struct thread_data *pdata, struct cmd_data *cdata);

#ifdef USE_OLD_MOLSEARCH
extern int CreateMolecules_v1_2(struct thread_data *pdata, struct cmd_data *cdata);
#else
extern int CreateMolecules_linear(struct thread_data *pdata, struct cmd_data *cdata);
#endif
// find nearest neighbors between two species using given distance measure and gamma.
//extern int FindNearest           (const struct thread_data *pdata, struct cmd_data *cdata);


#endif // _CELLS_AND_SPECIES
