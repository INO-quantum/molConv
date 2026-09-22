// helper threads header
// created 24/10/2023 by Andi
// last change 20/9/2026 by Andi

#ifndef _HELPER_THREADS
#define _HELPER_THREADS

// number of helper threads
// the given number of LSB bits are used to offset the seed of the RNG generators for different threads.
// note: if you change this you need to change also gray code!
#define THREADS_MAX_SEED_BITS   6
#define THREADS_MAX_NUM         (1<<THREADS_MAX_SEED_BITS)

#include "random.hpp"                   // random number generator (this needs THREADS_MAX_.. definitions)

#if defined(_WIN32) || defined(_WIN64) // windows

#include <windows.h>                    // windows default definitions
#include <stdint.h>
#include <limits.h>                     // numeric_limits class
#include <io.h>                         // needed for _setmode
#include <fcntl.h>                      // needed for _setmode

// threading
#define THREAD_HANDLE                   HANDLE
#define INVALID_THREAD                  NULL
#define ZEROMEMORY(address, size)       memset(address, 0, size)
//#define THREAD_ENTRY                    DWORD WINAPI
#define RETURN_FROM_THREAD(exit_code)   return(exit_code)
#define THREAD_JOIN(handle, res, err)   err = pthread_join(handle, &res);

// synchronization barrier
#define BARRIER                         SYNCHRONIZATION_BARRIER
#define BARRIER_INIT(bar,num)           InitializeSynchronizationBarrier(&bar, num, -1) // returns TRUE if ok
#define BARRIER_INIT_AND_OK(bar,num)    (InitializeSynchronizationBarrier(&bar, num, -1) != 0) // returns true if ok
#define BARRIER_WAIT(bar)               EnterSynchronizationBarrier(&bar, SYNCHRONIZATION_BARRIER_FLAGS_BLOCK_ONLY|SYNCHRONIZATION_BARRIER_FLAGS_NO_DELETE) // returns TRUE for last thread, otherwise FALSE
#define BARRIER_WAIT_AND_OK(bar,ret)    ( ((ret=EnterSynchronizationBarrier(&bar, SYNCHRONIZATION_BARRIER_FLAGS_BLOCK_ONLY|SYNCHRONIZATION_BARRIER_FLAGS_NO_DELETE)) == TRUE) || ( ret == FALSE ) ) // returns true if ok. ret = integer
#define BARRIER_DELETE(bar)             DeleteSynchronizationBarrier(&bar) // returns always TRUE
#define BARRIER_DELETE_AND_OK(bar)      (DeleteSynchronizationBarrier(&bar) != 0) // returns always true

// mutex
#define MUTEX                           HANDLE
#define MUTEX_INIT(mutex)               mutex = CreateMutexA(NULL, FALSE, NULL) // returns NULL on error
#define MUTEX_INIT_AND_OK(mutex)        ((mutex = CreateMutexA(NULL, FALSE, NULL)) != NULL) // returns true if ok
#define MUTEX_LOCK(mutex)               WaitForSingleObject(mutex, INFINITE) // returns WAIT_OBJECT_0 if ok
#define MUTEX_LOCK_AND_OK(mutex)        (WaitForSingleObject(mutex, INFINITE) == WAIT_OBJECT_0) // returns true if ok
#define MUTEX_UNLOCK(mutex)             ReleaseMutex(mutex) // returns nonzero if ok
#define MUTEX_UNLOCK_AND_OK(mutex)      (ReleaseMutex(mutex) != 0) // returns true if ok
#define MUTEX_DELETE(mutex)             CloseHandle(mutex) // returns nonzero on success
#define MUTEX_DELETE_AND_OK(mutex)      (CloseHandle(mutex) != 0) // returns true if ok

// semaphore
#define SEMAPHORE                       HANDLE
#define SEMAPHORE_INIT(sem, init)       sem = CreateSemaphore(NULL, init, 0x7fffffff, NULL) // returns NULL on error
#define SEMAPHORE_INIT_AND_OK(sem, init)    ( (sem = CreateSemaphore(NULL, init, 0x7fffffff, NULL)) != NULL)  // returns true if ok
#define SEMAPHORE_WAIT(sem)             WaitForSingleObject(sem, INFINITE) // returns WAIT_OBJECT_0 if ok
#define SEMAPHORE_WAIT_WITH_TIMEOUT(sem, timeout_ms)    ((int) WaitForSingleObject(sem, timeout_ms)) // wait with timeout. use SEMAPHORE_WAIT_OK to check return value is ok.
#define SEMAPHORE_WAIT_OK(ret)          ((ret) == WAIT_OBJECT_0)      // returns true if ok. ret = integer
#define SEMAPHORE_WAIT_AND_OK(sem)      (WaitForSingleObject(sem, INFINITE) == WAIT_OBJECT_0) // returns true if ok
#define SEMAPHORE_RELEASE(sem)          ReleaseSemaphore(sem, 1, NULL) // returns nonzero if ok
#define SEMAPHORE_RELEASE_AND_OK(sem)   (ReleaseSemaphore(sem, 1, NULL) != 0) // returns true if ok
#define SEMAPHORE_DELETE(sem)           CloseHandle(sem) // returns nonzero on success
#define SEMAPHORE_DELETE_AND_OK(sem)    (ReleaseSemaphore(sem, 1, NULL) != 0) // returns true if ok

// print
#define SCANF scanf_s

// get time difference in milliseconds
inline uint32_t get_ticks_delta(uint32_t start) {
    uint32_t end = GetTickCount();
    return (end >= start ? end - start : end + (0xffffffff - start) + 1);
}

// color output in console, windows needs console handle
#define GET_STDOUT_HANDLE               HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE)
#define GET_STDERR_HANDLE               HANDLE hStdErr = GetStdHandle(STD_ERROR_HANDLE)
#define STDOUT_HANDLE                   hStdOut
#define STDERR_HANDLE                   hStdErr
#define COLOR_BLACK                     0
#define COLOR_BLUE                      1
#define COLOR_GREEN                     2
#define COLOR_AQUA                      3
#define COLOR_RED                       4
#define COLOR_PURPLE                    5
#define COLOR_YELLOW                    6
#define COLOR_WHITE                     7
#define COLOR_GRAY                      8
#define COLOR_LIGHT_BLUE                9
#define COLOR_LIGHT_GREEN               10
#define COLOR_LIGHT_AQUA                11
#define COLOR_LIGHT_RED                 12
#define COLOR_LIGHT_PURPLE              13
#define COLOR_LIGHT_YELLOW              14
#define COLOR_BRIGHT_WHITE              15
#define COLOR_DEFAULT                   COLOR_WHITE
#define SET_CONSOLE_COLOR_ERROR(hcon)   SetConsoleTextAttribute(hcon, COLOR_RED    )
#define SET_CONSOLE_COLOR_WARN(hcon)    SetConsoleTextAttribute(hcon, COLOR_YELLOW )
#define SET_CONSOLE_COLOR_NOTE(hcon)    SetConsoleTextAttribute(hcon, COLOR_GREEN  )
#define SET_CONSOLE_COLOR_INFO(hcon)    SetConsoleTextAttribute(hcon, COLOR_GREEN  )
#define RESET_CONSOLE_COLOR(hcon)       SetConsoleTextAttribute(hcon, COLOR_DEFAULT)

#else // LINUX / UNIX

#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#include <iostream>
#include <semaphore.h>                  // semaphore
#include <unistd.h>                     // getcwd

// threading
#include <pthread.h>                    // POSIX threads
#define THREAD_HANDLE                   pthread_t    // this is actually thread ID but used as a windows handle
#define INVALID_THREAD                  0
#define ZEROMEMORY(address, size)       memset(address, 0, size)
//#define THREAD_ENTRY                    void *
#define RETURN_FROM_THREAD(exit_code)   return((void*)((long)exit_code)) // without (long) get warning "-Wint-to-pointer-cast"

// synchronization barrier
#define BARRIER                         pthread_barrier_t
#define BARRIER_INIT(barrier,num)       pthread_barrier_init(&barrier, NULL, num) // returns 0 if ok
#define BARRIER_INIT_AND_OK(barr,num)   (pthread_barrier_init(&barr, NULL, num) == 0) // returns true if ok
#define BARRIER_WAIT(barrier)           pthread_barrier_wait(&barrier) // returns PTHREAD_BARRIER_SERIAL_THREAD for one thread otherwise 0
#define BARRIER_WAIT_AND_OK(barr,ret)   ( ((ret=pthread_barrier_wait(&barr)) == PTHREAD_BARRIER_SERIAL_THREAD) || ( ret == 0 ) ) // returns true if ok. ret = integer.
#define BARRIER_DELETE(barrier)         pthread_barrier_destroy(&barrier) // returns 0 if ok
#define BARRIER_DELETE_AND_OK(barrier)  (pthread_barrier_destroy(&barrier) == 0) // returns true if ok

// mutex
#define MUTEX                           pthread_mutex_t
#define MUTEX_INIT(mutex)               pthread_mutex_init(&mutex, NULL) // returns 0 if ok
#define MUTEX_INIT_AND_OK(mutex)        (pthread_mutex_init(&mutex, NULL) == 0) // returns true if ok
#define MUTEX_LOCK(mutex)               pthread_mutex_lock(&mutex) // returns 0 if ok
#define MUTEX_LOCK_AND_OK(mutex)        (pthread_mutex_lock(&mutex) == 0) // returns true if ok
#define MUTEX_UNLOCK(mutex)             pthread_mutex_unlock(&mutex) // returns 0 if ok
#define MUTEX_UNLOCK_AND_OK(mutex)      (pthread_mutex_unlock(&mutex) == 0) // returns true if ok
#define MUTEX_DELETE(mutex)             pthread_mutex_destroy(&mutex) // returns 0 if ok
#define MUTEX_DELETE_AND_OK(mutex)      (pthread_mutex_destroy(&mutex) == 0) // returns true if ok

// semaphore
#define SEMAPHORE                       sem_t
#define SEMAPHORE_INIT(sem, init)       sem_init(&sem, 0, init) // returns 0 if ok 
#define SEMAPHORE_INIT_AND_OK(sem,init) (sem_init(&sem, 0, init) == 0) // returns true if ok
#define SEMAPHORE_WAIT(sem)             sem_wait(&sem) // return 0 on success
#define SEMAPHORE_WAIT_OK(ret)          ((ret) == 0)      // returns true if ok. ret = integer.
#define SEMAPHORE_WAIT_AND_OK(sem)      (sem_wait(&sem) == 0) // returns true if ok
#define SEMAPHORE_RELEASE(sem)          sem_post(&sem) // return 0 on success
#define SEMAPHORE_RELEASE_AND_OK(sem)   (sem_post(&sem) == 0) // return true on success
#define SEMAPHORE_DELETE(sem)           sem_destroy(&sem) // return 0 on success
#define SEMAPHORE_DELETE_AND_OK(sem)    (sem_destroy(&sem) == 0) // return true on success
int SEMAPHORE_WAIT_WITH_TIMEOUT(SEMAPHORE & sem, long timeout); // wait with timeout. use SEMAPHORE_WAIT_OK to check return value is ok.

// print
#define SCANF scanf

// color output in console, windows needs console handle
// ANSI escape color codes foreground + backround
#define GET_STDOUT_HANDLE
#define GET_STDERR_HANDLE
#define STDOUT_HANDLE                   std::cout
#define STDERR_HANDLE                   std::cerr
#define COLOR_BLACK                     "\033[30m"
#define COLOR_RED                       "\033[31m"
#define COLOR_GREEN                     "\033[32m"
#define COLOR_YELLOW                    "\033[33m"
#define COLOR_BLUE                      "\033[34m"
#define COLOR_MAGENTA                   "\033[35m"
#define COLOR_CYAN                      "\033[36m"
#define COLOR_WHITE                     "\033[37m"
#define COLOR_BRIGHT_BLACK              "\033[90m"
#define COLOR_BRIGHT_RED                "\033[91m"
#define COLOR_BRIGHT_GREEN              "\033[92m"
#define COLOR_BRIGHT_YELLOW             "\033[93m"
#define COLOR_BRIGHT_BLUE               "\033[94m"
#define COLOR_BRIGHT_MAGENTA            "\033[95m"
#define COLOR_BRIGHT_CYAN               "\033[96m"
#define COLOR_BRIGHT_WHITE              "\033[97m"
#define COLOR_DEFAULT                   "\033[0m"
#define SET_CONSOLE_COLOR_ERROR(hcon)   hcon << COLOR_RED
#define SET_CONSOLE_COLOR_WARN(hcon)    hcon << COLOR_YELLOW
#define SET_CONSOLE_COLOR_NOTE(hcon)    hcon << COLOR_GREEN
#define SET_CONSOLE_COLOR_INFO(hcon)    hcon << COLOR_GREEN
#define RESET_CONSOLE_COLOR(hcon)       hcon << COLOR_DEFAULT

////////////////////////////////////////////////////////////////////////////////////////////////////
// windows equivalents
////////////////////////////////////////////////////////////////////////////////////////////////////

// sleeping time in ms
// uses nanosleep defined in <time.h>
void sleep_ms(unsigned long ms);

// measure elapsed time in milliseconds
// wraps over every 4295s = 71'
inline uint32_t GetTickCount(void) {
    static struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts); //CLOCK_MONOTONIC or CLOCK_REALTIME
    return ts.tv_sec*1000 + (ts.tv_nsec / 1000000);
}

// get time difference in milliseconds
inline uint32_t get_ticks_delta(uint32_t start) {
    static struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts); //CLOCK_MONOTONIC or CLOCK_REALTIME
    uint32_t end = ts.tv_sec*1000 + (ts.tv_nsec / 1000000);
    return (end >= start ? end - start : end + (0xffffffff - start) + 1);
}

#endif // end LINUX / UNIX

////////////////////////////////////////////////////////////////////////////////////////////////////
// 128bit arithmetics
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
// thread starts and shutdown
////////////////////////////////////////////////////////////////////////////////////////////////////

// thread timeout in ms
#define THREAD_TIMEOUT      1000
#define WAIT_INFINITE       -1

// thread entry point
#if defined(_WIN32) || defined (_WIN64)
#define thread_func     PTHREAD_START_ROUTINE
// helper thread starting point
extern DWORD WINAPI helper_thread_func(void* data);
#else
typedef void*(*thread_func)(void*);
// helper thread starting point
extern void* helper_thread_func(void* data);
#endif

int thread_start(thread_func func, THREAD_HANDLE *handle, void *data);
int thread_shutdown(THREAD_HANDLE handle, unsigned long timeout);
int thread_wait_shutdown(THREAD_HANDLE handle);

////////////////////////////////////////////////////////////////////////////////////////////////////
// queue for messages between threads
////////////////////////////////////////////////////////////////////////////////////////////////////

// thread commands
#define THREAD_SHUTDOWN                 0   // shutdown thread
#define THREAD_START                    1   // stread start
#define THREAD_CREATE_BASIC             20  // create atomic distribution with basic algorithm
#define THREAD_CREATE_METROPOLIS        21  // create atomic distribution with metropolis algorithm
//#define THREAD_FIND_NEAREST             30  // find nearest neighbors
#define THREAD_CREATE_MOL               40  // create molecule

class queue;

// entries in queue
class queue_entry {
private:
    friend class queue;                     // allow queue to access next
    class queue_entry *next;                // next entry or NULL
public:
    int cmd;                                // THREAD_ command
    int id;                                 // thread id, -1 = any
    int result;                             // -1 = not finished. 0 = ok, >0 = other result, <0 = error
    void *data;                             // user data

    queue_entry() { cmd = -1; id = -1; result = -1; data = NULL; next = NULL; };
    queue_entry(int _cmd, void *_data) { cmd = _cmd; id = -1; result = -1; data = _data; next = NULL; };
    queue_entry(int _cmd, void *_data, int _id) { cmd = _cmd; id = _id; result = -1; data = _data; next = NULL; };
    ~queue_entry() {};
};

// simple queue class
#ifdef _DEBUG
#define QUEUE_STATUS_OK                     true
#endif
class queue {
private:
    SEMAPHORE qsem;                         // wait & counting
    MUTEX qmutex;                           // protection
    queue_entry *first, *last;              // pointer to first and last queue entries, NULL if empty
#ifdef _DEBUG
    bool status;                            // for debugging
#endif
public:
    queue();                                // constructor
    ~queue();                               // destructor
#ifdef _DEBUG
    bool put(queue_entry *entry);               // append entry (can be several) to queue. returns QUEUE_STATUS_OK if ok.
    queue_entry *get(int max, long timeout_ms); // remove maximum max entries from queue (<0 for all, 0 = peek first entry without removing)
    bool get_status(void);                      // returns QUEUE_STATUS_OK if ok
#else
    void put(queue_entry *entry);               // append entry (can be several) to queue
    queue_entry *get(int max, long timeout_ms); // remove maximum max entries from queue (<0 for all, 0 = peek first entry without removing)
#endif
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// helper threads
////////////////////////////////////////////////////////////////////////////////////////////////////

// global mutex and barrier
extern MUTEX   global_mutex;
extern BARRIER global_barrier;

extern uint8_t num_species;    // TODO: can this be avoided?

#endif // _HELPER_THREADS
