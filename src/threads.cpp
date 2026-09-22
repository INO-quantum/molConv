// helper threads
// created 24/10/2023 by Andi
// last change 24/10/2023 by Andi

#include <iostream>
#include <string.h>
#include "MolecularConversion.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// windows equivalents
////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(_WIN32) || defined (_WIN64) // Windows


#else // Linux / UNIX

// sleeping time in ms
// uses nanosleep defined in <time.h>
void sleep_ms(unsigned long ms) {
	struct timespec ts;
	ts.tv_sec = ms / 1000;
	ts.tv_nsec = (ms % 1000) * 1000000;
	nanosleep(&ts, NULL);
}

#endif // end linux

////////////////////////////////////////////////////////////////////////////////////////////////////
// thread starts and shutdown
////////////////////////////////////////////////////////////////////////////////////////////////////

// start thread
// return 0 and thread handle on success, otherwise error code
int thread_start(thread_func func, THREAD_HANDLE *handle, void *data) {
    int error = 0;
#if defined(_WIN32) || defined (_WIN64) // windows
    unsigned long id = 0;
    *handle = CreateThread(NULL,                        // security attributes
                0,                                      // stack size,0=default
                (LPTHREAD_START_ROUTINE)func,           // thread starting address
                (LPVOID)data,                           // parameters
                0,                                      // creation flags
                &id                                     // thread ID (unused)
                );
    if (*handle == INVALID_THREAD) {
        // error
        error = -1500;
    }
#else
    //std::cout << "create " << data << "..." << std::endl;
    error = pthread_create(handle,                      // thread id (used as handle)
                NULL,                                   // flags
                func,                                   // thread starting address
                data                                    // parameters
                );
    if (error != 0) { 
        // error
        *handle = INVALID_THREAD; 
        error = -1501;
    }
    /*else {
        std::cout << "create " << data << " ok" << std::endl;
    }*/
#endif
    return error;
}

int thread_shutdown(THREAD_HANDLE handle, unsigned long timeout) {
    int error = 0;
    // start shutdown of thread
    //err = shutdown(0);
    // for nonzero timeout wait until thread is terminated (regardless of error)
    if (timeout != 0) {
#if defined(_WIN32) || defined (_WIN64) // windows
        unsigned long result = WaitForSingleObject(handle, timeout);
        if (result == WAIT_TIMEOUT) {
            // thread still running after timeout ms: kill thread.
            error = -1600;
            TerminateThread(handle, error);
        }
        else {
            // thread terminated: get exit code
            unsigned long exitCode = 0;
            GetExitCodeThread(handle, &exitCode);
            error = (int)exitCode;
        }
#else
        void *exitCode;
        struct timespec ts;
        if (clock_gettime(CLOCK_REALTIME, &ts) == -1) error = -1601;
        else {
            ts.tv_sec += timeout/1000;
            error = pthread_timedjoin_np(handle, &exitCode, &ts);
            if(error == 0) {
                // all ok, return exit code of thread
                error = (long)exitCode;    // (long) avoids an error message
            }
            else if(error == ETIMEDOUT) {
                // timeout
                std::cout << "thread_shutdown: timeout!\n" << std::endl;
            }
            else {
                std::cout << "thread_shutdown: error " << error << std::endl;
                error = -1602;
            }
        }
#endif
    }
    return error;
}

// wait until server terminates
// returns exit code of thread or error code of pthread_join
// called from master thread
// TODO: not used here but Windows version is missing
int thread_wait_shutdown(THREAD_HANDLE handle) {
#if defined(_WIN32) || defined (_WIN64) // windows
    int err = 0;
#else
    int err = 0;
    void *exitCode;
    err = pthread_join(handle, &exitCode);
    if(err == 0) {
        // all ok, return exit code of thread
        err = (long)exitCode;    // (long) avoids an error message
    }
#endif
    return err;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// semaphore timed wait
////////////////////////////////////////////////////////////////////////////////////////////////////

// wait with timeout. use SEMAPHORE_WAIT_OK to check return value is ok.
// on timeout ret=-1 and errno == ETIMEDOUT
#if !defined(_WIN32) && !defined (_WIN64) // Linux/Unix. returns 0 if ok
int SEMAPHORE_WAIT_WITH_TIMEOUT(SEMAPHORE & sem, long timeout_ms) {
    int ret = -1;
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts) != -1) { 
        ts.tv_sec  += timeout_ms / 1000;
        ts.tv_nsec += (timeout_ms % 1000) * 1000000;
        // note: waiting might be interrupted by EINTR
        while ( ((ret = sem_timedwait(&sem, &ts)) == -1) && (errno == EINTR) );
    }
#ifdef _DEBUG    
    else {
        std::cout << "wait semaphore error get time!" << std::endl; 
    }
#endif
    return ret;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// queue for messages between threads
////////////////////////////////////////////////////////////////////////////////////////////////////

// default constructor
queue::queue() {
    first = last = NULL;
#ifdef _DEBUG
    status = SEMAPHORE_INIT_AND_OK(qsem, 0);
    status &= MUTEX_INIT_AND_OK(qmutex);
#else
    SEMAPHORE_INIT(qsem, 0);
    MUTEX_INIT(qmutex);
#endif
}

// destructor
queue::~queue() {
    if (first != NULL) {
        std::cout << std::endl << "queue: deleting non-empty queue! danger of memory leakage!" << std::endl << std::endl;
        while(first) {
            queue_entry *next = first->next; 
            delete first;
            first = next;
        }
    }
    last = NULL;
#ifdef _DEBUG
    status &= SEMAPHORE_DELETE_AND_OK(qsem);
    status &= MUTEX_DELETE_AND_OK(qmutex);
    if (status != QUEUE_STATUS_OK) {
        std::cout << std::endl << "queue: deleting error!" << std::endl << std::endl;
    }
#else
    SEMAPHORE_DELETE(qsem);
    MUTEX_DELETE(qmutex);
#endif
}

#ifdef _DEBUG
// returns QUEUE_STATUS_OK if ok
bool queue::get_status(void) { 
    if (status != QUEUE_STATUS_OK) {
        std::cout << std::endl << "queue: status error!" << std::endl;
        std::cout << std::endl << "queue: status    = " << status << std::endl;
        //std::cout << std::endl << "queue: semaphore = " << &qsem << std::endl;
        //std::cout << std::endl << "queue: mutex     = " << &mutex << std::endl;
    }
    return status; 
};
#endif

// append entry (can be several) to queue
// Attention: ensure entry->next == NULL for last entry.
#ifdef _DEBUG
bool queue::put(queue_entry *entry) {
    //int sval = 0;
    if (status == QUEUE_STATUS_OK) {
        status = MUTEX_LOCK_AND_OK(qmutex);
        if (status == QUEUE_STATUS_OK) {
            if (last == NULL) first = last = entry;
            else last = last->next = entry;
            while(last->next) {     // for several entries: increment semaphore for each of them and find last
                //std::cout << "put more than one!" << std::endl;
                status &= SEMAPHORE_RELEASE_AND_OK(qsem);
                last = last->next;
            }
            status &= SEMAPHORE_RELEASE_AND_OK(qsem);
            //sem_getvalue(&qsem, &sval);
            //std::cout << "new semaphore value = " << sval << std::endl;
            status &= MUTEX_UNLOCK_AND_OK(qmutex);
        }
    }
    return status;
}    
#else
void queue::put(queue_entry *entry) {
    MUTEX_LOCK(qmutex);
    if (last == NULL) first = last = entry;
    else last = last->next = entry;
    while(last->next) {     // for several entries: increment semaphore for each of them and find last
        SEMAPHORE_RELEASE(qsem);
        last = last->next;
    }
    MUTEX_UNLOCK(qmutex);
    SEMAPHORE_RELEASE(qsem);
}
#endif

// remove maximum max entries from queue (<0 for all, 0 = peek first entry without removing)
// if max = 0: does not wait but returns pointer to first entry without locking mutex or removing enrty. use only to check if queue is empty!
//             if this is NULL queue was at that moment empty, if its not NULL queue was not empty at that moment. 
//             Attention: do not use the returned pointer, since its not reliable if serveral threads can remove entries!
//             call get() with max!=0 to return the entry - this call will not wait unless another thread has already removed the entry.
// timeout_ms = maximum time the thread will wait in ms until the function returns with NULL. -1 = infinite.
// TODO: what happens with timeout_ms = 0?
#ifdef _DEBUG
class queue_entry *queue::get(int max, long timeout_ms) {
    class queue_entry *tmp = NULL;                 // returns NULL on timeout or error
    int s = -1;
    if (status == QUEUE_STATUS_OK) {
        if (max == 0) tmp = first;                      // peek first entry without locking or removing from queue! 
        else {
            if (timeout_ms == -1) {                     // wait without timeout
                s = SEMAPHORE_WAIT(qsem); 
            }
            else {                                      // wait with timeout
                s = SEMAPHORE_WAIT_WITH_TIMEOUT(qsem, timeout_ms);
            }
            if SEMAPHORE_WAIT_OK(s) {  // success
                status &= MUTEX_LOCK_AND_OK(qmutex);
                //int sval = -1;
                //sem_getvalue(&qsem, &sval);
                //std::cout << "wait with timeout ok first = " << first << std::endl;
                //std::cout << "new semaphore value = " << sval << std::endl;
                if (status == QUEUE_STATUS_OK) {
                    tmp = first;
                    if (max > 0) {                      // try to remove max entries
                        class queue_entry *l = tmp;
                        while(l->next && (--max > 0)) l = l->next;
                        first = l->next;
                        if (first == NULL) last = NULL; // all removed
                        else l->next = NULL;
                    }
                    else { // remove all entries
                        first = last = NULL;
                    }
                    status &= MUTEX_UNLOCK_AND_OK(qmutex);
                }
            }
        }
    }
    return tmp;
}
#else
class queue_entry *queue::get(int max, long timeout_ms) {
    class queue_entry *tmp = NULL;                 // returns NULL on timeout or error
    int s = -1;
    if (max == 0) tmp = first;                      // peek first entry without locking or removing from queue! 
    else {
        if (timeout_ms == -1) {                     // wait without timeout
            s = SEMAPHORE_WAIT(qsem); 
        }
        else {                                      // wait with timeout
            s = SEMAPHORE_WAIT_WITH_TIMEOUT(qsem, timeout_ms);
        }
        if SEMAPHORE_WAIT_OK(s) {  // success
            MUTEX_LOCK(qmutex);
            tmp = first;
            if (max > 0) {                      // try to remove max entries
                class queue_entry *l = tmp;
                while(l->next && (--max > 0)) l = l->next;
                first = l->next;
                if (first == NULL) last = NULL; // all removed
                else l->next = NULL;
            }
            else { // remove all entries
                first = last = NULL;
            }
            MUTEX_UNLOCK(qmutex);
        }
    }
    return tmp;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// helper thread
////////////////////////////////////////////////////////////////////////////////////////////////////

// helper thread entry
// returns 0 if ok, otherwise error code
#if defined(_WIN32) || defined (_WIN64)
DWORD WINAPI helper_thread_func(void *data) {
#else
void * helper_thread_func(void* data) {
#endif
    int error = 0;
    struct thread_data *pdata = reinterpret_cast<struct thread_data*>(data);
    uint8_t id          = pdata->id;  // thread id
    uint8_t num_threads = pdata->num_threads;
    class queue *recv   = pdata->to_helper;
    class queue *send   = pdata->from_helper;
    class queue_entry *entry;
    bool running = true;

    //MUTEX_LOCK(mutex);        
    //std::cout << "thread " << id << "/" << num_threads << " starting" << std::endl;
    //MUTEX_UNLOCK(mutex);
    
    // loop until shutdown
    while(running) {
        // recieve messages. returns NULL on timeout.
        entry = recv->get(1, THREAD_TIMEOUT);
        if (entry) {
            switch (entry->cmd) {
                case THREAD_SHUTDOWN:
                    //std::cout << "thread " << id << "/" << num_threads << " shutdown" << std::endl;
                    entry->result = 0;
                    running = false;
                    break;
                case THREAD_START:
                    //MUTEX_LOCK(mutex);        
                    //std::cout << "thread " << id << "/" << num_threads << " start" << std::endl;
                    //MUTEX_UNLOCK(mutex);        
                    entry->result = 0;
                    break;
                case THREAD_CREATE_BASIC:
                    //std::cout << "thread " << id << "/" << num_threads << " create atoms" << std::endl;
                    entry->result = CreateAtoms(pdata, reinterpret_cast<struct cmd_data*>(entry->data));
                    break;
                case THREAD_CREATE_METROPOLIS:
                    //std::cout << "thread " << id << "/" << num_threads << " create atoms" << std::endl;
                    entry->result = CreateAtoms_Metropolis(pdata, reinterpret_cast<struct cmd_data*>(entry->data));
                    break;
                //case THREAD_FIND_NEAREST:
                //    entry->result = FindNearest(pdata, reinterpret_cast<struct cmd_data*>(entry->data));
                //    break;
                case THREAD_CREATE_MOL:
                    //std::cout << "thread " << id << "/" << num_threads << " create molecule, seed " << pdata->seed << std::endl;
                    //entry->result = CreateMolecules_old(pdata);
#ifdef USE_OLD_MOLSEARCH
                    entry->result = CreateMolecules_v1_2(pdata, reinterpret_cast<struct cmd_data*>(entry->data));
#else                    
                    entry->result = CreateMolecules_linear(pdata, reinterpret_cast<struct cmd_data*>(entry->data));
#endif
                    //entry->result = CreateMolecules_deterministic(pdata);
                    break;
                default:
                    std::cout << "thread " << id << "/" << num_threads << " unknown command" << std::endl;
                    entry->result = -2;
                    running = false;
            }
            // each entry is replied back to primary thread after finished with thread id and result code.
            // this way we do not need to delete entry here.
            entry->id = id;
            send->put(entry);

            if (running && (num_threads > 1)) {
                // wait until all threads have finished last task
                // this ensures that each thread gets exactly 1x each shared task.
                int ret;
                if (!BARRIER_WAIT_AND_OK(::global_barrier, ret)) {
                    error = -1700;
                    break;
                }
            }
        }
    }

    if (error) { // print message only on error
        MUTEX_LOCK(::global_mutex);
        std::cout << "thread " << id << "/" << num_threads << " end with error " << error << std::endl;
        MUTEX_UNLOCK(::global_mutex);
    }
    //else {
    //    MUTEX_LOCK(::global_mutex);
    //    std::cout << "thread " << id << "/" << num_threads << " ending ok" << std::endl;
    //    MUTEX_UNLOCK(::global_mutex);
    //}
        
    // return error code, 0=ok
    RETURN_FROM_THREAD(error);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

