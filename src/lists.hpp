// template class double-linked list
// elements of class T must define: 
//  friend class double_linked_list<T, hash_type>;
//  private T *prev,*next; 
//  private hash_type hash_key;

#ifndef _DOUBLE_LINKED_LIST
#define _DOUBLE_LINKED_LIST

#define ERROR_LIST_ADD          -102    // add(first, last, num) with improper last or num.

#define LOOKUP_SIZE     100     // initial maximum size of lookup elements

template <class T, typename hash_type=uint64_t>  
class double_linked_list {
private:
    T* first;       // first element or NULL when empty
    T* last;        // last element of NULL when empty
    T** lookup;     // list of few in-between elements for faster sorting. NULL until add_sorted called.
    uint64_t num;   // number of elements in list. 0 initially
public:
    double_linked_list() { 
        first = last = NULL; 
        lookup = NULL;
        num = 0; 
    };
    ~double_linked_list() { 
        if (lookup) delete [] lookup;
        lookup = NULL;
    };

    // add elements at end of list
    void add(T* elements) {
        if (elements) { 
            ++num;
            T* last_element = elements;
            while(last_element->next) {
                last_element = last_element->next;
                ++num;
            }
            elements->prev = last;
            if (last == NULL) {
                first = elements;
                last  = last_element;
            }
            else {
                last->next = elements;
                last = last_element;
            }
        }
    };
    
    // add num elements from first until last at end of list.
    // this is more efficient than add(T* elements) above.
    // this does not need to count and find last element.
    // Attention: ensure that first and last are properly double-linked and contains num elements!
    //            otherwise, list will be corrupted and app will sooner or later crash!
    // in _DEBUG mode throws ERROR_LIST_ADD
    void add(T* first, T* last, uint64_t num) {
        if (first && last && (num>0)) { 
#ifdef _DEBUG
            uint64_t count = 1;
            T* test = first;
            while(test->next) {
                ++count;
                test = test->next;
            }
            if ((test != last) || (count != num)) {
                throw ERROR_LIST_ADD;
            }
#endif        
            first->prev = this->last;
            last ->next = NULL;
            if (this->last == NULL) {
                this->first      = first;
                this->last       = last;
            }
            else {
                this->last->next = first;
                this->last       = last;
            }
            this->num += num;
        }
    };

    // add elements at beginning of list
    void add_first(T* elements) {
        if (elements) {
            ++num;
            T* last_element = elements;
            while(last_element->next) {
                last_element = last_element->next;
                ++num;
            }
            elements->prev = NULL;
            if (first == NULL) {
                first = elements;
                last  = last_element;
            }
            else {
                last_element->next = first;
                first->prev = last_element;
                first = elements;
            }
        }
    }
    
    // insert elements after the given element
    void insert(T* elements, T *after) {
        if (elements && after) {
            T* last_element = elements;
            while(last_element->next) {
                last_element = last_element->next;
                ++num;
            }
            if (last_element) {
                ++num;
                elements->prev     = after;
                last_element->next = after->next;
                after->next        = elements;
                after->next->prev  = last_element;
            }
        }
    }
    
    // insert elements sorted by increasing hash_key.
    // hash_key does not need to be unique.
    // note: at the moment just linearly searches elements in list assuming it is already sorted.
    // TODO: improve the performance with lookup table. started but not used and not finished.
    //       have to find a way to insert elements into lookup table capturing the most dense regions in list.
    //       the index into lookup table could be searched fast using min and max of all hash_key.
    //       once the nearest element in lookup table is found the list can be searched starting from this.
    void add_sorted(T* elements) {
        /*if (lookup == NULL ) {
            // first call: populate lookup table and sort existing elements
            uint32_t n = num > LOOKUP_SIZE ? LOOKUP_SIZE : num;
            lookup = new T*[n];
            for (uint32_t i = 0; i < n; ++i) lookup[n] = NULL;
        }*/
        if (!elements) return;
        // go through elements
        T *e_new = elements, *e_list = NULL;
        T *old = elements;

        while(e_new) {
            elements = elements->next; // save next element
            hash_type hash_key = e_new->hash_key;

            // find starting list element in lookup table
            // lookup table contains in-between elements or NULL at end.
            if (lookup) {
                uint64_t i = 0;
                e_list = lookup[i];
                while(e_list) {
                    if (hash_key < e_list->hash_key) {
                        // insert e_new between e_list->prev and e_list
                        e_list = e_list->prev == NULL ? first : e_list->prev;
                        break;
                    }
                    else if (hash_key == e_list->hash_key) {
                        // insert e_new after e_list before first element with larger hash_key
                        break;
                    }
                    if (lookup[i+1] != NULL) e_list = lookup[++i];
                    else break; // insert e_new after last element of lookup table
                }
            }
            else e_list = first; // without lookup table start at first element

            // go through list
            while(e_list) {
                if (hash_key < e_list->hash_key) {
                    // insert e_new before e_list. 
                    // if e_new has same hash_key it is inserted after last e_list with same hash_key.
                    e_new->prev = e_list->prev;
                    e_new->next = e_list;
                    if (e_list->prev == NULL) first = e_list->prev = e_new;
                    else {
                        e_list->prev = e_list->prev->next = e_new;
                    }
                    break;
                }
                e_list = e_list->next;
            }
            if (e_list == NULL) {
                // hash_key is >= last lement hash_key: insert as new last element
                e_new->prev = last;
                e_new->next = NULL;
                if (last == NULL) first = last = e_new; // empty list
                else {
                    last = last->next = e_new;
                }
            }
            // new element inserted
            ++num;
            // get next element in elements list 
            e_new = elements;
        }

        e_list = first;
        uint64_t i = 0;
        while(e_list && (i < 10)) {
            e_new = old;
            while(e_new) {
                if (e_new == e_list) break;
                e_new = e_new->next;
            }
            //std::cout << i << " : " << e_list->hash_key << ((e_new == e_list) ? " *" : "") << std::endl; 
            e_list = e_list->next;
            ++i;
        }  

    }
    
    // remove this element from list
    void remove(T* element) {
        if (element) {
            if (element->prev == NULL) first               = element->next;
            else                       element->prev->next = element->next;
            if (element->next == NULL) last                = element->prev;
            else                       element->next->prev = element->prev;
            element->prev = element->next = NULL;
            --num;
        }
    }

    // remove all elements including element from list
    // returns first element of removed elements
    T* remove_until(T* element) {
        if (!element) return NULL;
        // count number of elements to be removed
        T* r_list = first;
        --num;
        while((r_list) && (r_list != element)) {
            --num;
            r_list = r_list->next;
        }
        r_list = first;
        // remove elements
        first = element->next;
        if (element->next) element->next->prev = NULL;
        else               last = NULL;
        element->next = NULL;
        // return removed elements
        return r_list;
    }

    // remove all elementsfrom list
    // returns first element
    T* remove_all(void) {
        T* first = this->first;
        this->first = this->last = NULL;
        this->num = 0;
        return first;
    }

    // return number of elements in list
#ifdef _DEBUG
    uint64_t get_num(void) { 
        if ((num == 0) && ( (first != NULL) || (last != NULL) ) ) {
            std::cout << std::endl << "list error: 0 entries but first = " << first << ", last = " << last << std::endl << std::endl;
        }
        return num; 
    }
#else
    uint64_t get_num(void) { return num; }
#endif
    
    // get first element in list
    T* get_first(void) { return first; }
    
    // get last element in list
    T* get_last(void)  { return last; }

    // get element with given offset index from first element. 
    // index = 0 gives first. if element does not exists returns NULL.
    T* get_first(uint64_t index) { 
        T* element = first;
        while ( (index != 0) && (element != NULL) ) {
            --index;
            element = element->next;
        }
        return element;
    }

};

#endif // _DOUBLE_LINKED_LIST
