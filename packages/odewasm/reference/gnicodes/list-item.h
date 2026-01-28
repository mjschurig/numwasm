#include <string>
#include <iostream>

// Class to construct a list ... even if the list-items
// are spread across multiple source files, and even if
// the list needs to be constructed before main() begins.
// See menu-list.c and menu-list-orange.c for examples
// of usage.

// Principles of operation:  For list_item instances
// that are declared outside of any function, the
// rootptr should be declared as
//      list_item* my_root_ptr(0);
// where we rely on the fact that scalar types (including
// pointers) are initialized _before_ things that require
// constructors are initialized.  That is, we need the
// rootptr to be initialized very early.
//
// To say the same thing another way, we cannot use an
// instance of std::list, because the list has to be
// initialized before it can be used, and we have no
// control over the order in which constructors will be
// called, for things declared outside of any function.
//
// In this scenario, an instance of class list_guard
// should also be declared, outside of any function,
// preferably right after the rootptr is declared.
// This guarantees that rootptr will always be a valid
// pointer, and rootptr->next will safely evaluate to
// zero (when the list is empty) or to a pointer to the
// first element of the list (when that exists).
//
// Also: the list_guard destructor will clean up the dummy
// root-item if necessary.  This upholds the rule that a
// "new" in a constructor should be matched by a "delete"
// in the corresponding destructor ... and indeed compensates
// for the fact that this rule is /not/ upheld by the
// list_item constructor.
//
// Initialization of the rootptr cannot be built into
// the constructor for the ordinary list item.  The whole
// point of the list_guard class is to cope with the case
// where there /aren't/ any list_item constructors, i.e.
// the empty list.
//
// After initialization, rootptr does not point to the
// first "real" item of the list.  Instead it points to
// a dummy root-item, which is an instance of list_item,
// and that in turn has a pointer to the first real
// (non-dummy) item in the list, if any.  The rationale
// is simple:  We want every non-dummy item to have a
// parent (of type list_item) that can be cleanly updated
// when the item is destructed.  Updating the rootptr
// itself would not be nearly so clean.
//
// In contrast, for instances of list_item that are declared
// within some function, the order of initialization is
// under control, so the dummy root-item can be declared
// directly and then used as the root of the list, with
// no need for a root-pointer.

class list_item{
public:
  list_item* prev;      // conventional FLP-linked list
  list_item* next;      // ditto
  std::string name;

  void attach(list_item& root){
    prev = &root;       // we always add ourself at the head
    next = prev->next;

    prev->next = this;          // tell our parent about us
    if(next) next->prev = this; // tell our kid about us
  }

// Normal constructor #1, for attaching via a root-pointer
// that may need initializing:
  list_item(list_item*& rootptr, const std::string& _name) :
    name(_name)
  {
    if (!rootptr) {
      rootptr = new list_item(_name);
    }                   // rootptr is now guaranteed valid
    attach(*rootptr);
  }

// Normal constructor #2, for attaching to an existing dummy root-item:
  list_item(list_item& root, const std::string& _name) :
    name(_name)
  {
    attach(root);
  }

// Similar to the above:
// Rarely-used constructor, for attaching via a root-pointer
// that has already been initialized, perhaps by taking the
// address of an existing dummy root-item:
  list_item(list_item* const& rootptr, const std::string& _name) :
    name(_name)
  {
    attach(*rootptr);
  }

// Special constructor for creating a dummy root-item.
// Must not be used outside of a function;  use
//   list_item* rootptr(0) and list_guard ....
// for lists declared outside of any function.
  list_item(const std::string _name) :
    prev(0),
    next(0),
    name("dummy: root of list: " + _name) {}

  void dump (const int skip = 0) const {
    const list_item* item(this);
    for (int ii = 0; ii < skip; ii++) {
      if (item) item = item->next;
    }
    while (item) {
      std::cout << item->name << std::endl;
      item = (list_item*)item->next;
    }
    std::cout << "===========" << std::endl;
  }

// Destructor:
  ~list_item(){
#ifdef DEBUG_LIST_ITEM
    using namespace std;
    cout << "   Destructor: " << name << endl;
#endif
    if (next) next->prev = prev;
    if (prev) prev->next = next;
  }
};


// Special class used mainly for the side effects of
// its constructor and destructor, namely to initialize
// rootptr (if necessary) to make sure it points to a
// valid dummy root-item ... and conversely to delete
// the dummy root-item when the time comes.
//
// Normally every time you declare a rootptr you declare
// a list_guard to watch over it.
//
class list_guard{
  list_item** ptrptr;
public:
  list_guard(list_item*& rootptr, const std::string _name)
    : ptrptr(&rootptr)
  {
    if (!rootptr) {
      rootptr = new list_item(_name);
    }
  }

  ~list_guard() {
    if (*ptrptr) {
      delete *ptrptr;
      *ptrptr = 0;
    }
  }
};

// Use this within a function to create a root-pointer
// with a built-in guard.
// That is, we automagically delete the dummy root-item
// when we go out of scope.
//
// You must *not* use this for any lists that get constructed
// before main() begins.  Use plain list_item* and list_guard
// instead.  See the Principles of Operation above.
//
// I don't know whether this is more elegant or less elegant
// than just creating your own dummy root-item and using it
// directly.  (Outside of any function you can't do that,
// either.)  The advantage is that the list_root::root can
// be used in exactly the same ways as a guarded list_item*.
// The disadvantage is that it creates a pointer variable
// that wouldn't otherwise be needed.

class list_root {
public:
  list_item* root;

  list_root(const string& _name = "whatever") :
    root(new list_item(_name))
  {}

  ~list_root() {
    if (root) delete root;
  }
};