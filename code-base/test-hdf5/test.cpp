#include <iostream>
#include "/cygdrive/c/hdf5/include/hdf5.h"

using namespace std;

int main (int argc, char** argv) {
    hid_t id = H5Tcreate(H5T_COMPOUND, 8); 

    cout << id;

    return 0;
}

