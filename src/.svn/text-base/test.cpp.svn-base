#include "SquRand.h"
#include "CubRand.h"
//#include "Cubic4d.h"
#include <iostream>
//#include <vector>
//#include <algorithm>

using namespace std;

void test1(size_t len, double fb) {
    CubRand L(len, 1, fb);

    size_t num_sites = L.getNumSites();
    bool error = true;
    int i, k;

    cout << num_sites << endl;

    for(i=0; i<num_sites; ++i) {
        cout << i << ": ";
        for(size_t j=0; j<6; ++j) {
            error = true;
            k = L.getNbr(i, j);
            cout << L.getNbr(i, j) << " ";

            for(int l=0; l<6; ++l)
            {
                if( L.getNbr(k, l) == i ) {
                    error = false;
                    break;
                }
            }
            if(error) {
                cout << "Malformed bond list." << endl;
                exit(0);
            }
        }
        cout << endl;
    }
    
}

void test2(size_t len, double fb) {

    size_t nsites = len*len*len;
    size_t nbonds = 3*nsites;

    CubRand L(len, 1, fb);
    size_t s, nidx;
    
    for(int i=0; i < nbonds; ++i) {
        nidx = i%3;
        s = (i-nidx)/3;
        cout << s << ' ' << L.getNbr(s, nidx) << endl;
    }
}

//void test3(size_t len, double fb) {
//    Cubic4d L(len);
//
//    size_t num_sites = L.getNumSites();
//    bool error = true;
//    int i, k;
//
//    cout << num_sites << endl;
//
//    for(i=0; i<num_sites; ++i) {
//        cout << i << ": ";
//        for(size_t j=0; j<8; ++j) {
//            error = true;
//            k = L.getNbr(i, j);
//            cout << L.getNbr(i, j) << " ";
//
//            for(int l=0; l<8; ++l)
//            {
//                if( L.getNbr(k, l) == i ) {
//                    error = false;
//                    break;
//                }
//            }
//            if(error) {
//                cout << "Malformed bond list." << endl;
//                exit(0);
//            }
//        }
//        cout << endl;
//    }
//    
//}

void test4(size_t len, double fb) {

    size_t nsites = len*len;
    size_t nbonds = 2*nsites;

    SquRand L(len, 1, fb);
    size_t s, nidx;
    
    for(int i=0; i < nbonds; ++i) {
        nidx = i%2;
        s = (i-nidx)/2;
        cout << s << ' ' << L.getNbr(s, nidx) << endl;
    }
}


int main(int argc, char **argv) {

    size_t len = 128;
    double fb = 0;

    if( argc >= 2 ) {
        len = atoi(argv[1]);
    }
    else if( argc >= 3 ) {
        fb = atof(argv[2]);
    }

    test4(len, fb);
}

