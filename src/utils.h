#ifndef UTILS_H
#define UTILS_H

inline void swap( int &a, int &b ) {
    int temp = a;
    a = b;
    b = temp;
}

inline void swap( size_t &a, size_t &b ) {
    size_t temp = a;
    a = b;
    b = temp;
}

inline void swap( double &a, double &b ) {
    double temp = a;
    a = b;
    b = temp;
}

#endif
