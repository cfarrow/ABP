/* Point information used in lattice calculations. 

These classes are responsible for translating lattice indexes to and from points
in 2, 3, and 4 dimensions.

*/

#ifndef POINT_H
#define POINT_H


/* Utility for putting a value between 0 (inclusive) and len (exclusive) */
static inline int bounded(int v, int len)
{
    v = v < 0 ? v + len : v;
    v = v >= len ? v % len : v;
    return v;
}


struct Point2d {

    Point2d() {
        update(0, 1);
    }

    Point2d(size_t index, size_t len) {
        update(index, len);
    }

    void update(size_t index, size_t len) {
        x = index % len;
        y = (index - x) / len;
        length = len;
    }

    /* Shift x and/or y and return the new index. Does not modify internal data.
    */
    size_t shift(int sx, int sy) {
        size_t x_, y_;
        x_ = bounded(x + sx, length);
        y_ = bounded(y + sy, length);
        return x_ + y_ * length;
    }

    size_t x, y, length;
};


struct Point3d {

    Point3d() {
        update(0, 1);
    }

    Point3d(size_t index, size_t len) {
        update(index, len);
    }

    void update(size_t index, size_t len) {
        length = len;
        length2 = len * len;
        x = index % len;
        y = ((index - x) / len) % len;
        z = (index - x - y * len) / length2;
    }

    /* Shift x, y, z and return the new index. Does not modify internal data.
    */
    size_t shift(int sx, int sy, int sz) {
        size_t x_, y_, z_;
        x_ = bounded(x + sx, length);
        y_ = bounded(y + sy, length);
        z_ = bounded(z + sz, length);
        return x_ + y_ * length + z_ * length2;
    }

    size_t x, y, z, length, length2;
};


struct Point4d {

    Point4d() {
        update(0, 1);
    }

    Point4d(size_t index, size_t len) {
        update(index, len);
    }

    void update(size_t index, size_t len) {
        length = len;
        length2 = len * len;
        length3 = len * length2;
        w = index % length;
        x = ((index - w) % length2) / length;
        y = ((index - w - x * length) % length3) / length2;
        z = (index - w - x * length - y * length2) / length3;
    }

    /* Shift w, x, y, z and return the new index. Does not modify internal data.
    */
    size_t shift(int sw, int sx, int sy, int sz) {
        size_t w_, x_, y_, z_;
        w_ = bounded(w + sw, length);
        x_ = bounded(x + sx, length);
        y_ = bounded(y + sy, length);
        z_ = bounded(z + sz, length);
        return w_ + x_ * length + y_ * length2 + z_ * length3;
    }

    size_t w, x, y, z, length, length2, length3;
};

#endif