#ifndef CPP_UTIL_H
#define CPP_UTIL_H

#include <vector>
#include <string>

/**
 *  Scales the outer edges of the grid with 1/2. Grid is passed in v, as a vector (or iterable).
 *  VW is the grid width, and size the total grid size. Result is calculated in place.
 */
template<typename T>
void scale(T & v, int VW, int size){
    for (int i = 0; i < size; i++) {
        if (i < VW){ v[i] *= 0.5; }
        else if (i > size - VW - 1){ v[i] *= 0.5; }

        if (i % VW == 0){ v[i] *= 0.5; }
        else if ((i + 1) % VW == 0) {v[i] *= 0.5; }
    }
};

/**
 * returns the minimum of two values of type T
 */
template<typename  T>
T min(T v1, T v2) { return v1 < v2 ? v1 : v2; };

/**
 * Returns the maximum of two values of type T
 */
template<typename  T>
T max(T v1, T v2) { return v1 > v2 ? v1 : v2; };

/**
 * Returns the maximum element of a vector v
 */
double max(std::vector<double> v);
double abs_v(double x);
/**
 * Returns the infinity norm of the difference of v and v_old.
 */
double inf_norm_diff(std::vector<double> const &v, std::vector<double> const & v_old);
/**
 * Fancy print functions for different types.
 */
void Print(const std::vector<double>& v);
void Print(double v);
void Print(double v, std::string name);

/**
 * Interpolates a vector v, interpreted as a grid. New grid points are linearly interpolated.
 */
std::vector<double> interpolate(std::vector<double> const & v, double s);
#endif //CPP_UTIL_H
