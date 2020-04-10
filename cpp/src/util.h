#ifndef CPP_UTIL_H
#define CPP_UTIL_H

#include <vector>
#include <string>

template<typename T>
void scale(T & v, int VW, int size){
    for (int i = 0; i < size; i++) {
        if (i < VW){ v[i] *= 0.5; }
        else if (i > size - VW - 1){ v[i] *= 0.5; }

        if (i % VW == 0){ v[i] *= 0.5; }
        else if ((i + 1) % VW == 0) {v[i] *= 0.5; }
    }
};

template<typename  T>
T min(T v1, T v2) { return v1 < v2 ? v1 : v2; };

template<typename  T>
T max(T v1, T v2) { return v1 > v2 ? v1 : v2; };

double max(std::vector<double> v);
double abs_v(double x);
double inf_norm_diff(std::vector<double> const &v, std::vector<double> const & v_old);
void Print(const std::vector<double>& v);
void Print(double v);
void Print(double v, std::string name);


#endif //CPP_UTIL_H
