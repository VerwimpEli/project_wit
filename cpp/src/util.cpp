#ifndef CPP_UTIL_H
#define CPP_UTIL_H

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

template<typename T>
void scale(T & v, int VW, int size){
    for (int i = 0; i < size; i++) {
        if (i < VW){ v[i] *= 0.5; }
        else if (i > size - VW - 1){ v[i] *= 0.5; }

        if (i % VW == 0){ v[i] *= 0.5; }
        else if ((i + 1) % VW == 0) {v[i] *= 0.5; }
    }
}

template<typename  T>
T min(T v1, T v2){
    return v1 < v2 ? v1 : v2;
}

template<typename  T>
T max(T v1, T v2){
    return v1 > v2 ? v1 : v2;
}

double max(std::vector<double> v){
    double result = v[0];
    for (size_t i = 1; i < v.size(); i++){
        result = v[i] > result ? v[i] : result;
    }
    return result;
}

double abs_v(double x){
    return x > 0 ?  x : -1*x;
}

double inf_norm_diff(std::vector<double> const &v, std::vector<double> const & v_old){
    double ch = 0.0;
    for (size_t i=0; i < v.size(); i++){
        ch = max(ch,abs_v(v[i] - v_old[i]));
    }
    return ch;
}

void Print(const std::vector<double>& v) {
    std::cout.precision(16);
    std::cout << "[" << v.size() << "] (\t";
    for (size_t i=0; i<v.size();i++){
        std::cout << v[i] << " ";
    }
    std::cout << ")" << std::endl;
}

void Print(double v){
    std::cout.precision(16);
    std::cout << v << std::endl;
}

void Print(double v, std::string name){
    std::cout.precision(16);
    std::cout << name << ":\t" << v << std::endl;
}

/**
 * Interpolates v as a 2D grid by averaging the values of the new grid points
 */
std::vector<double> interpolate(std::vector<double> const & v, double s) {
    std::vector<double> v_new(s*s*4, 0);
    for (int j=0; j<(s-1); j++){
        for (int i=0; i<(s-1); i++){
            v_new[2*i + 4*j*s] = v[i + j*s];
            v_new[2*i+1 + 4*j*s] = (v[i + j*s] + v[i+1 + j*s]) / 2;
            v_new[2*i + (4*j+2)*s] = (v[i + j*s] + v[i+(j+1)*s]) / 2;
            v_new[2*i+1 + (4*j+2)*s] = (v[i + j*s] + v[i+1 + j*s] + v[i + (j+1)*s] + v[i+1+(j+1)*s]) / 4;

        }
        v_new[2 * (s-1) + 2 * j * 2*s] = v[(s-1) + j * s];
        v_new[2 * (s-1) + 1 + 2 * j * 2*s] = v[(s-1) + j * s];
        v_new[2 * (s-1) + (2 * j + 1) * 2*s] = (v[(s-1) + j * s] + v[(s-1) + (j+1) * s]) / 2;
        v_new[2 * (s-1) + 1 + (2 * j + 1) * 2*s] = (v[(s-1) + j * s] + v[(s-1) + (j+1) * s]) / 2;
    }

    for (int i=0; i<s; i++) {
        v_new[2*i + 2*(s-1) * 2*s] = v[i + (s-1) * s];
        v_new[2*i+1 + 2 * (s-1) * 2*s] = (v[i + (s-1) * s] + v[i + 1 + (s-1) * s]) / 2;
        v_new[2*i + (2 * (s-1) + 1) * 2*s] = v[i + (s-1) * s];
        v_new[2*i+1 + (2 * (s-1) + 1) * 2*s] = (v[i + (s-1) * s] + v[i + (s-1) * s]) / 2;
    }
    return v_new;
}

#endif

