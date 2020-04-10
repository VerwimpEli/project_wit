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
    for (int i = 1; i < v.size(); i++){
        result = v[i] > result ? v[i] : result;
    }
}

double abs_v(double x){
    return x > 0 ?  x : -1*x;
}

double inf_norm_diff(std::vector<double> const &v, std::vector<double> const & v_old){
    double ch = 0.0;
    for (int i=0; i < v.size(); i++){
        ch = max(ch,abs_v(v[i] - v_old[i]));
    }
    return ch;
}

void Print(const std::vector<double>& v) {
    std::cout.precision(16);
    std::cout << "[" << v.size() << "] (\t";
    for (int i=0; i<v.size();i++){
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

#endif

