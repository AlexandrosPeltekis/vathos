#pragma once
#include <iostream>
#include <type_traits>

// Base class: only holds the value
template<typename T>
struct DualBase {
    T val;

    DualBase() : val(0) {}
    DualBase(const T& v) : val(v) {}


    friend std::ostream& operator<<(std::ostream& os, const DualBase& x) {
        os << x.val;
        return os;
    }
};