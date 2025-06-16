# AutoDiff Library (Working Title)

## Overview

This project is a C++20 automatic differentiation (AD) library designed for scientific and engineering applications. It provides forward-mode AD using dual numbers, supporting first and second derivatives, as well as mixed partials in multiple dimensions. The library is modular, extensible, and tested with Catch2 v2.

---

## Features

- **Forward-mode automatic differentiation** via dual numbers.
- **First-order (`Dual1`) and second-order (`Dual2`) dual classes** for scalar and complex types.
- **Multi-dimensional support** with `Dual2_3d` for mixed partial derivatives.
- **Complex number support** for all dual types.
- **Operator overloading** for seamless mathematical expressions.
- **Test suite** using Catch2 v2.
- **CMake-based build system** (CMake ≥ 3.8, tested with 3.31.6-msvc6, Ninja generator).
- **C++20 standard** for modern language features.

---

## Goals

- Provide a robust, extensible, and easy-to-use AD library for scientific computing.
- Support higher-order and mixed derivatives for advanced applications (e.g., physics, optimization, machine learning).
- Enable both real and complex differentiation.
- Facilitate integration into larger C++ projects with modern build systems.

---

## Library Structure

- `dual1_class.h` — First-order dual numbers (forward-mode AD).
- `dual2_class.h` — Second-order dual numbers (forward-mode AD).
- `dual2_3d_class.h` — Second-order dual numbers with 3D mixed partials.
- `dualbase_class.h` — Common base functionality for dual classes.
- `reverse_ad.h` — (Planned) Reverse-mode AD support.
- `main.cpp` — Example usage and comprehensive operator tests.
- `tests/test_dual_classes.cpp` — Unit tests (Catch2 v2).

---

## Building

### Prerequisites

- CMake ≥ 3.8 (tested with 3.31.6-msvc6)
- C++20 compatible compiler (MSVC, GCC, Clang)
- [Ninja](https://ninja-build.org/) (recommended)
- [Catch2 v2](https://github.com/catchorg/Catch2) (for tests)

### Build Instructions
