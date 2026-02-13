
# Initial Orbit Determination: Extended Kalman Filter (GEOS-3)

![C++](https://img.shields.io/badge/C++-17-blue.svg)
![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

## 🛰️ Project Overview
This project implements a high-performance **Extended Kalman Filter (EKF)** for satellite orbit determination, specifically targeted at the **GEOS-3** satellite.

The primary objective was to migrate a legacy **MATLAB** codebase (2,800+ lines) to **Modern C++**, focusing on:
* **Memory Safety & Efficiency**: Replacing high-level abstractions with optimized memory management.
* **Performance Optimization**: Reducing execution latency through binary serialization and compiler optimizations.
* **Numerical Precision**: Ensuring sub-decimal accuracy compared to the original mathematical model.

## ⚡ Key Performance Metrics
By migrating from MATLAB to C++ and applying static optimization techniques (`-O3`), the system achieved an **~80% reduction in execution time**.

| Metric | Legacy MATLAB | Optimized C++ (`-O3`) | Improvement |
| :--- | :--- | :--- | :--- |
| **Execution Time** | 5.56 s | 1.06 s | **~5.2x Faster** |
| **I/O Strategy** | Text Parsing (`.txt`) | Binary Reading (`.bin`) | Reduced Latency |
| **Memory Usage** | Dynamic (High Overhead) | Optimized (Custom Matrix Class) | Efficient Allocation |

> **Note:** File loading time was reduced by approx. 3-4 seconds by converting large coefficient files (`DE430Coeff`) from text to binary format.

## 🛠️ Technical Features

### 1. Custom Linear Algebra Engine (`Matrix` Class)
Instead of relying on heavy external libraries, a custom `Matrix` class was engineered to handle vector/matrix operations.
* **Operator Overloading**: Implemented `+`, `-`, `*`, `/` for scalars and matrices to mimic mathematical syntax.
* **Memory Management**: Handled dynamic memory allocation/deallocation to prevent leaks.
* **Vector Operations**: Support for dot products, norms, and transpositions.

### 2. Physics & Math Implementation
The system implements complex astrodynamics algorithms, including:
* **Legendre Polynomials** for gravitational potential.
* **Runge-Kutta Integration (DEInteg)** for orbital propagation.
* **Chebyshev Polynomials** for ephemeris interpolation.
* **Coordinate Transformations**: ECI (Earth-Centered Inertial) to ECEF (Earth-Centered, Earth-Fixed).

### 3. Engineering Practices
* **Test-Driven Development (TDD)**: Unit tests were created for every module (e.g., `Mjday`, `AccelPointMass`) to validate C++ output against MATLAB baselines.
* **Profiling & Analysis**: Used **gprof** for performance profiling and **Understand** for static code analysis and complexity metrics.

## 📂 Project Structure
The project follows a standard C++ engineering structure:

```text
├── src/            # Source files (.cpp) - Main logic and algorithms
├── include/        # Header files (.h) - Declarations and global constants
├── data/           # Input data (Earth orientation parameters, coefficients)
├── tests/          # Unit tests for individual math functions
└── docs/           # Documentation generated via Doxygen

```

## 🚀 How to Build & Run

### Prerequisites

* C++ Compiler (GCC/G++ recommended)
* CMake (Optional, but recommended)

### Compilation (Manual)

You can compile the project using `g++` with optimization flags as demonstrated in the performance analysis:

```bash
# Compile with O3 optimization
g++ -O3 EKF_GEOS3.cpp ./src/*.cpp -o iod_ekf

# Run the executable
./iod_ekf

```

### Compilation (CMake)

```bash
mkdir build && cd build
cmake ..
make
./iod_ekf

```

## 📊 Documentation

Full documentation of functions and dependencies was generated using **Doxygen** and **GraphViz**. You can find the dependency graphs and call trees in the `docs/` folder.

## 👤 Author

**Martín Hernández Tonzán**

* [LinkedIn](https://www.google.com/search?q=https://www.linkedin.com/in/martin-hernandez-tonzan)
* [GitHub](https://www.google.com/search?q=https://github.com/mahernto)

---


Project originally developed for "Taller Transversal I" course at Universidad de La Rioja.
