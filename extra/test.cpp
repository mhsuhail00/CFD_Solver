#include <iostream>
#include <chrono>
#include <blitz/array.h>
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>    // for elementwise operations like +
#include <xtensor/xio.hpp> 
using namespace std;
using namespace blitz;

const int NX = 1000;

void testBlitz() {
    Array<double,3> a{NX, NX, NX}, b{NX, NX, NX}, c{NX, NX, NX};

    // Initialize arrays
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j)
            for (int k = 0; k < NX; ++k){
                a(i,j,k) = random();
                b(i,j,k) = random();
            }

    auto start = chrono::high_resolution_clock::now();
    c = a + b;
    auto end = chrono::high_resolution_clock::now();

    double seconds = chrono::duration<double>(end - start).count();
    cout << "Vector addition (Blitz) took " << seconds << " seconds" << endl;
}

double*** allocate3D(int depth, int rows, int cols) {
    double*** arr = new double**[depth];
    for (int i = 0; i < depth; i++) {
        arr[i] = new double*[rows];
        for (int j = 0; j < rows; j++) {
            arr[i][j] = new double[cols];
        }
    }
    return arr;
}

void deallocate3D(double*** arr, int depth, int rows) {
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < rows; j++) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}

void testRaw() {
    // Allocate contiguous memory
    double ***a = allocate3D(NX, NX, NX);
    double ***b = allocate3D(NX, NX, NX);
    double ***c = allocate3D(NX, NX, NX);

    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j)
            for (int k = 0; k < NX; ++k){
                a[i][j][k] = random();
                b[i][j][k] = random();
            }

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j)
            for (int k = 0; k < NX; ++k)
                c[i][j][k] = a[i][j][k] + b[i][j][k];
    auto end = chrono::high_resolution_clock::now();
    double seconds = chrono::duration<double>(end - start).count();
    cout << "Vector addition (Raw) took " << seconds << " seconds" << endl;

    deallocate3D(a, NX, NX);
    deallocate3D(b, NX, NX);
    deallocate3D(c, NX, NX);
}

void testXtensor() {
    cout << "\n--- xtensor (row-major, SIMD) ---" << endl;

    // Force row-major (C order)
    xt::xarray<double, xt::layout_type::row_major> a = xt::random::rand<double>({NX, NX, NX});
    xt::xarray<double, xt::layout_type::row_major> b = xt::random::rand<double>({NX, NX, NX});
    xt::xarray<double, xt::layout_type::row_major> c({NX, NX, NX});

    auto start = chrono::high_resolution_clock::now();
    c = a + b;   // SIMD-fused
    auto end = chrono::high_resolution_clock::now();

    double seconds = chrono::duration<double>(end - start).count();
    cout << "Vector addition (xtensor row-major) took " << seconds << " seconds\n";
}
int main() {
    // testBlitz();
    testRaw();
    testXtensor();
    return 0;
}
