#ifndef __CPUCTLIB_H
#define __CPUCTLIB_H

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>

#ifndef MAX
#define MAX(x, y) ((x)>=(y)?(x):(y))
#endif

#ifndef MIN
#define MIN(x, y) ((x)<=(y)?(x):(y))
#endif

#define IMG_Header_Size 64
#define LEN    2048
#define KERNEL_SIZE 3
#define Z_THRESHOLD 2.0
#define EPS    1e-6

//static const double PI = acos(-1.0);


#ifndef MAX
#define MAX(x, y) ((x)>=(y)?(x):(y))
#endif

#ifndef MIN
#define MIN(x, y) ((x)<=(y)?(x):(y))
#endif

#define IMG_Header_Size 64
#define LEN    2048
#define KERNEL_SIZE 3
#define Z_THRESHOLD 2.0
#define EPS    1e-6

static const double PI = acos(-1.0);

typedef long long           int64;
typedef unsigned short      T_half;
typedef float               T;





struct Matrix2dReal3x3 {
    Matrix2dReal3x3();
    virtual ~Matrix2dReal3x3();

    Matrix2dReal3x3(double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8);
    Matrix2dReal3x3 operator*(const Matrix2dReal3x3& m);
    Matrix2dReal3x3& operator=(const Matrix2dReal3x3& m);
#ifdef _WIN32 
    __declspec(align(16)) double data[9];
#else
    double data[9]  __attribute__((aligned(16)));
#endif
};

Matrix2dReal3x3 CalculateProjMat(int nx, int ny, double dx, double dy, int nu, double du, double angle_rad, double cor);

int findNextPowerOf2(int n);
void Filtering(T* input_buffer, int64 np, int64 nb, int64 nv, int64 nu, double* pReal_global, double* pImag_global, double* fft_kernel_, int64 fft_kernel_size);
void GetFilter(double* fft_kernel_, int64 fft_kernel_size, double frequency_max = 1.0);
void GetProjectionMatrix(Matrix2dReal3x3* mat_proj, int64 np, int64 nb, int direction, double scan_angle, double volume_postion_angle, int64 index, int64 nx, int64 ny, double dx, double dy, int64 nu, double du, double cor);
void BackProjection(T* input_buffer, T* output_buffer, int64 nb, int64 nv, int64 nu, int64 nx, int64 ny, int64 nz, Matrix2dReal3x3* mat_proj, int64 offset_CT, float RC);








#endif //__CPUCTLIB_H