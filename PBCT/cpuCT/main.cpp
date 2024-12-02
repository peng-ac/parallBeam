#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include "../cpuCTlib/cpuCTlib.h"



void generate_gaussian_kernel_2D(T kernel[KERNEL_SIZE][KERNEL_SIZE], T sigma = 1.0) {
    T sum = 0.0;
    int half_size = KERNEL_SIZE / 2;

    for (int x = -half_size; x <= half_size; x++) {
        for (int y = -half_size; y <= half_size; y++) {
            T value = exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
            kernel[x + half_size][y + half_size] = value;
            sum += value;
        }
    }

    for (int x = 0; x < KERNEL_SIZE; x++) {
        for (int y = 0; y < KERNEL_SIZE; y++) {
            kernel[x][y] /= sum;
        }
    }
}

struct IMG_Header {
    char            head[2];
    short            comment_length;
    short            width;
    short            height;
    short            x_offset;
    short            y_offset;
    short            type;
    char            reserved[63 - 13];
    char* comment;
};
typedef struct IMG_Header IMGHeader;

void loadImage(char* filePath, T* image, T_half* loader, int64 nv, int64 nu, const T* dark_flat = NULL, const double scale = 0.0) {
    IMGHeader        img_head;
    char* comment;
    char             buffer[IMG_Header_Size];

    FILE* fp = fopen(filePath, "rb");

    fread(buffer, sizeof(char), IMG_Header_Size, fp);
    memcpy(&img_head.head, buffer, 2);    // char x 2
    memcpy(&img_head.comment_length, buffer + 2, 2);    // short
    memcpy(&img_head.width, buffer + 4, 2);    // short
    memcpy(&img_head.height, buffer + 6, 2);    // short
    memcpy(&img_head.x_offset, buffer + 8, 2);    // short
    memcpy(&img_head.y_offset, buffer + 10, 2);    // short
    memcpy(&img_head.type, buffer + 12, 2);    // short
    memcpy(&img_head.reserved, buffer + 14, 50);    // char x 50
    comment = (char*)malloc(img_head.comment_length + 1);
    fread(comment, sizeof(char), img_head.comment_length, fp);
    comment[img_head.comment_length] = '\0';

    free(comment);

    fread(loader, sizeof(T_half), nu * nv, fp);

    fclose(fp);

    T gauss_kernel[KERNEL_SIZE][KERNEL_SIZE];
    generate_gaussian_kernel_2D(gauss_kernel, 1.0);
    int half_size = KERNEL_SIZE / 2;

#pragma omp parallel for collapse(2)
    for (int64 i = 0; i < nv; i++) {
        for (int64 j = 0; j < nu; j++) {
            double sum = 0.0;
            for (int kx = -half_size; kx <= half_size; kx++) {
                for (int ky = -half_size; ky <= half_size; ky++) {
                    int ix = i + kx;
                    int iy = j + ky;

                    if (ix >= 0 && ix < nv && iy >= 0 && iy < nu) {
                        sum += (T)loader[ix * nu + iy] * gauss_kernel[kx + half_size][ky + half_size];
                    }
                }
            }

            if (dark_flat == NULL) {
                image[i * nu + j] = sum;
            }
            else {
                const double denom = MAX(EPS, (1.0 - scale) * dark_flat[1 * nv * nu + i * nu + j] + scale * dark_flat[2 * nv * nu + i * nu + j] - dark_flat[i * nu + j]);
                const double img = MAX(EPS, sum - dark_flat[i * nu + j]);
                image[i * nu + j] = MAX(0, -log(img / denom));
            }
        }
    }
}

void LoadData(T* input_buffer, int64 np, int64 nb, int64 nv, int64 nu, T* loader, T* dark_flat, int64 index, char* format) {
    for (int64 cnt = 0; cnt < nb; cnt++) {
        char inFilePath[LEN];
        sprintf(inFilePath, format, index + cnt + 2);
        loadImage(inFilePath, input_buffer + cnt * nv * nu, (T_half*)loader, nv, nu, dark_flat, 1.0 * (index + cnt) / np);
    }
}

int main(int argc, char* argv[]) {

    std::string input_dir = argv[1];
    std::string output_dir = argv[2];
    const int64 nv = atoi(argv[3]);
    const int64 nu = atoi(argv[4]);
    const int64 np = atoi(argv[5]);
    const int64 offset_CT = atoi(argv[6]);
    const float RC = atoi(argv[7]); // position of rotation axis
    const int64 resolution = atoi(argv[8]);
    const int64 digit = atoi(argv[9]);
    const int64 nb_max = atoi(argv[10]);

    const double du = 12.06, dv = 12.06;
    const int64  nx = resolution, ny = resolution, nz = nv;
    const double scale = (offset_CT == 1) ? 2.0 * RC / nu : 1.0;
    const double dx = (1.0 * du * nu * scale) / nx, dy = (1.0 * du * nu * scale) / ny, dz = dv;
    const int    direction = -1;
    const double scan_angle = (offset_CT == 1) ? 360 : 180;
    const double volume_postion_angle = 0.0;
    const double cor = nu / 2 - RC;

    const int fft_kernel_size = findNextPowerOf2(MAX(resolution, nu));
    const double frequency_max = 0.95; // cutoff frequency

    double* fft_kernel_ = static_cast<double*>(malloc(fft_kernel_size * sizeof(double)));
    GetFilter(fft_kernel_, fft_kernel_size, frequency_max);

    T* loader = static_cast<T*>(malloc((int64)nv * nu * sizeof(T)));
    T* dark_flat = static_cast<T*>(malloc((int64)3 * nv * nu * sizeof(T)));
    T* input_buffer = static_cast<T*>(malloc((int64)nb_max * nv * nu * sizeof(T)));
    T* output_buffer = static_cast<T*>(malloc((int64)nx * ny * nz * sizeof(T)));
    memset(output_buffer, 0, nx * ny * nz * sizeof(T));

    double* pReal_global = static_cast<double*>(malloc(nb_max * nv * fft_kernel_size * sizeof(double)));
    double* pImag_global = static_cast<double*>(malloc(nb_max * nv * fft_kernel_size * sizeof(double)));
    Matrix2dReal3x3* mat_proj = static_cast<Matrix2dReal3x3*>(malloc(nb_max * sizeof(Matrix2dReal3x3)));

    char format[LEN];
    sprintf(format, "%s/q%%0%lldlld.img", input_dir.c_str(), digit);
    {
        char inFilePath[LEN];
        sprintf(inFilePath, "%s/dark.img", input_dir.c_str());
        loadImage(inFilePath, dark_flat, (T_half*)loader, nv, nu);

        sprintf(inFilePath, format, 1);
        loadImage(inFilePath, dark_flat + nv * nu, (T_half*)loader, nv, nu);

        sprintf(inFilePath, format, np + 2);
        loadImage(inFilePath, dark_flat + 2 * nv * nu, (T_half*)loader, nv, nu);
    }

    for (int64 index = 0; index < np; index += nb_max) {
        const int64 nb = MIN(nb_max, np - index);

        LoadData(input_buffer, np, nb, nv, nu, loader, dark_flat, index, format);
        Filtering(input_buffer, np, nb, nv, nu, pReal_global, pImag_global, fft_kernel_, fft_kernel_size);
        GetProjectionMatrix(mat_proj, np, nb, direction, scan_angle, volume_postion_angle, index, nx, ny, dx, dy, nu, du, cor);
        BackProjection(input_buffer, output_buffer, nb, nv, nu, nx, ny, nz, mat_proj, offset_CT, RC);

        printf("[%lld/%lld]\n", index + nb, np);
    }

#pragma omp parallel for
    for (int64 k = 0; k < nz; k++) {
        char outFilePath[LEN];
        sprintf(outFilePath, "%s/volume_%05lld.raw", output_dir.c_str(), k);
        FILE* fp = fopen(outFilePath, "wb");
        if (fp) {
            fwrite(output_buffer + k * nx * ny, sizeof(*output_buffer), nx * ny, fp);
            fclose(fp);
        }
        else {
            printf("Store Error\n");
        }
    }

    free(mat_proj);
    free(fft_kernel_);
    free(dark_flat);
    free(loader);
    free(input_buffer);
    free(output_buffer);
    free(pReal_global);
    free(pImag_global);

    return 0;
}