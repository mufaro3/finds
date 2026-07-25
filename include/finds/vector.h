/*
 * \file vector.h
 *
 * \brief Basic vector and matrix computations.
 *
 * This file contains routines for operating on 3-dimensional vectors and
 * 3-by-3 matrices. It's focused mostly on speed, and as such, is incredibly
 * small.
 *
 * \author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#ifndef VECTOR_HEADER
#define VECTOR_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"

typedef union {
    struct { double x, y, z; };
    double data[3];
} double_3d_t;

typedef struct {
    double data[3][3];
} mat3_t;

#define D3(x,y,z) ((double_3d_t) {x,y,z})
#define D3_FILL(a) D3(a,a,a)
#define D3_ONE D3_FILL(1)
#define D3_ZERO D3_FILL(0)

#define MAT3(m11,m12,m13,m21,m22,m23,m31,m32,m33)   \
    ((mat3_t){ .data = {                            \
            {m11, m12, m13},                        \
            {m21, m22, m23},                        \
            {m31, m32, m33}                         \
        }})

static inline mat3_t rotation_matrix_x(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        1, 0,  0,
        0, c, -s,
        0, s,  c);
}

static inline mat3_t rotation_matrix_y(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        c, 0, s,
        0, 1, 0,
        -s, 0, c);
}

static inline mat3_t rotation_matrix_z(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        c, -s, 0,
        s, c, 0,
        0, 0, 1);
}

static inline mat3_t mat3_mult(const mat3_t a, const mat3_t b)
{
    mat3_t result = {0};

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result.data[i][j] = 0.0;
            for (size_t k = 0; k < 3; ++k)
                result.data[i][j] += a.data[i][k] * b.data[k][j];
        }
    }

    return result;
}

static inline mat3_t rotation_matrix(
    const double roll, const double pitch, const double yaw)
{
    mat3_t R_x = rotation_matrix_x(roll);
    mat3_t R_y = rotation_matrix_y(pitch);
    mat3_t R_z = rotation_matrix_z(yaw);

    return mat3_mult(R_z, mat3_mult(R_y, R_x));
}

static inline void d3_print(const double_3d_t v)
{
    printf("[ %6.2e %6.2e %6.2e ]", v.x, v.y, v.z);
}

static inline double_3d_t d3_add(const double_3d_t a, const double_3d_t b)
{
    return (double_3d_t) {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

static inline double_3d_t d3_mult(const double_3d_t vec, const double scalar)
{
    return (double_3d_t) {
        vec.x * scalar,
        vec.y * scalar,
        vec.z * scalar
    };
}

static inline double_3d_t d3_sub(const double_3d_t a, const double_3d_t b)
{
    return d3_add( a, d3_mult( b, -1 ) );
}

static inline double_3d_t d3_div(const double_3d_t a, const double scalar)
{
    return d3_mult( a, 1 / scalar );
}

static inline double d3_norm(const double_3d_t v)
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
}

static inline double_3d_t d3_mat3_mult(const mat3_t mat, const double_3d_t v)
{
    return D3(
        mat.data[0][0] * v.x + mat.data[0][1] * v.y + mat.data[0][2] * v.z,
        mat.data[1][0] * v.x + mat.data[1][1] * v.y + mat.data[1][2] * v.z,
        mat.data[2][0] * v.x + mat.data[2][1] * v.y + mat.data[2][2] * v.z
    );
}

static inline double_3d_t d3_min_each_dim(
    const double_3d_t a, const double_3d_t b)
{
    return D3(MIN(a.x,b.x), MIN(a.y,b.y), MIN(a.z,b.z));
}

static inline double_3d_t d3_max_each_dim(
    const double_3d_t a, const double_3d_t b)
{
    return D3(MAX(a.x,b.x), MAX(a.y,b.y), MAX(a.z,b.z));
}

static inline double d3_max_component(const double_3d_t v)
{
    return MAX(v.z, MAX(v.y, v.x));
}

static inline double_3d_t d3_rotate(
    const double_3d_t v,
    const double roll,
    const double pitch,
    const double yaw)
{
    return d3_mat3_mult(rotation_matrix(roll, pitch, yaw), v);
}

static inline bool d3_is_close(const double_3d_t a, const double_3d_t b)
{
    return IS_CLOSE(a.x, b.x) && IS_CLOSE(a.y, b.y) && IS_CLOSE(a.z, b.z);
}

static inline double d3_dot(const double_3d_t a, const double_3d_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* dealing with lists of vectors */

static inline double_3d_t *d3_list_allocate(size_t N) {
    return calloc(N, sizeof(double_3d_t));
}

static inline double_3d_t *d3_list_add(
    const double_3d_t *a,
    const double_3d_t *b,
    size_t N)
{
    double_3d_t *sum = d3_list_allocate(N);
    for (size_t i = 0; i < N; ++i)
        sum[i] = d3_add(a[i], b[i]);
    return sum;
}

static inline void d3_list_print(const double_3d_t *lov, size_t N)
{
    printf("[ ");
    d3_print(lov[0]);

    for (size_t i = 1; i < N; ++i) {
        printf("\n  ");
        d3_print(lov[i]);
    }

    printf(" ]\n");
}

static inline void d3_copy(double_3d_t *dest, const double_3d_t src)
{
    dest->x = src.x;
    dest->y = src.y;
    dest->z = src.z;
}

#endif
