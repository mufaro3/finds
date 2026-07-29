/**
 * @file vector.h
 *
 * @brief Basic vector and matrix computations.
 *
 * This file contains routines for operating on 3-dimensional vectors and
 * 3-by-3 matrices. It's focused mostly on speed, and as such, is incredibly
 * small.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
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

/** @brief 3-dimensional cartesian vector implementation */
typedef union {
    struct { double x, y, z; };
    double data[3];
} double_3d_t;

/** @brief 3x3 matrix implementation */
typedef struct {
    double data[3][3];
} mat3_t;

/** @brief Shorthand for making D3 vectors */
#define D3(x,y,z) ((double_3d_t) {x,y,z})

/** @brief Fills a D3 with one value on each dimension */
#define D3_FILL(a) D3(a,a,a)

/** @brief A D3 of only 1s (1,1,1). */
#define D3_ONE D3_FILL(1)

/** @brief A D3 of only 0s (0,0,0)*/
#define D3_ZERO D3_FILL(0)

/** @brief Shorthand for creating 3x3 matrices */
#define MAT3(m11,m12,m13,m21,m22,m23,m31,m32,m33)   \
    ((mat3_t){ .data = {                            \
            {m11, m12, m13},                        \
            {m21, m22, m23},                        \
            {m31, m32, m33}                         \
        }})

/**
 * @brief Produces an x-rotation matrix of angle \f$\theta\f$.
 * @param[in] theta The x-rotation angle.
 * @return The x-rotation matrix.
 */
static inline mat3_t rotation_matrix_x(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        1, 0,  0,
        0, c, -s,
        0, s,  c);
}

/**
 * @brief Produces an y-rotation matrix of angle \f$\theta\f$.
 * @param[in] theta The y-rotation angle.
 * @return The y-rotation matrix.
 */
static inline mat3_t rotation_matrix_y(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        c, 0, s,
        0, 1, 0,
        -s, 0, c);
}

/**
 * @brief Produces an z-rotation matrix of angle \f$\theta\f$.
 * @param[in] theta The z-rotation angle.
 * @return The z-rotation matrix.
 */
static inline mat3_t rotation_matrix_z(const double theta)
{
    double s = sin(theta);
    double c = cos(theta);

    return MAT3(
        c, -s, 0,
        s, c, 0,
        0, 0, 1);
}

/**
 * @brief Multiplies two matrices a and b.
 * @param[in]  a  Matrix \f$A\f$.
 * @param[in]  b  Matrix \f$B\f$.
 * @return The multiplied matrix \f$AB\f$.
 */
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

/**
 * @brief Produces a 3-D rotation matrix for the roll, pitch, and yaw angles.
 * @param[in]  roll   The x-rotation angle, \f$\theta_x\f$.
 * @param[in]  pitch  The y-rotation angle, \f$\theta_y\f$.
 * @param[in]  yaw    The z-rotation angle, \f$\theta_z\f$.
 * @return The full rotation matrix,
 *         \f$R_z(\theta_z)R_y(\theta_y)R_x(\theta_x)\f$
 */
static inline mat3_t rotation_matrix(
    const double roll, const double pitch, const double yaw)
{
    mat3_t R_x = rotation_matrix_x(roll);
    mat3_t R_y = rotation_matrix_y(pitch);
    mat3_t R_z = rotation_matrix_z(yaw);

    return mat3_mult(R_z, mat3_mult(R_y, R_x));
}

/**
 * @brief Prints a D3 vector.
 * @param[in] v The vector.
 */
static inline void d3_print(const double_3d_t v)
{
    printf("[ %6.2e %6.2e %6.2e ]", v.x, v.y, v.z);
}

/**
 * @brief Adds together two D3 vectors.
 * @param[in] a Vector \f$\mathbf{a}\f$.
 * @param[in] b Vector \f$\mathbf{b}\f$.
 * @return \f$\mathbf{a} + \mathbf{b}\f$
 */
static inline double_3d_t d3_add(const double_3d_t a, const double_3d_t b)
{
    return (double_3d_t) {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

/**
 * @brief Scales a D3 vector by a scalar.
 * @param[in] vec    The vector to scale.
 * @param[in] scalar The scalar.
 * @return The vector scaled by the scalar.
 */
static inline double_3d_t d3_mult(const double_3d_t vec, const double scalar)
{
    return (double_3d_t) {
        vec.x * scalar,
        vec.y * scalar,
        vec.z * scalar
    };
}

/**
 * @brief Subtracts vector a by vector b.
 * @param[in] a Vector a.
 * @param[in] b Vector b.
 * @return a - b.
 */
static inline double_3d_t d3_sub(const double_3d_t a, const double_3d_t b)
{
    return d3_add( a, d3_mult( b, -1 ) );
}

/**
 * @brief Scales the vector a by 1 / scalar.
 * @param[in] a      The vector.
 * @param[in] scalar The scalar denominator.
 * @return a / scalar.
 */
static inline double_3d_t d3_div(const double_3d_t a, const double scalar)
{
    return d3_mult( a, 1 / scalar );
}

/**
 * @brief Computes the norm of a vector.
 * @param[in] v The vector.
 * @return The norm of the vector.
 */
static inline double d3_norm(const double_3d_t v)
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
}

/**
 * @brief Multiplies a 3x3 matrix by a D3 vector.
 * @param[in] mat The 3x3 matrix.
 * @param[in] v   The vector to multiply.
 * @return The matrix multiplied by the vector.
 */
static inline double_3d_t d3_mat3_mult(const mat3_t mat, const double_3d_t v)
{
    return D3(
        mat.data[0][0] * v.x + mat.data[0][1] * v.y + mat.data[0][2] * v.z,
        mat.data[1][0] * v.x + mat.data[1][1] * v.y + mat.data[1][2] * v.z,
        mat.data[2][0] * v.x + mat.data[2][1] * v.y + mat.data[2][2] * v.z
    );
}

/**
 * @brief Computes a vector containing the minimum values for each dimension
 *        of the two vectors.
 * @param[in] a Vector a.
 * @param[in] b Vector b.
 * @return The vector of minimums.
 */
static inline double_3d_t d3_min_each_dim(
    const double_3d_t a, const double_3d_t b)
{
    return D3(MIN(a.x,b.x), MIN(a.y,b.y), MIN(a.z,b.z));
}

/**
 * @brief Computes a vector containing the maximum values for each dimension
 *        of the two vectors.
 * @param[in] a Vector a.
 * @param[in] b Vector b.
 * @return The vector of maximums.
 */
static inline double_3d_t d3_max_each_dim(
    const double_3d_t a, const double_3d_t b)
{
    return D3(MAX(a.x,b.x), MAX(a.y,b.y), MAX(a.z,b.z));
}

/**
 * @brief Checks if a D3 vector contains any NaN values.
 * @param[in] v The vector.
 * @return Whether this vector contains any NaN values.
 */
static inline bool d3_has_nan(const double_3d_t v)
{
    return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

/**
 * @brief Returns the maximum of the components of a D3 vector.
 * @param[in] v The vector.
 * @return The maximum value of all of the components of a D3 vector.
 */
static inline double d3_max_component(const double_3d_t v)
{
    return MAX(v.z, MAX(v.y, v.x));
}

/**
 * @brief Rotates a D3 vector in 3-D.
 * @param[in]  v      The input vector \f$\mathbf{v}\f$
 * @param[in]  roll   The x-rotation angle, \f$\theta_x\f$.
 * @param[in]  pitch  The y-rotation angle, \f$\theta_y\f$.
 * @param[in]  yaw    The z-rotation angle, \f$\theta_z\f$.
 * @return The full rotation matrix,
 *         \f$R_z(\theta_z)R_y(\theta_y)R_x(\theta_x)\mathbf{v}\f$
 */
static inline double_3d_t d3_rotate(
    const double_3d_t v,
    const double roll,
    const double pitch,
    const double yaw)
{
    return d3_mat3_mult(rotation_matrix(roll, pitch, yaw), v);
}

/**
 * @brief Compares two D3 vectors to check if they are both close (essentially
 *        checking if they are equal within ABOUT_ZERO).
 * @param[in] a Vector a
 * @param[in] b Vector b
 * @return a == b within ABOUT_ZERO.
 */
static inline bool d3_is_close(const double_3d_t a, const double_3d_t b)
{
    return IS_CLOSE(a.x, b.x) && IS_CLOSE(a.y, b.y) && IS_CLOSE(a.z, b.z);
}

/**
 * @brief Computes the dot product between two vectors.
 * @param[in] a Vector a
 * @param[in] b Vector b
 * @return The dot product of a and b.
 */
static inline double d3_dot(const double_3d_t a, const double_3d_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* dealing with lists of vectors */

/**
 * @brief Allocates a list of D3 vectors.
 * @param[in] N The number of vectors to allocate.
 * @return A list of allocated vectors.
 */
static inline double_3d_t *d3_list_allocate(size_t N) {
    return calloc(N, sizeof(double_3d_t));
}

/**
 * @brief Adds two lists of D3 vectors together (indexwise).
 * @param[in] a The first list of D3 vectors.
 * @param[in] b The second list of D3 vectors.
 * @param[in] N The length of A and B to sum from.
 * @return a + b.
 */
static inline double_3d_t *d3_list_add(
    const double_3d_t *restrict a,
    const double_3d_t *restrict b,
    const size_t N)
{
    double_3d_t *sum = d3_list_allocate(N);
    for (size_t i = 0; i < N; ++i)
        sum[i] = d3_add(a[i], b[i]);
    return sum;
}

/**
 * @brief Prints a list of D3 vectors.
 * @param[in] lov The list of D3 vectors.
 * @param[in] N   The size of the list of vectors.
 */
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

/**
 * @brief Copies a D3 vector from src to dest.
 * @param[in] dest  Pointer to the destination vector.
 * @param[in] src   The source vector.
 */
static inline void d3_copy(double_3d_t *dest, const double_3d_t src)
{
    dest->x = src.x;
    dest->y = src.y;
    dest->z = src.z;
}

#endif
