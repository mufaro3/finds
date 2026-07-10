#ifndef VECTOR_HEADER
#define VECTOR_HEADER

#include <stdio.h>

typedef double double_3d_t[3];

#define DOUBLE_3D_ONE (double_3d_t) { 1, 1, 1 };

inline void double_3d_print(const double_3d_t vec)
{
    printf("%6.2e %6.2e %6.2e", vec[0], vec[1], vec[2]);
}

inline void double_3d_list_print(const double_3d_t *lov, size_t N)
{
    printf("[ ");
    double_3d_print(lov[0])

    for (size_t i = 1; i < N; ++i) {
        printf("\n  ")
        double_3d_print(lov[i]);
    }

    printf(" ]\n");
}

inline double_3d_t double_3d_add(const double_3d_t a, const double_3d_t b)
{
    return (double_3d_t) {
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2]
    }
}

inline double_3d_t double_3d_mult(const double_3d_t vec, double scalar)
{
    return (double_3d_t) {
        vec[0] * scalar,
        vec[1] * scalar,
        vec[2] * scalar
    }
}

inline double_3d_t double_3d_subtract(const double_3d_t a, const double_3d_t b)
{
    return double_3d_add( a, double_3d_mult( b, -1 ) );
}

inline double_3d_t double_3d_did(const double_3d_t a, double scalar)
{
    return double_3d_mult( a, 1 / scalar );
}

/* dealing with lists of vectors */

double_3d_t *double_3d_list_allocate(size_t N);
void *double_3D_list_free();

#endif
