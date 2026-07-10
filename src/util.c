#include <gsl/gsl_vector.h>
#include "util.h"

void gsl_vector_print(const gsl_vector *v)
{
    printf("[");

    for (size_t i = 0; i < v->size; ++i) {
        printf("%g", v->data[i * v->stride]);

        if (i + 1 != v->size)
            printf(", ");
    }

    printf("]\n");
}

void gsl_vector_add_scalar(gsl_vector *v, double c)
{
    for (size_t i = 0; i < v->size; ++i)
        v->data[i * v->stride] += c;
}
