#ifndef UTILITY_HEADER
#define UTILITY_HEADER

#include <gsl/gsl_vector.h>

void gsl_vector_print(const gsl_vector *v);
void gsl_vector_add_scalar(gsl_vector *v, double c);

#endif /* UTILITY_HEADER */
