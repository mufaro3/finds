/**
 * @file linker_test.c
 *
 * @brief Just a linker test that includes all of the libraries used in FINDS.
 *
 * @author Mufaro Machaya
 *
 * License: MIT
 */
#include <argp.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <error.h>

#include <CL/cl.h>
#include <gsl/gsl_fit.h>
#include <gnuplot_i/gnuplot_i.h>
#include <progressbar/progressbar.h>
#include <tomlc17.h>
#include <hdf5.h>

#include <finds/constants.h>
#include <finds/datastream.h>
#include <finds/derivative.h>
#include <finds/evaluation.h>
#include <finds/integration.h>
#include <finds/simulation.h>
#include <finds/system.h>
#include <finds/util.h>
#include <finds/vector.h>

int main(void)
{
    puts("Linking successful!\n");
    return 0;
}
