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

#define FILENAME "dset.h5"

int
main(void)
{
    const char *filename = "test.h5";

    /*
     * Create file
     */
    hid_t file = H5Fcreate(
        filename,
        H5F_ACC_TRUNC,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
        return 1;
    }


    /*
     * Data to write
     *
     * 4 particles, each with x,y,z
     */
    double data[4][3] = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0},
        {10.0, 11.0, 12.0}
    };


    /*
     * Create dataspace
     */
    hsize_t dims[2] = {4, 3};

    hid_t dataspace = H5Screate_simple(
        2,
        dims,
        NULL
    );

    if (dataspace < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
        H5Fclose(file);
        return 1;
    }


    /*
     * Create dataset
     */
    hid_t dataset = H5Dcreate(
        file,
        "/position",
        H5T_NATIVE_DOUBLE,
        dataspace,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    if (dataset < 0) {
        fprintf(stderr, "Failed to create dataset\n");
        H5Sclose(dataspace);
        H5Fclose(file);
        return 1;
    }


    /*
     * Write data
     */
    herr_t status = H5Dwrite(
        dataset,
        H5T_NATIVE_DOUBLE,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        data
    );

    if (status < 0) {
        fprintf(stderr, "Failed to write data\n");
    }


    /*
     * Cleanup
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);

    return 0;
}
