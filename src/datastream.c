#include <errno.h>
#include <limits.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <hdf5.h>

#include <finds/vector.h>
#include <finds/error.h>
#include <finds/datastream.h>

int mkdir_p(const char *path, mode_t mode)
{
    char tmp[PATH_MAX];
    char *p;
    size_t len;

    if (path == NULL || *path == '\0') {
        errno = EINVAL;
        return -1;
    }

    len = strlen(path);
    if (len >= sizeof(tmp)) {
        errno = ENAMETOOLONG;
        return -1;
    }

    strcpy(tmp, path);

    /* Remove trailing '/' (except for "/") */
    if (len > 1 && tmp[len - 1] == '/')
        tmp[len - 1] = '\0';

    /* Skip leading '/' so absolute paths work */
    p = tmp + 1;

    for (; *p; ++p) {
        if (*p == '/') {
            *p = '\0';

            if (mkdir(tmp, mode) == -1 && errno != EEXIST)
                return -1;

            *p = '/';
        }
    }

    if (mkdir(tmp, mode) == -1 && errno != EEXIST)
        return -1;

    return 0;
}

static hid_t create_dataset_tensor(
    hid_t file, const char *path, hid_t space, hid_t dcpl)
{
    return H5Dcreate(
        file,
        path,
        H5T_NATIVE_DOUBLE,
        space,
        H5P_DEFAULT,
        dcpl,
        H5P_DEFAULT);
}

static error_e create_datasets(datastream_t *stream)
{
    error_e errcode = ERR_OK;

    const hsize_t dims1[1] = { 0 };
    const hsize_t dims2[2] = { 0, stream->n_particles    };
    const hsize_t dims3[3] = { 0, stream->n_particles, 3 };

    const hsize_t chunk1[1] = { 1 };
    const hsize_t chunk2[2] = { 1, stream->n_particles    };
    const hsize_t chunk3[3] = { 1, stream->n_particles, 3 };

    const hsize_t maxdims1[1] = { H5S_UNLIMITED };
    const hsize_t maxdims2[2] = { H5S_UNLIMITED, stream->n_particles };
    const hsize_t maxdims3[3] = { H5S_UNLIMITED, stream->n_particles, 3};

    hid_t space1 = H5Screate_simple(1, dims1, maxdims1);
    hid_t space2 = H5Screate_simple(2, dims2, maxdims2);
    hid_t space3 = H5Screate_simple(3, dims3, maxdims3);

    /* dimension lists for each tensor rank */
    hid_t dcpl1 = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dcpl2 = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dcpl3 = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_chunk(dcpl1, 1, chunk1);
    H5Pset_chunk(dcpl2, 2, chunk2);
    H5Pset_chunk(dcpl3, 3, chunk3);

    stream->time_dataset = create_dataset_tensor(
        stream->file, "/time", space1, dcpl1);

    stream->position_dataset = create_dataset_tensor(
        stream->file, "/position", space3, dcpl3);

    stream->orientation_dataset = create_dataset_tensor(
        stream->file, "/orientation", space3, dcpl3);

    stream->length_dataset = create_dataset_tensor(
        stream->file, "/length", space2, dcpl2);

    stream->sigma_dataset = create_dataset_tensor(
        stream->file, "/volumetric_flow_rate", space2, dcpl2);

    H5Sclose(space1);
    H5Sclose(space2);
    H5Sclose(space3);

    H5Pclose(dcpl1);
    H5Pclose(dcpl2);
    H5Pclose(dcpl3);

    return errcode;
}

error_e datastream_create_file(
    datastream_t *stream,
    const char *filename,
    const size_t N)
{
    error_e errcode = ERR_OK;

    stream->file = -1;
    stream->time_dataset = -1;
    stream->position_dataset = -1;
    stream->orientation_dataset = -1;
    stream->length_dataset = -1;
    stream->sigma_dataset = -1;

    stream->file = H5Fcreate(
        filename,
        H5F_ACC_TRUNC,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    stream->open = true;

    if (stream->file < 0)
        return RAISE_ERROR(
            ERR_DATASTREAM_OPEN,
            "couldn't create HDF5 file"
        );

    stream->open = true;
    stream->n_frames = 0;
    stream->n_particles = N;

    return create_datasets(stream);
}

error_e datastream_close(datastream_t *stream)
{
    if (!stream->open)
        return ERR_OK;

    if (stream->position_dataset >= 0)
        H5Dclose(stream->position_dataset);

    if (stream->orientation_dataset >= 0)
        H5Dclose(stream->orientation_dataset);

    if (stream->length_dataset >= 0)
        H5Dclose(stream->length_dataset);

    if (stream->sigma_dataset >= 0)
        H5Dclose(stream->sigma_dataset);

    if (stream->time_dataset >= 0)
        H5Dclose(stream->time_dataset);

    if (stream->file >= 0)
        H5Fclose(stream->file);

    stream->open = false;

    return ERR_OK;
}

#define WRAP_HDF5_CHECK(errcode, jmpto, result) \
    ({ if (result < 0) {                        \
            errcode = ERR_DATASTREAM_WRITE;     \
            goto jmpto;                         \
        } (result); })

#define WRAP_CHECK(errcode, jmpto, prev_errcode)    \
    ({ if (prev_errcode != ERR_OK) {                \
            errcode = prev_errcode;                 \
            goto jmpto;                             \
        } (prev_errcode); })

static error_e write_hdf5_frame(hid_t dataset, const void *data,
    const int rank, const hsize_t *start, const hsize_t *count)
{
    error_e errcode = ERR_OK;
    hid_t file_space = H5I_INVALID_HID;
    hid_t mem_space = H5I_INVALID_HID;

    file_space = H5Dget_space(dataset);

    // Assuming WRAP_HDF5_CHECK handles negative errors correctly
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
            start, NULL, count, NULL));

    /* FIX: Use the 'rank' variable instead of the broken sizeof() */
    mem_space = WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Screate_simple(rank, count, NULL));

    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_space,
            file_space, H5P_DEFAULT, data));

jmp_err:
    if (mem_space >= 0)
        H5Sclose(mem_space);
    if (file_space >= 0)
        H5Sclose(file_space);

    return errcode;
}

error_e datastream_write_system(
    datastream_t *stream,
    const fish_system_t *system,
    const double time)
{
    error_e errcode = ERR_OK;
    herr_t status;

    const size_t N = system->size;
    if (N != stream->n_particles) {
        errcode = RAISE_ERROR(ERR_STATE_CHANGE, "number of fish not constant");
        goto jmp_err;
    }

    /* setting up some temporary buffers for each data slot */
    double_3d_t *positions_buf    = calloc(N, sizeof(double_3d_t));
    double_3d_t *orientations_buf = calloc(N, sizeof(double_3d_t));

    double *lengths_buf = calloc(N, sizeof(double));
    double *sigmas_buf  = calloc(N, sizeof(double));

    /* if anything breaks, free everything */
    if (positions_buf == NULL || orientations_buf == NULL ||
        lengths_buf == NULL   || sigmas_buf == NULL)
    {
        errcode = RAISE_ERROR(ERR_ALLOC,
            "Could not allocate HDF5 buffers!");
        goto jmp_err;
    }

    /* Convert fish_system_t into contiguous arrays */
    for (size_t i = 0; i < N; ++i)
    {
        const swimmer_t *swimmer = &system->swimmers[i];

        d3_copy(&positions_buf[i], swimmer->position);
        d3_copy(&orientations_buf[i], swimmer->orientation);

        lengths_buf[i] = swimmer->length;
        sigmas_buf[i] = swimmer->volumetric_flow_rate;
    }

    /* Append the time */

    /* 1. Expand the dataset to accommodate the new frame */
    const hsize_t time_size[1] = { stream->n_frames + 1 };
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dset_extent(stream->time_dataset, time_size));

    /* 2. Get the updated file dataspace AFTER the extension */
    hid_t time_space = WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dget_space(stream->time_dataset));

    /* 3. Select the exact slot where the new data goes */
    const hsize_t time_start[1] = { stream->n_frames };
    const hsize_t time_count[1] = { 1 };
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Sselect_hyperslab(time_space, H5S_SELECT_SET,
            time_start, NULL, time_count, NULL));

    /* 4. Create the memory dataspace (1 element large) */
    hid_t time_memspace = WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Screate_simple(1, time_count, NULL));

    /* 5. Write data into the specific hyperslab selection */
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dwrite(stream->time_dataset, H5T_NATIVE_DOUBLE,
            time_memspace, time_space, H5P_DEFAULT, &time));

    /* 6. Clean up resources */
    H5Sclose(time_memspace);
    H5Sclose(time_space);

    /* Append the position and orientation data */
    const hsize_t dims3[3] = { stream->n_frames + 1, N, 3};
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dset_extent(stream->position_dataset, dims3));
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dset_extent(stream->orientation_dataset, dims3));

    const hsize_t start3[3] = { stream->n_frames, 0, 0 };
    const hsize_t count3[3] = { 1, N, 3 };

    WRAP_CHECK(errcode, jmp_err, write_hdf5_frame(stream->position_dataset,
            positions_buf, 3, start3, count3));

    WRAP_CHECK(errcode, jmp_err,
        write_hdf5_frame(stream->orientation_dataset,
            orientations_buf, 3, start3, count3));

    /* Append the length and volumetric flow rate data */
    const hsize_t dims2[2] = { stream->n_frames + 1, N };
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dset_extent(stream->length_dataset, dims2));
    WRAP_HDF5_CHECK(errcode, jmp_err,
        H5Dset_extent(stream->sigma_dataset, dims2));

    const hsize_t start2[2] = { stream->n_frames, 0 };
    const hsize_t count2[2] = { 1, N };
    WRAP_CHECK(errcode, jmp_err,
        write_hdf5_frame(stream->length_dataset,
            lengths_buf, 2, start2, count2));
    WRAP_CHECK(errcode, jmp_err,
        write_hdf5_frame(stream->sigma_dataset,
            sigmas_buf, 2, start2, count2));

jmp_skip:

    ++stream->n_frames;

    /* free everything */
jmp_err:
    if (positions_buf != NULL)
        free(positions_buf);
    if (orientations_buf != NULL)
        free(orientations_buf);
    if (lengths_buf != NULL)
        free(lengths_buf);
    if (sigmas_buf != NULL)
        free(sigmas_buf);

    return errcode;
}
