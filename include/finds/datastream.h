#ifndef IO_HEADER
#define IO_HEADER

#include <hdf5.h>
#include <stdbool.h>
#include "system.h"
#include "error.h"

typedef struct {
    bool open;
    hid_t file;

    hid_t time_dataset;
    hid_t position_dataset;
    hid_t orientation_dataset;
    hid_t length_dataset;
    hid_t sigma_dataset;

    size_t n_particles;
    size_t n_frames;
} datastream_t;


error_e mkdir_p(const char *path, const mode_t mode);

error_e datastream_open_file(
    datastream_t *stream,
    const char *filename);

error_e datastream_read_frame(
    datastream_t *stream,
    const size_t frame,
    fish_system_t **dest_ptr,
    double *time);

error_e datastream_create_file(
    datastream_t *stream,
    const char *filename,
    const size_t N);

error_e datastream_close(datastream_t *stream);

error_e datastream_write_system(
    datastream_t *stream,
    const fish_system_t *system,
    const double time);

#endif /* IO_HEADER */
