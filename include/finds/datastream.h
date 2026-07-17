#ifndef IO_HEADER
#define IO_HEADER

#include "system.h"
#include "error.h"

typedef enum {
    DATA_MODE_WRITING,
    DATA_MODE_READING
} ds_data_mode_e;

typedef struct {
    FILE *file;
} datastream_t;

error_e datastream_open_file(
    datastream_t *stream,
    const char *filename,
    const ds_data_mode_e mode);

error_e datastream_close(datastream_t *stream);

error_e datastream_write_system(
    const datastream_t stream,
    const fish_system_t *system,
    const double time);

#endif /* IO_HEADER */
