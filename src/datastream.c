#include <stdio.h>
#include <finds/error.h>
#include <finds/datastream.h>

error_e datastream_open_file(
    datastream_t *stream,
    const char *filename,
    const ds_data_mode_e mode)
{
    switch (mode) {
        case DATA_MODE_WRITING:
            stream->file = fopen(filename, "w");
            if (stream->file == NULL) {
                return RAISE_ERROR_ERRNO(ERR_DATASTREAM_OPEN,
                    "couldn't open file for writing");
            }
            break;
        case DATA_MODE_READING:
            stream->file = fopen(filename, "r");
            if (stream->file == NULL) {
                return RAISE_ERROR_ERRNO(ERR_DATASTREAM_OPEN,
                    "couldn't open file for reading");
            }
            break;
        default:
            return RAISE_ERROR(ERR_INVALID_ARG, "invalid data mode selection");
    }

    return ERR_OK;
}

error_e datastream_close(datastream_t *stream)
{
    fclose(stream->file);
    return ERR_OK;
}

error_e datastream_write_system(
    const datastream_t stream,
    const fish_system_t *system,
    const double time)
{
    fprintf(stream.file, "%lf", time);
    for (size_t i = 0; i < system->size; ++i) {
        /* position orientation */
        double_3d_t pos = system->swimmers[i].position;
        double_3d_t ori = system->swimmers[i].orientation;
        fprintf(stream.file, ",%lf,%lf,%lf,%lf,%lf,%lf",
            pos.x, pos.y, pos.z,
            ori.x, ori.y, ori.z);
    }
    fprintf(stream.file, "\n");
    return ERR_OK;
}
