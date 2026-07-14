#ifndef IO_HEADER
#define IO_HEADER

#include "system.h"

void *write_system_to_file(
    const FILE *output_file,
    const fish_system_t *system);

#endif /* IO_HEADER */
