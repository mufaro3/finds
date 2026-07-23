#ifndef ERROR_HEADER
#define ERROR_HEADER

#include <errno.h>

typedef enum {
    ERR_DATASTREAM_OPEN,
    ERR_DATASTREAM_CLOSED,
    ERR_DATASTREAM_WRITE,
    ERR_INVALID_ARG,
    ERR_ALLOC,
    ERR_STATE_CHANGE,
    ERR_INVALID_STATE,
    ERR_FAILURE,
    ERR_OK,

    ERR_TYPE_COUNT
} error_e;

/** Wrapping function for HDF5 errors */
#define WRAP_HDF5_CHECK(errcode, jmpto, result) \
    ({ if (result < 0) {                        \
            errcode = ERR_DATASTREAM_WRITE;     \
            goto jmpto;                         \
        } (result); })

/** Wrapping function for native errors */
#define WRAP_CHECK(errcode, jmpto, prev_errcode)    \
    if (prev_errcode != ERR_OK) {                   \
        errcode = prev_errcode;                     \
        goto jmpto;                                 \
    }

#define RAISE_ERROR(errtype, errmsg) ({                     \
            fprintf(stderr, "%s:%d %s: %s\n",               \
                __FILE__, __LINE__, __func__, (errmsg));    \
            (errtype);                                      \
                })

#define RAISE_ERROR_ERRNO(errtype, errmsg) ({           \
            fprintf(stderr, "%s:%d %s: %s (%s)\n",      \
                __FILE__, __LINE__,                     \
                __func__, (errmsg), strerror(errno));   \
            (errtype);                                  \
                })

#define GOTO_IF_NULL(errcode, jmpto, ptr)       \
    { if (ptr == NULL) {                        \
            errcode = ERR_FAILURE;              \
            goto jmpto;                         \
        }}

#endif /* ERROR_HEADER */
