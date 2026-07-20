#ifndef ERROR_HEADER
#define ERROR_HEADER

#include <errno.h>

typedef enum {
    ERR_DATASTREAM_OPEN,
    ERR_INVALID_ARG,
    ERR_ALLOC,
    ERR_STATE_CHANGE,
    ERR_INVALID_STATE,
    ERR_DATASTREAM_WRITE,
    ERR_FAILURE,
    ERR_OK,

    ERR_TYPE_COUNT
} error_e;

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

#endif /* ERROR_HEADER */
