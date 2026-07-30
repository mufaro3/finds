/**
 * @file test_datastream.c
 *
 * @brief Unit tests for the Datastream I/O module.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include "tests.h"
#include <check.h>
#include <finds/datastream.h>
#include <finds/constants.h>

#define SYSTEM_SIZE 100
#define TEST_OUTPUT_FILE "test_datafile.h5"

#define __WRAP(code) ck_assert(code == ERR_OK)

/**
 * @brief Unit test for reading/writing datastream HDF5 files.
 *
 *  1. The deserialized version of the data should be identical
 *     to the data before serialization. (We should be able to
 *     fully retrieve the data afterward.)
 */
START_TEST(test_read_write_hdf5)
{
    /* write system */
    fish_system_t *input_system = \
        fish_system_generate_random(SYSTEM_SIZE, false);

    datastream_t stream;
    __WRAP(datastream_create_file(&stream, TEST_OUTPUT_FILE, SYSTEM_SIZE));
    __WRAP(datastream_write_system(&stream, input_system, 0));
    __WRAP(datastream_close(&stream));

    /* read system out */
    fish_system_t *output_system = NULL;
    double time;

    __WRAP(datastream_open_file(&stream, TEST_OUTPUT_FILE));
    __WRAP(datastream_read_frame(&stream, 0, &output_system, &time));
    __WRAP(datastream_close(&stream));

    /* check that the systems are exactly the same */
    ck_assert_uint_eq(input_system->size, output_system->size);
    for (size_t i = 0; i < input_system->size; ++i)
    {
        swimmer_t in  = input_system->swimmers[i];
        swimmer_t out = output_system->swimmers[i];

        ck_assert(d3_is_close(in.position, out.position));
        ck_assert(d3_is_close(in.orientation, out.orientation));
        ck_assert_double_eq(in.volumetric_flow_rate, out.volumetric_flow_rate);
        ck_assert_double_eq(in.length, out.length);
    }
}
END_TEST

Suite *test_suite_datastream(void)
{
    Suite *s = suite_create("datastream");

    TCase *tc = tcase_create("I/O");
    tcase_add_test(tc, test_read_write_hdf5);
    suite_add_tcase(s, tc);

    return s;
}
