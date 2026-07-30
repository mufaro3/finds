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

/**
 * @brief Unit test for reading/writing datastream HDF5 files.
 *
 * For an input shape of :math:`(N,6)`:
 *
 *  1. Each write should create only one row.
 *  2. Each row should have exactly :math:`6N+1` columns.
 *  3. The deserialized version of the data should be identical
 *     to the data before serialization. (We should be able to
 *     fully retrieve the data afterward.)
 *  4. Reading/writing should also work as-intended
 */
START_TEST(test_read_write_hdf5)
{}
END_TEST

Suite *test_suite_datastream(void)
{
    Suite *s = suite_create("datastream");

    TCase *tc = tcase_create("I/O");
    tcase_add_test(tc, test_read_write_hdf5);
    suite_add_tcase(s, tc);

    return s;
}
