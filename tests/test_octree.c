/**
 * @file test_octree.c
 *
 * @brief Unit tests for the linear octree generation.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <check.h>
#include <finds/barneshut.h>

#include "tests.h"

START_TEST(test_octree_generation)
{}
END_TEST

Suite *test_suite_octree(void)
{
    Suite *s = suite_create("octree");
    TCase *tc = tcase_create("Linear Octree Computation Testing");
    tcase_add_test(tc, test_octree_generation);
    suite_add_tcase(s, tc);

    return s;
}
