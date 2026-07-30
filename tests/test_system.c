/**
 * @file test_system.c
 *
 * @brief Unit tests for the fish system generation module.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <check.h>
#include <finds/system.h>

#include "tests.h"

START_TEST(test_system_random)
{
    ck_assert_double_eq(-1, 0);
}
END_TEST

/* todo: add additional tests */

Suite *test_suite_system(void)
{
    Suite *s = suite_create("system");
    TCase *tc = tcase_create("System Tests");
    tcase_add_test(tc, test_system_random);
    suite_add_tcase(s, tc);

    return s;
}
