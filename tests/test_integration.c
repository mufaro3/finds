/**
 * @file test_integration.c
 *
 * @brief Unit tests for the numeric integration module.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <check.h>
#include <finds/integration.h>

#include "tests.h"

START_TEST(test_step_euler)
{}
END_TEST

START_TEST(test_step_rk4)
{}
END_TEST

/* TODO: ADD TESTS FOR ALL INTEGRATORS */

Suite *test_suite_integration(void)
{
    Suite *s = suite_create("integration");
    TCase *tc = tcase_create("Integration Computation Testing");
    tcase_add_test(tc, test_step_euler);
    tcase_add_test(tc, test_step_rk4);
    suite_add_tcase(s, tc);

    return s;
}
