/**
 * @file test_derivative.c
 *
 * @brief Unit tests for the derivative and interaction calculation module.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <stdlib.h>
#include <check.h>
#include <finds/derivative.h>

#include "tests.h"

START_TEST(test_compute_derivative_brute_force)
{}
END_TEST

START_TEST(test_compute_derivative_barnes_hut)
{}
END_TEST

START_TEST(test_compute_derivative_fast_multipole_method)
{}
END_TEST

START_TEST(test_regularization)
{}
END_TEST

Suite *test_suite_derivative(void)
{
    Suite *s = suite_create("derivative");

    TCase *tc = tcase_create("Derivative Computation Testing");
    tcase_add_test(tc, test_compute_derivative_brute_force);
    tcase_add_test(tc, test_compute_derivative_barnes_hut);
    tcase_add_test(tc, test_compute_derivative_fast_multipole_method);
    tcase_add_test(tc, test_regularization);
    suite_add_tcase(s, tc);

    return s;
}
