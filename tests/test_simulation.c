/**
 * @file test_simulation.c
 *
 * @brief Unit tests for the main simulation module.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <check.h>
#include <finds/simulation.h>

#include "tests.h"

START_TEST(test_simulation_non_adaptive_rk4)
{}
END_TEST

START_TEST(test_simulation_adaptive_rk45)
{}
END_TEST

Suite *test_suite_simulation(void)
{
    Suite *s = suite_create("simulation");
    TCase *tc = tcase_create("Adaptive/Non-Adaptive Simulation Testing");
    tcase_add_test(tc, test_simulation_non_adaptive_rk4);
    tcase_add_test(tc, test_simulation_adaptive_rk45);
    suite_add_tcase(s, tc);

    return s;
}
