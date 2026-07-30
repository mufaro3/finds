/**
 * @file main.c
 *
 * @brief Entry-point for running unit tests.
 *
 * @author Mufaro Machaya <mufaro2@student.ubc.ca>
 *
 * License: MIT
 */
#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "tests.h"

int main(void)
{
    SRunner *runner = srunner_create(test_suite_system());

    srunner_add_suite(runner, test_suite_datastream());
    srunner_add_suite(runner, test_suite_derivative());
    srunner_add_suite(runner, test_suite_integration());
    srunner_add_suite(runner, test_suite_octree());
    srunner_add_suite(runner, test_suite_simulation());

    srunner_run_all(runner, CK_NORMAL);

    int failures = srunner_ntests_failed(runner);

    srunner_free(runner);

    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
