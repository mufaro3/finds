#include <check.h>
#include <math.h>

double square(double x)
{
    return x * x;
}

double distance(double x1, double y1,
                double x2, double y2)
{
    double dx = x2 - x1;
    double dy = y2 - y1;

    return sqrt(dx * dx + dy * dy);
}

START_TEST(test_square_positive)
{
    ck_assert_double_eq(square(5.0), 25.0);
}
END_TEST


START_TEST(test_square_negative)
{
    ck_assert_double_eq(square(-3.0), 9.0);
}
END_TEST


START_TEST(test_distance)
{
    double d = distance(0.0, 0.0, 3.0, 4.0);

    ck_assert_double_eq_tol(d, 5.0, 1e-10);
}
END_TEST


Suite *math_suite(void)
{
    Suite *s = suite_create("Math");

    TCase *tc = tcase_create("Core");

    tcase_add_test(tc, test_square_positive);
    tcase_add_test(tc, test_square_negative);
    tcase_add_test(tc, test_distance);

    suite_add_tcase(s, tc);

    return s;
}


int main(void)
{
    int failures;

    Suite *s = math_suite();
    SRunner *runner = srunner_create(s);

    srunner_run_all(runner, CK_NORMAL);

    failures = srunner_ntests_failed(runner);

    srunner_free(runner);

    return failures == 0 ? 0 : 1;
}
