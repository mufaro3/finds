/**
 * @file integration.c
 *
 * @brief Integrator implementations.
 *
 * This file contains the numerical integration schemes to be used in
 * performing a simulation over some time domain, including forward Euler and
 * embedded Runge-Kutta of variable order.
 *
 * @author Mufaro J. Machaya
 *
 * License: MIT
 */
#include <finds/derivative.h>
#include <finds/integration.h>
#include <finds/error.h>

/** Named lookup table for integration methods */
const named_enum_t INTEGRATION_METHODS_TABLE[INTEGRATION_METHODS_COUNT] = {
    { EULER,          "euler" },
    { RUNGE_KUTTA_2,  "RK2"   },
    { RUNGE_KUTTA_4,  "RK4"   },
    { RUNGE_KUTTA_5,  "RK5"   },
    { RUNGE_KUTTA_6,  "RK6"   },

    /* adaptive */
    { RUNGE_KUTTA_32, "RK32"  },
    { RUNGE_KUTTA_54, "RK54"  }
};

/**
 * @brief Computes \f$y + y' \delta t\f$.
 *
 * @param[in] state      \f$y\f$
 * @param[in] derivative \f$y'\f$
 * @param[in] time_step  \f$\delta t\f$
 */
static fish_system_t *euler_advance(
    const fish_system_t *state,
    const system_derivative_t *derivative,
    const double time_step)
{
    fish_system_t *advanced_state = fish_system_copy(state);
    for (size_t i = 0; i < state->size; ++i)
    {
        swimmer_t *swimmer = &advanced_state->swimmers[i];
        double_3d_t trans_deriv = derivative->translational[i];
        double_3d_t rot_deriv = derivative->rotational[i];

        swimmer->position = d3_add(swimmer->position,
            d3_mult(trans_deriv, time_step));
        swimmer->orientation = d3_add(swimmer->orientation,
            d3_mult(rot_deriv, time_step));
    }
    fish_system_normalize_orientation(advanced_state);
    return advanced_state;
}

/**
 * @brief Naive Euler-Stepping without Error Estimation
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state \f$y(t + \delta t)\f$.
 */
fish_system_t *step_euler(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    system_derivative_t *derivative = \
        compute_system_derivative(state, dc_opts);
    *error = 0;
    fish_system_t *advanced = euler_advance(state, derivative, time_step);
    derivative_destroy(&derivative);
    return advanced;
}

/**
 * @brief Naive 2nd-order Runge-Kutta
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state \f$y(t + \delta t)\f$.
 */
fish_system_t *step_rk2(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    *error = 0;
    error_e errcode = ERR_OK;

    system_derivative_t *k1 = NULL;
    system_derivative_t *k2 = NULL;
    fish_system_t *mid = NULL;
    fish_system_t *next = NULL;

    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    mid = euler_advance(state, k1, time_step / 2.0);
    GOTO_IF_NULL(errcode, jmp_err, mid);

    k2 = compute_system_derivative(mid, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    next = euler_advance(state, k2, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);

    fish_system_destroy(&mid);

    if (errcode != ERR_OK) {
        fish_system_destroy(&next);
        return NULL;
    }

    return next;
}

/**
 * @brief Naive 4th-order Runge-Kutta
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state \f$y(t + \delta t)\f$.
 */
fish_system_t *step_rk4(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    *error = 0;
    error_e errcode = ERR_OK;

    system_derivative_t *k1, *k2, *k3, *k4, *k_avg;
    fish_system_t *y1, *y2, *y3, *next;

    k1    = NULL;
    k2    = NULL;
    k3    = NULL;
    k4    = NULL;
    k_avg = NULL;
    y1    = NULL;
    y2    = NULL;
    y3    = NULL;
    next  = NULL;

    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    y1 = euler_advance(state, k1, time_step / 2.0);
    GOTO_IF_NULL(errcode, jmp_err, y1);

    k2 = compute_system_derivative(y1, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    y2 = euler_advance(state, k2, time_step / 2.0);
    GOTO_IF_NULL(errcode, jmp_err, y2);

    k3 = compute_system_derivative(y2, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k3);

    y3 = euler_advance(state, k3, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y3);

    k4 = compute_system_derivative(y3, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k4);

    const derivative_weight_t k_avg_terms[] = {
        {1.0/6.0,  k1}, {1.0/3.0, k2}, {1.0/3.0, k3}, {1.0/6.0, k4}
    };
    k_avg = derivative_average(k_avg_terms, state->size, 4);
    GOTO_IF_NULL(errcode, jmp_err, k_avg);

    next = euler_advance(state, k_avg, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);
    derivative_destroy(&k3);
    derivative_destroy(&k4);
    derivative_destroy(&k_avg);

    fish_system_destroy(&y1);
    fish_system_destroy(&y2);
    fish_system_destroy(&y3);

    if (errcode != ERR_OK) {
        fish_system_destroy(&next);
        return NULL;
    }

    return next;
}

/**
 * @brief Naive 5th-order Runge-Kutta
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state \f$y(t + \delta t)\f$.
 */
fish_system_t *step_rk5(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    *error = 0;
    error_e errcode = ERR_OK;
    system_derivative_t *k1, *k2, *k3, *k4, *k5, *k6;
    system_derivative_t *s2, *s3, *s4, *s5, *s6, *k_avg;
    fish_system_t *y2, *y3, *y4, *y5, *y6, *next;
    k1    = NULL; k2 = NULL; k3 = NULL; k4 = NULL; k5 = NULL; k6 = NULL;
    s2    = NULL; s3 = NULL; s4 = NULL; s5 = NULL; s6 = NULL;
    k_avg = NULL;
    y2    = NULL; y3 = NULL; y4 = NULL; y5 = NULL; y6 = NULL;
    next  = NULL;

    /* k1: c1 = 0 */
    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    /* k2: c2 = 1/4 */
    const derivative_weight_t s2_terms[] = { {1.0/4.0, k1} };
    s2 = derivative_average(s2_terms, state->size, 1);
    GOTO_IF_NULL(errcode, jmp_err, s2);
    y2 = euler_advance(state, s2, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y2);
    k2 = compute_system_derivative(y2, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    /* k3: c3 = 1/4 */
    const derivative_weight_t s3_terms[] = {
        {1.0/8.0, k1}, {1.0/8.0, k2}
    };
    s3 = derivative_average(s3_terms, state->size, 2);
    GOTO_IF_NULL(errcode, jmp_err, s3);
    y3 = euler_advance(state, s3, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y3);
    k3 = compute_system_derivative(y3, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k3);

    /* k4: c4 = 1/2 */
    const derivative_weight_t s4_terms[] = {
        {0.0, k1}, {0.0, k2}, {1.0/2.0, k3}
    };
    s4 = derivative_average(s4_terms, state->size, 3);
    GOTO_IF_NULL(errcode, jmp_err, s4);
    y4 = euler_advance(state, s4, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y4);
    k4 = compute_system_derivative(y4, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k4);

    /* k5: c5 = 3/4 */
    const derivative_weight_t s5_terms[] = {
        {3.0/16.0, k1}, {-3.0/8.0, k2}, {3.0/8.0, k3}, {9.0/16.0, k4}
    };
    s5 = derivative_average(s5_terms, state->size, 4);
    GOTO_IF_NULL(errcode, jmp_err, s5);
    y5 = euler_advance(state, s5, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y5);
    k5 = compute_system_derivative(y5, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k5);

    /* k6: c6 = 1 */
    const derivative_weight_t s6_terms[] = {
        {-3.0/7.0, k1}, {8.0/7.0, k2}, {6.0/7.0, k3}, {-12.0/7.0, k4}, {8.0/7.0, k5}
    };
    s6 = derivative_average(s6_terms, state->size, 5);
    GOTO_IF_NULL(errcode, jmp_err, s6);
    y6 = euler_advance(state, s6, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y6);
    k6 = compute_system_derivative(y6, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k6);

    /* final weighted average (Boole's rule weights) */
    const derivative_weight_t k_avg_terms[] = {
        {7.0/90.0,  k1}, {0.0,       k2}, {32.0/90.0, k3},
        {12.0/90.0, k4}, {32.0/90.0, k5}, {7.0/90.0,  k6}
    };
    k_avg = derivative_average(k_avg_terms, state->size, 6);
    GOTO_IF_NULL(errcode, jmp_err, k_avg);
    next = euler_advance(state, k_avg, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);
    derivative_destroy(&k3);
    derivative_destroy(&k4);
    derivative_destroy(&k5);
    derivative_destroy(&k6);
    derivative_destroy(&s2);
    derivative_destroy(&s3);
    derivative_destroy(&s4);
    derivative_destroy(&s5);
    derivative_destroy(&s6);
    derivative_destroy(&k_avg);
    fish_system_destroy(&y2);
    fish_system_destroy(&y3);
    fish_system_destroy(&y4);
    fish_system_destroy(&y5);
    fish_system_destroy(&y6);
    if (errcode != ERR_OK) {
        fish_system_destroy(&next);
        return NULL;
    }
    return next;
}

/**
 * @brief Naive 6th-order Runge-Kutta
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state \f$y(t + \delta t)\f$.
 */
fish_system_t *step_rk6(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    *error = 0;
    error_e errcode = ERR_OK;
    system_derivative_t *k1, *k2, *k3, *k4, *k5, *k6, *k7;
    system_derivative_t *s2, *s3, *s4, *s5, *s6, *s7, *k_avg;
    fish_system_t *y2, *y3, *y4, *y5, *y6, *y7, *next;
    k1    = NULL; k2 = NULL; k3 = NULL; k4 = NULL;
    k5    = NULL; k6 = NULL; k7 = NULL;
    s2    = NULL; s3 = NULL; s4 = NULL; s5 = NULL; s6 = NULL; s7 = NULL;
    k_avg = NULL;
    y2    = NULL; y3 = NULL; y4 = NULL; y5 = NULL; y6 = NULL; y7 = NULL;
    next  = NULL;

    /* k1: c1 = 0 */
    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    /* k2: c2 = 1/2 */
    const derivative_weight_t s2_terms[] = { {1.0/2.0, k1} };
    s2 = derivative_average(s2_terms, state->size, 1);
    GOTO_IF_NULL(errcode, jmp_err, s2);
    y2 = euler_advance(state, s2, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y2);
    k2 = compute_system_derivative(y2, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    /* k3: c3 = 2/3 */
    const derivative_weight_t s3_terms[] = {
        {2.0/9.0, k1}, {4.0/9.0, k2}
    };
    s3 = derivative_average(s3_terms, state->size, 2);
    GOTO_IF_NULL(errcode, jmp_err, s3);
    y3 = euler_advance(state, s3, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y3);
    k3 = compute_system_derivative(y3, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k3);

    /* k4: c4 = 1/3 */
    const derivative_weight_t s4_terms[] = {
        {7.0/36.0, k1}, {2.0/9.0, k2}, {-1.0/12.0, k3}
    };
    s4 = derivative_average(s4_terms, state->size, 3);
    GOTO_IF_NULL(errcode, jmp_err, s4);
    y4 = euler_advance(state, s4, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y4);
    k4 = compute_system_derivative(y4, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k4);

    /* k5: c5 = 5/6 */
    const derivative_weight_t s5_terms[] = {
        {-35.0/144.0, k1}, {-55.0/36.0, k2}, {35.0/48.0, k3}, {15.0/8.0, k4}
    };
    s5 = derivative_average(s5_terms, state->size, 4);
    GOTO_IF_NULL(errcode, jmp_err, s5);
    y5 = euler_advance(state, s5, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y5);
    k5 = compute_system_derivative(y5, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k5);

    /* k6: c6 = 1/6 */
    const derivative_weight_t s6_terms[] = {
        {-1.0/360.0, k1}, {-11.0/36.0, k2}, {-1.0/8.0, k3},
        {1.0/2.0, k4}, {1.0/10.0, k5}
    };
    s6 = derivative_average(s6_terms, state->size, 5);
    GOTO_IF_NULL(errcode, jmp_err, s6);
    y6 = euler_advance(state, s6, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y6);
    k6 = compute_system_derivative(y6, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k6);

    /* k7: c7 = 1 */
    const derivative_weight_t s7_terms[] = {
        {-41.0/260.0, k1}, {22.0/13.0, k2}, {43.0/156.0, k3},
        {-118.0/39.0, k4}, {32.0/195.0, k5}, {80.0/39.0, k6}
    };
    s7 = derivative_average(s7_terms, state->size, 6);
    GOTO_IF_NULL(errcode, jmp_err, s7);
    y7 = euler_advance(state, s7, time_step);
    GOTO_IF_NULL(errcode, jmp_err, y7);
    k7 = compute_system_derivative(y7, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k7);

    /* final weighted average */
    const derivative_weight_t k_avg_terms[] = {
        {13.0/200.0, k1}, {0.0,        k2}, {11.0/40.0, k3},
        {11.0/40.0,  k4}, {4.0/25.0,   k5}, {4.0/25.0,  k6}, {13.0/200.0, k7}
    };
    k_avg = derivative_average(k_avg_terms, state->size, 7);
    GOTO_IF_NULL(errcode, jmp_err, k_avg);
    next = euler_advance(state, k_avg, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);
    derivative_destroy(&k3);
    derivative_destroy(&k4);
    derivative_destroy(&k5);
    derivative_destroy(&k6);
    derivative_destroy(&k7);
    derivative_destroy(&s2);
    derivative_destroy(&s3);
    derivative_destroy(&s4);
    derivative_destroy(&s5);
    derivative_destroy(&s6);
    derivative_destroy(&s7);
    derivative_destroy(&k_avg);
    fish_system_destroy(&y2);
    fish_system_destroy(&y3);
    fish_system_destroy(&y4);
    fish_system_destroy(&y5);
    fish_system_destroy(&y6);
    fish_system_destroy(&y7);
    if (errcode != ERR_OK) {
        fish_system_destroy(&next);
        return NULL;
    }
    return next;
}

/**
 * @brief 3rd-order Runge-Kutta (Bogacki-Shampine) with 2nd order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk32(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    error_e errcode = ERR_OK;

    system_derivative_t *k1, *k2, *k3, *k4, *k_avg_2, *k_avg_3;
    fish_system_t *y1, *y2, *next_2, *next_3, *diff;

    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    y1 = euler_advance(state, k1, time_step / 2.0);
    GOTO_IF_NULL(errcode, jmp_err, y1);

    k2 = compute_system_derivative(y1, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    y2 = euler_advance(state, k2, time_step * 3.0 / 4);
    GOTO_IF_NULL(errcode, jmp_err, y2);

    k3 = compute_system_derivative(y2, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k3);

    const derivative_weight_t k_avg_3_terms[] = {
        {2.0f/9,  k1}, {1.0f/3, k2}, {4.0f/9, k3}
    };
    k_avg_3 = derivative_average(k_avg_3_terms, state->size, 3);
    GOTO_IF_NULL(errcode, jmp_err, k_avg_3);

    next_3 = euler_advance(state, k_avg_3, time_step);
    GOTO_IF_NULL(errcode, jmp_err, k_avg_3);

    k4 = compute_system_derivative(next_3, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k4);

    const derivative_weight_t k_avg_2_terms[] = {
        {7.0f/24, k1}, {1.0f/4, k2}, {1.0f/3, k3}, {1.0f/8, k4}
    };
    k_avg_2 = derivative_average(k_avg_2_terms, state->size, 4);
    GOTO_IF_NULL(errcode, jmp_err, k_avg_2);

    next_2 = euler_advance(state, k_avg_2, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next_2);

    diff = fish_system_difference(next_3, next_2);
    *error = fish_system_norm(diff);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);
    derivative_destroy(&k3);
    derivative_destroy(&k4);
    derivative_destroy(&k_avg_2);
    derivative_destroy(&k_avg_3);

    fish_system_destroy(&y1);
    fish_system_destroy(&y2);
    fish_system_destroy(&next_2);
    fish_system_destroy(&diff);

    if (errcode != ERR_OK) {
        fish_system_destroy(&next_3);
        return NULL;
    }

    return next_3;
}

/**
 * @brief 5th-order Runge-Kutta (Dormand-Prince) with 4th order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk54(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{
    error_e errcode = ERR_OK;

    system_derivative_t *k1, *k2, *k3, *k4, *k5, *k6, *k7, *k_avg_4, *k_avg_5;
    fish_system_t *y1, *y2, *y3, *y4, *y5, *next_4, *next_5, *diff;

    /* Stage 1 */
    k1 = compute_system_derivative(state, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k1);

    /* Stage 2 */
    {
        const derivative_weight_t terms[] = {
            {1.0/5.0, k1}
        };
        system_derivative_t *avg = derivative_average(terms, state->size, 1);
        GOTO_IF_NULL(errcode, jmp_err, avg);

        y1 = euler_advance(state, avg, time_step);
        derivative_destroy(&avg);
    }
    GOTO_IF_NULL(errcode, jmp_err, y1);

    k2 = compute_system_derivative(y1, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k2);

    /* Stage 3 */
    {
        const derivative_weight_t terms[] = {
            { 3.0/40.0, k1 },
            { 9.0/40.0, k2 }
        };
        system_derivative_t *avg = derivative_average(terms, state->size, 2);
        GOTO_IF_NULL(errcode, jmp_err, avg);

        y2 = euler_advance(state, avg, time_step * 3.0/10.0);
        derivative_destroy(&avg);
    }
    GOTO_IF_NULL(errcode, jmp_err, y2);

    k3 = compute_system_derivative(y2, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k3);

    /* Stage 4 */
    {
        const derivative_weight_t terms[] = {
            {  44.0/45.0, k1 },
            { -56.0/15.0, k2 },
            {  32.0/9.0,  k3 }
        };
        system_derivative_t *avg = derivative_average(terms, state->size, 3);
        GOTO_IF_NULL(errcode, jmp_err, avg);

        y3 = euler_advance(state, avg, time_step * 4.0/5.0);
        derivative_destroy(&avg);
    }
    GOTO_IF_NULL(errcode, jmp_err, y3);

    k4 = compute_system_derivative(y3, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k4);

    /* Stage 5 */
    {
        const derivative_weight_t terms[] = {
            {  19372.0/6561.0, k1 },
            { -25360.0/2187.0, k2 },
            {  64448.0/6561.0, k3 },
            { -212.0  /729.0,  k4 }
        };
        system_derivative_t *avg = derivative_average(terms, state->size, 4);
        GOTO_IF_NULL(errcode, jmp_err, avg);

        y4 = euler_advance(state, avg, time_step * 8.0/9.0);
        derivative_destroy(&avg);
    }
    GOTO_IF_NULL(errcode, jmp_err, y4);

    k5 = compute_system_derivative(y4, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k5);

    /* Stage 6 */
    {
        const derivative_weight_t terms[] = {
            {  9017.0/3168.0,  k1 },
            { -355.0/33.0,     k2 },
            {  46732.0/5247.0, k3 },
            {  49.0/176.0,     k4 },
            { -5103.0/18656.0, k5 }
        };
        system_derivative_t *avg = derivative_average(terms, state->size, 5);
        GOTO_IF_NULL(errcode, jmp_err, avg);

        y5 = euler_advance(state, avg, time_step);
        derivative_destroy(&avg);
    }
    GOTO_IF_NULL(errcode, jmp_err, y5);

    k6 = compute_system_derivative(y5, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k6);

    /* Fifth-order solution */
    {
        const derivative_weight_t terms[] = {
            {  35.0/384.0,    k1 },
            {  500.0/1113.0,  k3 },
            {  125.0/192.0,   k4 },
            { -2187.0/6784.0, k5 },
            {  11.0/84.0,     k6 }
        };

        k_avg_5 = derivative_average(terms, state->size, 5);
    }
    GOTO_IF_NULL(errcode, jmp_err, k_avg_5);

    next_5 = euler_advance(state, k_avg_5, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next_5);

    /* FSAL stage */
    k7 = compute_system_derivative(next_5, dc_opts);
    GOTO_IF_NULL(errcode, jmp_err, k7);

    /* Embedded fourth-order solution */
    {
        const derivative_weight_t terms[] = {
            { 5179.0/57600.0,   k1 },
            { 7571.0/16695.0,   k3 },
            { 393.0/640.0,      k4 },
            {-92097.0/339200.0, k5 },
            { 187.0/2100.0,     k6 },
            { 1.0/40.0,         k7 }
        };

        k_avg_4 = derivative_average(terms, state->size, 6);
    }
    GOTO_IF_NULL(errcode, jmp_err, k_avg_4);

    next_4 = euler_advance(state, k_avg_4, time_step);
    GOTO_IF_NULL(errcode, jmp_err, next_4);

    diff = fish_system_difference(next_5, next_4);
    *error = fish_system_norm(diff);

jmp_err:
    derivative_destroy(&k1);
    derivative_destroy(&k2);
    derivative_destroy(&k3);
    derivative_destroy(&k4);
    derivative_destroy(&k5);
    derivative_destroy(&k6);
    derivative_destroy(&k7);
    derivative_destroy(&k_avg_4);
    derivative_destroy(&k_avg_5);

    fish_system_destroy(&y1);
    fish_system_destroy(&y2);
    fish_system_destroy(&y3);
    fish_system_destroy(&y4);
    fish_system_destroy(&y5);
    fish_system_destroy(&next_4);
    fish_system_destroy(&diff);

    if (errcode != ERR_OK) {
        fish_system_destroy(&next_5);
        return NULL;
    }

    return next_5;
}
