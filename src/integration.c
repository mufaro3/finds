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
    { RUNGE_KUTTA_23, "RK23"  },
    { RUNGE_KUTTA_45, "RK45"  },
    { RUNGE_KUTTA_54, "RK54"  },
    { RUNGE_KUTTA_65, "RK65"  },
    { RUNGE_KUTTA_78, "RK78"  },
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
 * @brief 2nd-order Runge-Kutta with 3rd order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk23(
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
 * @brief 4nd-order Runge-Kutta with 5th order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk45(
    UNUSED const fish_system_t *state,
    UNUSED const double time_step,
    UNUSED const derivative_computation_opts_t dc_opts,
    UNUSED double *error)
{
    NOT_IMPLEMENTED();
}

/**
 * @brief 5th-order Runge-Kutta with 4th order error estimation
 *        (Dormand-Prince Method).
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk54(
    UNUSED const fish_system_t *state,
    UNUSED const double time_step,
    UNUSED const derivative_computation_opts_t dc_opts,
    UNUSED double *error)
{
    NOT_IMPLEMENTED();
}

/**
 * @brief 6th-order Runge-Kutta with 5th order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk65(
    UNUSED const fish_system_t *state,
    UNUSED const double time_step,
    UNUSED const derivative_computation_opts_t dc_opts,
    UNUSED double *error)
{
    NOT_IMPLEMENTED();
}

/**
 * @brief 7th-order Runge-Kutta with 8th order error estimation.
 *
 * @param[in]  state     The system state for this time.
 * @param[in]  time_step The time-step to advance to.
 * @param[out] error     The RMS error-norm for this state.
 *
 * @return The advanced state y(t + delta t).
 */
fish_system_t *step_rk78(
    UNUSED const fish_system_t *state,
    UNUSED const double time_step,
    UNUSED const derivative_computation_opts_t dc_opts,
    UNUSED double *error)
{
    NOT_IMPLEMENTED();
}
