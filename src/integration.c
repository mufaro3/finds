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
    const system_derivative_t *derivative = \
        compute_system_derivative(state, dc_opts);
    *error = 0;
    return euler_advance(state, derivative, time_step);
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
    UNUSED const fish_system_t *state,
    UNUSED const double time_step,
    UNUSED const derivative_computation_opts_t dc_opts,
    UNUSED double *error)
{
    NOT_IMPLEMENTED();
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
