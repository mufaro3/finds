#include "derivative.h"
#include "integration.h"

/*
 * \brief Computes y + y' * dt.
 *
 * \param[in] state      y
 * \param[in] derivative y'
 * \param[in] time_step  dt
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

/*
 * \brief Naive Euler-Stepping without Error Estimation
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).

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


/*
 * \brief 2nd-order Runge-Kutta with 3rd order error estimation.
 *
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).
 */
fish_system_t *step_rk23(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{}

/*
 * \brief 4nd-order Runge-Kutta with 5th order error estimation.
 *
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).
 */
fish_system_t *step_rk45(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{}

/*
 * \brief 5th-order Runge-Kutta with 4th order error estimation
 *        (Dormand-Prince Method).
 *
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).
 */
fish_system_t *step_rk54(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{}

/*
 * \brief 6th-order Runge-Kutta with 5th order error estimation.
 *
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).
 */
fish_system_t *step_rk65(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{}

/*
 * \brief 7th-order Runge-Kutta with 8th order error estimation.
 *
 * \param[in]  state     The system state for this time.
 * \param[in]  time_step The time-step to advance to.
 * \param[out] error     The RMS error-norm for this state.
 *
 * \return The advanced state y(t + delta t).
 */
fish_system_t *step_rk78(
    const fish_system_t *state,
    const double time_step,
    const derivative_computation_opts_t dc_opts,
    double *error)
{}
