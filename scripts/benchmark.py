# disable any numba speedhacking
import os
os.environ["NUMBA_DISABLE_JIT"] = "1"

import time
from typing import Optional

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from finds.calculations import calculate_system_derivative
from finds.constants import BENCHMARK_OUTPUT_PATH
from finds.fish import generate_system

COMPARISON_FIGURE_OUTPUT_NAME='comparison.png'


def perform_time_test(n: int,
                      bh_ratio: Optional[float] = None,
                      repeats: int = 3) -> float:
    r"""
    Times a singular calculation for the state derivative of a random
    system of size :math:`n`.

    :param n: The system size.
    :type  n: int

    :param bh_ratio: The Barnes-Hut Ratio, if Barnes-Hut approximation is
      to be used.
    :type  bh_ratio: Optional[float]

    :param repeats: The number of times to repeat the experiment for this
      configuration (default 3).
    :type  repeats: int

    :returns: The mean runtime for computing :math:`d\mathbf{X}/dt` in
      milliseconds.
    :rtype: float
    """
    # generate the random system
    example_system = generate_system(
        distribution='random',
        orientation='random',
        n_random=n,
        debug_print=False
    )

    use_barnes_hut = bh_ratio is not None

    # a warmup to exclude numba compilation times
    calculate_system_derivative(
        example_system,
        use_barnes_hut,
        bh_ratio if bh_ratio is not None else 0.0
    )

    times = np.zeros(repeats)

    for i in range(repeats):
        # begin the timer
        start_time = time.perf_counter()

        # compute the derivative of the system matrix
        calculate_system_derivative(
            example_system,
            use_barnes_hut,
            bh_ratio if bh_ratio is not None else 0.0
        )

        # end the timer
        end_time = time.perf_counter()
        times[i] = end_time - start_time

    # return only the mean
    return np.mean(times)


def generate_comparison_figure(
        min_log_n: int,
        max_log_n: int) -> None:
    r"""
    Generates the timing comparison between Brute-Force and Barnes-Hut.

    :param min_log_n: The starting :math:`\log_2 N` (default 1).
    :type  min_log_n: int

    :param max_log_n: The ending :math:`\log_2 N` (default 5).
    :type  max_log_n: int

    The relationship between the input size :math:`N` and the output time
    :math:`t` is defined as

    .. math::
        t_{BF} = \mathcal{O}(N^2) = k_{BF} N^2

    for brute-force and

    .. math::
        t_{BH} = \mathcal{O}(N \log_2 N) = k_{BH} N \log_2 N

    for the Barnes-Hut approximation. The proportionality constants,
    :math:`k_{BF}` and :math:`k_{BH}`, represent the inherent, non-algorithmic
    system load from how I've written the program across the two systems, and
    as a result of performing object-based recursion for Barnes-Hut and
    utilizing Numba (and parallelization) for Brute-Force, we can be fairly
    confident that :math:`k_{BF} << k_{BH}` (which makes this comparison
    somewhat difficult at low :math:`N`, as the proportionality coefficients
    will dominate there).

    Given these relationships, the produced plot is made logarithmic, i.e.,

    .. math::
        \log_2 t_{BF} = \log_2 k_{BF} + 2 \log_2 N

    for brute-force and

    .. math::
        \log_2 t_{BH} = \log_2 k_{BH} + \log_2 N + \log_2 \log_2 N

    for the barnes-hut approximation. These functions are essentially linear,
    given :math:`y := \log_2 t, b := \log_2 k`, :math:`x := \log_2 N`, and the
    fact that :math:`\log_2 \log_2 N` is essentially constant at high
    :math:`N`, thus we can perform linear regression to compute the
    proportionality constants as the :math:`y`-intercept, and validate the
    functional forms of each.

    In addition to this form, the benchmarker also produces a plot of
    :math:`k` as a function of :math:`N` for both Brute-Force and Barnes-Hut
    by normalizing the non-logarithmic runtime by the algorithmic growth
    function. For Brute-Force, this is

    .. math::
        y := k_{BF} = \frac{t_{BF}}{N^2}

    and for Barnes-Hut, this is

    .. math::
        y := k_{BH} = \frac{t_{BH}}{N \log_2 N},

    and this should produce either a flat/constant or an asymptotic curve as
    :math:`k` should be roughly constant (and ideally, quite small). This acts
    as a verification measure, because if this curve shows that :math:`k` is
    a function of :math:`N` in a nonnegligible way, then the algorithm is
    likely contains a bug.
    """
    tqdm.write('Beginning benchmark.')

    # Generate the system sizes from the input
    lognvalues = np.flip(np.arange(min_log_n, max_log_n + 1))
    nvalues = np.exp2(lognvalues)

    # Barnes-Hut opening angles
    bh_ratios = np.arange(0.25, 1.0, step=0.25)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))

    log_plot_ax = axes[0]
    k_asympt_ax = axes[1]

    # brute-force baseline
    brute_force_times = []

    for logn in tqdm(lognvalues, desc='Brute Force Benchmark'):
        time=perform_time_test(2 ** logn, bh_ratio=None)
        brute_force_times.append(np.log2(time))

    slope_bf, intercept_bf = np.polyfit(lognvalues, brute_force_times, 1)
    prop_const_bf = 2 **  intercept_bf

    # formats exponents
    def fe(x):
        mantissa, exponent = f"{x:.6e}".split("e")
        return f"{float(mantissa):.2f} \\times 10^{{{int(exponent)}}}"

    log_plot_ax.plot(
        lognvalues, brute_force_times,
        label=f"Brute Force, $k={fe(prop_const_bf)}$, slope={slope_bf:.1f}",
        linewidth=3
    )

    k_asympt_ax.plot(
        nvalues,
        np.exp2(brute_force_times) / (nvalues ** 2),
        label=r'Brute Force, $T(N)/N^2$'
    )

    # Barnes-Hut curves
    for bh_ratio in tqdm(bh_ratios, leave=False, desc='Barnes-Hut'):

        # compute the times
        bh_times = []
        for logn in tqdm(lognvalues, desc=f'Barnes-Hut, Theta={bh_ratio:.2f}'):
            time = perform_time_test(2 ** logn, bh_ratio=bh_ratio)
            bh_times.append(np.log2(time))

        slope_bh, intercept_bh = np.polyfit(lognvalues, bh_times, 1)
        prop_const_bh = 2 ** intercept_bh

        # plot
        log_plot_ax.plot(
            lognvalues, bh_times,
            label=f'Barnes-Hut, $\\theta={bh_ratio:.2f}$, '+\
            f'$k={fe(prop_const_bh)}$, Slope$={slope_bh:.1f}$',
            linewidth=3
        )

        k_asympt_ax.plot(
            nvalues,
            np.exp2(bh_times) / (nvalues * np.log2(nvalues)),
            label=f'Barnes-Hut, $\\theta={bh_ratio:.2f}$, $T(N)/(N \\log_2 N)$'
        )

    # set labels and title
    log_plot_ax.set_xlabel(r"Logarithmic Number of Fish, $\log_2 N$")
    log_plot_ax.set_ylabel(r"Logarithmic Runtime ($\log_2$ seconds)")
    log_plot_ax.set_title("Barnes–Hut vs Brute Force\n"+\
                          r"$\dot{\mathbf{X}}$ Computation Runtime")
    log_plot_ax.legend()
    log_plot_ax.grid(True)

    k_asympt_ax.set_xlabel(r'Number of Fish, $N$')
    k_asympt_ax.set_ylabel(
        r'Proportionality Cofficient $k$\\(seconds per fish)')
    k_asympt_ax.set_title('Normalized Algorithmic Runtime Growth vs. $N$')
    k_asympt_ax.legend()
    k_asympt_ax.grid(True)

    BENCHMARK_OUTPUT_PATH.mkdir(parents=True, exist_ok=True)
    output_path = BENCHMARK_OUTPUT_PATH / COMPARISON_FIGURE_OUTPUT_NAME
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    tqdm.write(f"Saved comparison figure to {output_path}")


if __name__ == '__main__':
    generate_comparison_figure(1, 6)
