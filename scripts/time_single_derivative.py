import argparse
import time
from tqdm import tqdm

from finds.fish import generate_system
from finds.calculations import calculate_system_derivative


def format_duration(seconds: float) -> str:
    seconds = int(round(seconds))

    days, seconds = divmod(seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)

    parts = []

    if days:
        parts.append(f"{days} day{'s' if days != 1 else ''}")
    if hours:
        parts.append(f"{hours} hour{'s' if hours != 1 else ''}")
    if minutes:
        parts.append(f"{minutes} minute{'s' if minutes != 1 else ''}")

    parts.append(f"{seconds} second{'s' if seconds != 1 else ''}")

    return ", ".join(parts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Benchmark derivative calculation."
    )
    parser.add_argument(
        "N",
        type=int,
        help="Number of fish to generate.",
    )
    parser.add_argument(
        "--theta",
        type=float,
        default=None,
        help="Barnes-Hut theta. If omitted or 0, Barnes-Hut is disabled.",
    )

    args = parser.parse_args()

    use_barnes_hut = args.theta is not None and args.theta > 0.0

    if use_barnes_hut:
        tqdm.write(f"Using Barnes-Hut (theta = {args.theta})")
    else:
        tqdm.write("Using direct calculation (Barnes-Hut disabled)")

    system = generate_system(
        distribution="random",
        orientation="random",
        n_random=args.N,
        debug_print=True,
    )

    start_time = time.perf_counter()

    calculate_system_derivative(
        system,
        use_barnes_hut=use_barnes_hut,
        bh_ratio=args.theta if use_barnes_hut else 0.0,
        show_progress=True,
    )

    elapsed = time.perf_counter() - start_time

    tqdm.write(f"Finished in {format_duration(elapsed)} ({elapsed:.3f} s)")
