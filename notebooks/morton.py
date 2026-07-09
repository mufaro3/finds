# Source - https://stackoverflow.com/a/18528775
# Posted by Gabriel, modified by community. See post 'Timeline' for change history.
# Retrieved 2026-07-08, License - CC BY-SA 4.0

from __future__ import annotations

import random


def pretty_bin_string(
    x: int,
    width: int = 32,
    group: int = 4,
    sep: str = ".",
    fill: str = "0",
) -> str:
    """Return a grouped binary string."""

    bits = format(x, "b")

    zeros = width - len(bits)

    if zeros <= 0:
        zeros = 0
        k = group - (len(bits) % group)
    else:
        k = group - (width % group)

    out = []

    for i in range(zeros):
        if k % group == 0 and i != 0:
            out.append(sep)
        out.append(fill)
        k += 1

    for i, bit in enumerate(bits):
        if (
            (k % group == 0 and i != 0 and zeros == 0)
            or (k % group == 0 and zeros != 0)
        ):
            out.append(sep)
        out.append(bit)
        k += 1

    return "".join(out)


def bin_str(x: int) -> str:
    return pretty_bin_string(x, 32, 4, " ", "0")


def compute_bit_mask_pattern_and_code(
    number_of_bits: int,
    number_of_empty_bits: int,
) -> str:
    bit_distances = [
        i * number_of_empty_bits
        for i in range(number_of_bits)
    ]
    print(f"Bit Distances: {bit_distances}")

    bit_distances_binary = [
        format(dist, "b")
        for dist in bit_distances
    ]
    print(f"Bit Distances (binary): {bit_distances_binary}")

    move_bits: list[list[int]] = []

    max_length = max(map(len, bit_distances_binary))

    for i in range(max_length):
        move_bits.append([])

        for idx, bits in enumerate(bit_distances_binary):
            if len(bits) - 1 >= i and bits[-i - 1] == "1":
                move_bits[i].append(idx)

    for i, bits in enumerate(move_bits):
        print(f"Shifting bits by {2**i}\t for bits idx: {bits}")

    bit_positions = list(range(number_of_bits))
    print(f"BitPositions: {bit_positions}")

    mask_old = (1 << number_of_bits) - 1

    code_lines = [
        f"x &= {hex(mask_old)}"
    ]

    for idx in range(len(move_bits) - 1, -1, -1):
        if not move_bits[idx]:
            continue

        shifted = 0

        for bit_idx in move_bits[idx]:
            shifted |= 1 << bit_positions[bit_idx]
            bit_positions[bit_idx] += 2**idx

        nonshifted = (~shifted) & mask_old

        print(
            f"Shifted bef.:\t{bin_str(shifted)} "
            f"hex: {hex(shifted)}"
        )

        shifted <<= 2**idx

        print(
            f"Shifted:\t{bin_str(shifted)} "
            f"hex: {hex(shifted)}"
        )

        print(
            f"NonShifted:\t{bin_str(nonshifted)} "
            f"hex: {hex(nonshifted)}"
        )

        mask_new = shifted | nonshifted

        print(
            f"Bitmask is now:\t{bin_str(mask_new)} "
            f"hex: {hex(mask_new)}\n"
        )

        code_lines.append(
            f"x = (x | (x << {2**idx})) & {hex(mask_new)}"
        )

        mask_old = mask_new

    return "\n".join(code_lines)


number_of_bits = 10
number_of_empty_bits = 2

code_string = compute_bit_mask_pattern_and_code(
    number_of_bits,
    number_of_empty_bits,
)

print(code_string)


def partition_by_2(x: int) -> int:
    namespace = {"x": x}
    exec(code_string, {}, namespace)
    return namespace["x"]


def check_partition(x: int) -> bool:
    print(f"Check partition for:\t{bin_str(x)}")

    part = partition_by_2(x)
    print(f"Partition is:\t\t{bin_str(part)}")

    # Construct the partition manually.
    part_check = 0

    for bit_idx in range(number_of_bits):
        part_check |= (
            (x & (1 << bit_idx))
            << (number_of_empty_bits * bit_idx)
        )

    print(f"Partition check is:\t{bin_str(part_check)}")

    return part_check == part


check_error = False

for _ in range(20):
    x = random.getrandbits(number_of_bits)

    if not check_partition(x):
        check_error = True
        break

if not check_error:
    print("CHECK PARTITION SUCCESSFUL!!!!!!!!!!!!!!!!...")
else:
    print("checkPartition has ERROR!!!!")
