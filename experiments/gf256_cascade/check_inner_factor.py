"""Enumerate a small inner factor to audit the structured decoder equation."""

from itertools import product
from random import Random


def normalized(values):
    total = sum(values)
    return [value / total for value in values]


def local_message_bruteforce(channel, incoming, target):
    q = len(channel[0])
    outgoing = [0.0] * q
    for tail in range(q):
        for projected in product(range(q), repeat=len(incoming)):
            weight = channel[0][tail]
            for a, value in enumerate(projected):
                weight *= channel[a + 1][tail ^ value]
                if a != target:
                    weight *= incoming[a][value]
            outgoing[projected[target]] += weight
    return normalized(outgoing)


def local_message_factorized(channel, incoming, target):
    q = len(channel[0])
    weight = channel[0][:]
    for a in range(len(incoming)):
        if a == target:
            continue
        for tail in range(q):
            weight[tail] *= sum(
                incoming[a][value] * channel[a + 1][tail ^ value]
                for value in range(q)
            )
    outgoing = [
        sum(weight[tail] * channel[target + 1][tail ^ value]
            for tail in range(q))
        for value in range(q)
    ]
    return normalized(outgoing)


rng = Random(20261022)
q = 4
channel = [normalized([rng.random() + 0.1 for _ in range(q)]) for _ in range(6)]
incoming = [normalized([rng.random() + 0.1 for _ in range(q)]) for _ in range(5)]
maximum_difference = 0.0
for target in range(5):
    brute = local_message_bruteforce(channel, incoming, target)
    factored = local_message_factorized(channel, incoming, target)
    maximum_difference = max(
        maximum_difference,
        *(abs(x - y) for x, y in zip(brute, factored))
    )
assert maximum_difference < 1e-12, maximum_difference
print(f"inner factor enumeration: all 5 messages match, max difference={maximum_difference:.3g}")
