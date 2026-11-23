"""Enumerate a small CSS cascade to audit structured stabilizer membership."""

from itertools import product


def span(rows):
    result = {0}
    for row in rows:
        result |= {value ^ row for value in tuple(result)}
    return result


def x_check(error):
    blocks = [error[6*j:6*j+6] for j in range(2)]
    if any(block[0] != (sum(block[1:]) & 1) for block in blocks):
        return False
    return all(blocks[0][a+1] == blocks[1][a+1] for a in range(5))


def z_check(error):
    blocks = [error[6*j:6*j+6] for j in range(2)]
    return all((blocks[0][a] ^ blocks[0][5]) ==
               (blocks[1][a] ^ blocks[1][5]) for a in range(5))


# H1=[1,1], G2=[1|I5], selector=[I5|0], H2=all ones.
hx_rows = []
hz_rows = []
for a in range(5):
    hx_rows.append(sum(1 << (6*j+c) for j in range(2) for c in (0,a+1)))
    hz_rows.append(sum(1 << (6*j+a) for j in range(2)))
hz_rows.extend(sum(1 << (6*j+c) for c in range(6)) for j in range(2))
hx_span = span(hx_rows)
hz_span = span(hz_rows)

for value in range(1 << 12):
    bits = [(value >> bit) & 1 for bit in range(12)]
    assert x_check(bits) == (value in hx_span), value
    assert z_check(bits) == (value in hz_span), value

print("structured stabilizer audit: 4096/4096 vectors agree for X and Z")
