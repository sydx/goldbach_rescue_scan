#!/usr/bin/env python3
"""Reproduce the finite computations in Twin--Cousin Goldbach Representations.

Python 3.9+; standard library only. This is a finite verification,
not a proof of either universal conjecture.
"""
import argparse
from collections import Counter
from statistics import mean


def sieve(limit):
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0:2] = b'\x00\x00'
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p*p::p] = b'\x00' * ((limit - p*p) // p + 1)
    return is_prime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=int, default=15000,
                        help='Largest cumulative stage (default: 15000)')
    parser.add_argument('--blocks', type=int, default=5000,
                        help='Largest certified block index (default: 5000)')
    args = parser.parse_args()
    assert args.stage >= 7 and args.blocks >= 2

    limit = max(2 * args.stage + 10, 6 * (args.blocks - 1) + 5)
    is_prime = sieve(limit)
    twins = [a for a in range(1, (limit - 1) // 6 + 1)
             if is_prime[6*a - 1] and is_prime[6*a + 1]]
    cousins = [b for b in range(1, (limit - 5) // 6 + 1)
               if is_prime[6*b + 1] and is_prime[6*b + 5]]

    # Exact certificate counts and earliest arrival times.
    max_m = twins[-1] + cousins[-1] + 2
    INF = 10**15
    counts = [0] * (max_m + 1)
    arrivals = [INF] * (max_m + 1)
    first_witness = [None] * (max_m + 1)
    for a in twins:
        for b in cousins:
            m = a + b
            counts[m] += 1
            arrival = max(3*a + 2, 3*b + 4)
            if arrival < arrivals[m]:
                arrivals[m] = arrival
                first_witness[m] = (a, b)

    # Check via a second, direct scan of all candidate a at each block index.
    for m in range(2, args.blocks + 1):
        direct_count = sum(
            bool(is_prime[6*a-1] and is_prime[6*a+1]
                 and is_prime[6*(m-a)+1] and is_prime[6*(m-a)+5])
            for a in range(1, m)
        )
        assert direct_count == counts[m], (m, direct_count, counts[m])
        assert direct_count > 0, f'TCB failed at index {m}'
        assert arrivals[m] <= 3*m + 1, ('Deadline failed', m)

    selected = {22, 23, 25, 32, 37, 38, 40, 43, 55, 100, 500,
                1000, 3000, 5000, 10000, 15000}
    occupancy = Counter()
    frontier = 2
    minimum_buffer = (INF, -1)
    print('stage n | first uncertified m | 6m+2 | buffer D')
    for n in range(3, args.stage + 1):
        while arrivals[frontier] <= n:
            frontier += 1
        buffer = 3*frontier + 1 - n
        if n >= 7:
            assert buffer > 0, (n, frontier, buffer)
            minimum_buffer = min(minimum_buffer, (buffer, n))
        if 50 <= n <= 10000:
            occupancy[frontier % 5] += 1
        if n in selected:
            print(f'{n:>7} | {frontier:>19} | {6*frontier+2:>4} | {buffer:>8}')

    print(f'\nVerified every regular block 2 <= m <= {args.blocks}')
    print(f'Checked positive rescue buffer through stage n={args.stage}')
    print('Minimum buffer (n>=7):', minimum_buffer)
    if args.blocks >= 5000:
        print('\nMean R(m) over 50 <= m <= 5000 by m mod 5:')
        for r in range(5):
            print(f'{r}: {mean(counts[m] for m in range(50,5001) if m%5==r):.3f}')
    if args.stage >= 10000:
        print('Frontier occupancy for 50 <= n <= 10000: ',
              {r: occupancy[r] for r in range(5)})

    # Independently reconstruct all individual twin/cousin-prime sums at a stage.
    odd_primes = [p for p in range(3, 2*args.stage, 2) if is_prime[p]]
    for n in [22, 23, 25, 32, 37, 38, 40]:
        if n > args.stage:
            continue
        cutoff = 2*n - 3
        T = [p for p in odd_primes if p <= cutoff and
             ((p >= 5 and is_prime[p-2]) or
              (p+2 <= cutoff and is_prime[p+2]))]
        C = [p for p in odd_primes if p <= cutoff and
             ((p >= 7 and is_prime[p-4]) or
              (p+4 <= cutoff and is_prime[p+4]))]
        sums = {t+c for t in T for c in C}
        maximum = max(sums)
        gaps = [v for v in range(6, maximum + 1, 2) if v not in sums]
        print(f'n={n}: maximum={maximum}, actual gaps={gaps}')


if __name__ == '__main__':
    main()
