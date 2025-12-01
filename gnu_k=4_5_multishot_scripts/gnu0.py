import collections

print(">>>>")

def is_gnu4(factorization: dict) -> bool:
    """
    Checks if the number of groups of order n is 4.
    The conditions are derived from a paper by G. A. Miller.
    Factorization is a dict {p1: e1, p2: e2, ...}.
    """
    primes = list(factorization.keys())
    exponents = list(factorization.values())

    # Helper to count congruences p % q == 1
    def get_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and p1 % p2 == 1:
                    graph[p1].add(p2)
        return graph

    # Case 1, 2, 3: n is square-free
    if all(e == 1 for e in exponents):
        graph = get_congruence_graph(primes)

        # Case 1 (sq-free, odd)
        if 2 not in primes:
            counts = {p: len(targets) for p, targets in graph.items()}
            primes_with_2_congs = [p for p, count in counts.items() if count == 2]
            primes_with_gt2_congs = [p for p, count in counts.items() if count > 2]
            if len(primes_with_2_congs) == 1 and not primes_with_gt2_congs:
                return True

        # Case 2 (sq-free, odd)
        if 2 not in primes:
            counts = {p: len(targets) for p, targets in graph.items()}
            primes_with_1_cong = [p for p, count in counts.items() if count == 1]
            primes_with_gt1_cong = [p for p, count in counts.items() if count > 1]
            if len(primes_with_1_cong) == 2 and not primes_with_gt1_cong:
                p1, p2 = primes_with_1_cong
                q1 = list(graph[p1])[0]
                q2 = list(graph[p2])[0]
                if q1 != q2:
                    return True

        # Case 3 (sq-free, even, n=2pq)
        if factorization.get(2, 0) == 1 and len(primes) == 3:
            odd_primes = [p for p in primes if p != 2]
            p, q = min(odd_primes), max(odd_primes)
            if q % p != 1:
                return True

        return False

    # Cases with square factors
    primes_with_sq = [p for p, e in factorization.items() if e >= 2]

    # Case 4: n = p1^2 * p2^2 * ..., no congruences
    if len(primes_with_sq) >= 2 and all(factorization[p] == 2 for p in primes_with_sq) and all(e == 1 for p, e in factorization.items() if p not in primes_with_sq):
        graph = get_congruence_graph(primes)
        if not any(graph.values()):
            return True

    # Cases with exactly one square factor, p1^2
    if len(primes_with_sq) == 1 and factorization[primes_with_sq[0]] == 2 and all(e == 1 for p, e in factorization.items() if p not in primes_with_sq):
        p1 = primes_with_sq[0]
        p_rest = [p for p in primes if p != p1]

        # Case 5: n = p1^2 * q * ..., q % p1 == 1 is primary feature
        congs_to_p1 = [q for q in p_rest if q % p1 == 1]
        if len(congs_to_p1) == 1:
            q_special = congs_to_p1[0]
            if q_special % (p1 * p1) != 1:
                graph_rest = get_congruence_graph(p_rest)
                if not any(graph_rest.values()):
                    if not any((p1 + 1) % p == 0 for p in p_rest):
                        return True

        # Case 6: n = p1^2 * ..., exactly one congruence p%q=1, where q!=p1
        graph = get_congruence_graph(primes)
        total_congs = sum(len(v) for v in graph.values())
        if total_congs == 1:
            p_cong = [p for p, targets in graph.items() if targets][0]
            q_cong = list(graph[p_cong])[0]
            if q_cong != p1:
                if not any(p % p1 == 1 for p in p_rest):
                    if not any((p1 + 1) % p == 0 for p in p_rest):
                        return True

    return False

def is_gnu5(factorization: dict) -> bool:
    """
    Checks if the number of groups of order n is 5.
    The conditions are derived from a paper by G. A. Miller.
    Factorization is a dict {p1: e1, p2: e2, ...}.
    """
    if factorization == {2: 2, 3: 1}: # n=12 is a special case
        return True

    primes = list(factorization.keys())
    exponents = list(factorization.values())

    # Helper to get inverse congruence graph {q: {p | p%q==1}}
    def get_inverse_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and p1 % p2 == 1:
                    graph[p2].add(p1)
        return graph

    # Square-free cases
    if all(e == 1 for e in exponents):
        if 2 in primes: return False
        inv_graph = get_inverse_congruence_graph(primes)
        total_congs = sum(len(v) for v in inv_graph.values())

        # Case SF1: n=3qr..., exactly two primes are 1 mod 3, no other congruences
        if 3 in primes and len(inv_graph.get(3, [])) == 2 and total_congs == 2:
            return True

        # Case SF2/3: n=p1 p2 p3 p4..., p2%p1=1, p3%p1=1, and (p2%p4=1 or p3%p4=1)
        if len(primes) >= 4:
            for p1 in primes:
                if len(inv_graph.get(p1, [])) == 2:
                    p2, p3 = tuple(inv_graph[p1])
                    p_rest = [p for p in primes if p not in [p1, p2, p3]]
                    for p4 in p_rest:
                        is_p2_mod_p4 = p2 % p4 == 1
                        is_p3_mod_p4 = p3 % p4 == 1
                        if is_p2_mod_p4 ^ is_p3_mod_p4:
                           other_congs_count = total_congs - 3 # p2%p1, p3%p1, p_x%p4
                           if other_congs_count == 0:
                               return True

    # Cube case
    primes_with_cube = [p for p, e in factorization.items() if e >= 3]
    if len(primes_with_cube) == 1 and factorization[primes_with_cube[0]] == 3 and all(e == 1 for p_other, e in factorization.items() if p_other != primes_with_cube[0]):
        p = primes_with_cube[0]
        if 2 in primes: return False
        inv_graph = get_inverse_congruence_graph(primes)
        if not any(inv_graph.values()):
             p_rest = [pr for pr in primes if pr != p]
             if all((p*p-1)%pr!=0 and (p*p+p+1)%pr!=0 for pr in p_rest):
                 return True

    # Square cases
    primes_with_sq = [p for p, e in factorization.items() if e >= 2]
    if len(primes_with_sq) == 1 and factorization[primes_with_sq[0]] == 2 and all(e == 1 for p_other, e in factorization.items() if p_other != primes_with_sq[0]):
        p1 = primes_with_sq[0]
        p_rest = [p for p in primes if p != p1]

        # Case SQ1: g=p1^2 q..., q%p1^2==1 is only congruence
        total_congs = sum(len(v) for v in get_inverse_congruence_graph(primes).values())
        if total_congs == 1:
            for q_cand in p_rest:
                if q_cand % (p1*p1) == 1:
                    if q_cand % p1 == 1:
                        if not any((p1*p1 - 1) % pr == 0 for pr in p_rest):
                            return True

        # Case SQ3: g=2*p^2 (p odd)
        if p1 != 2 and len(primes) == 2 and 2 in primes:
            return True

        # Cases with no p%q==1 congruences
        if total_congs == 0:
            divisors_of_p1_plus_1 = [p for p in p_rest if (p1 + 1) % p == 0]
            # Case SQ4: |{q in p_rest | (p1+1)%q==0}| == 2
            if len(divisors_of_p1_plus_1) == 2:
                return True

            # Case SQ5: |{q | (p1+1)%q==0}|=1, |{r | (q-1)%r==0}|=1
            if len(divisors_of_p1_plus_1) == 1:
                q1 = divisors_of_p1_plus_1[0]
                p_rest_no_q1 = [p for p in p_rest if p != q1]
                divisors_of_q1_minus_1 = [p for p in p_rest_no_q1 if (q1 - 1) % p == 0]
                if len(divisors_of_q1_minus_1) == 1:
                    return True

    return False
