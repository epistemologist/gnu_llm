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
        # counts[p] = k means p is congruent to 1 mod k other primes.
        counts = {p: len(targets) for p, targets in graph.items() if targets}

        if 2 not in primes:
            # Case 1 (sq-free, odd): Exactly one prime p is congruent to 1 mod q
            # for exactly two distinct q's. No other congruences.
            if len(counts) == 1 and list(counts.values())[0] == 2:
                return True

            # Case 2 (sq-free, odd): Exactly two primes p1, p2 are each congruent to 1 mod q
            # for exactly one q (q1, q2). q1 != q2. No other congruences exist.
            # And no congruence between p1 and p2 themselves.
            if len(counts) == 2 and sorted(list(counts.values())) == [1, 1]:
                p1, p2 = list(counts.keys())
                q1 = list(graph[p1])[0]
                q2 = list(graph[p2])[0]
                # The condition "the larger of these two is not thus congruent with respect to the smaller"
                # implies no congruence between p1 and p2. E.g., if p1 > p2 and p1 % p2 == 1,
                # then p2 would be in graph[p1], meaning p2 = q1.
                if q1 != q2 and p1 != q2 and p2 != q1:
                    return True

        # Case 3 (sq-free, even, n=2pq)
        if factorization.get(2, 0) == 1 and len(primes) == 3:
            odd_primes = [p for p in primes if p != 2]
            p, q = min(odd_primes), max(odd_primes)
            if q % p != 1:
                return True

        return False

    # Cases with square factors
    
    # Case 4: n = p1^2 * p2^2 * ..., where all groups are abelian.
    # Exponents must be 1 or 2, with at least two 2s.
    primes_with_exp2 = [p for p, e in factorization.items() if e == 2]
    if len(primes_with_exp2) >= 2 and all(e in [1, 2] for e in exponents):
        # Condition for all groups to be abelian:
        # For any two distinct primes p1, p2, p2 must not divide |Aut(G)| for any group G of order p1^e1.
        is_abelian = True
        for p1 in primes:
            for p2 in primes:
                if p1 == p2: continue
                e1 = factorization[p1]
                # If e1=1, Aut order is p1-1. We need (p1-1)%p2 != 0.
                if e1 == 1 and (p1 - 1) % p2 == 0:
                    is_abelian = False; break
                # If e1=2, Aut orders are p1(p1-1) and p1(p1-1)^2(p1+1).
                # Need p2 to not divide p1-1 or p1+1, i.e., (p1^2-1)%p2 != 0.
                if e1 == 2 and (p1 * p1 - 1) % p2 == 0:
                    is_abelian = False; break
            if not is_abelian: break
        if is_abelian:
            return True

    # Cases with exactly one square factor, p1^2
    primes_with_sq = [p for p, e in factorization.items() if e >= 2]
    if len(primes_with_sq) == 1 and factorization[primes_with_sq[0]] == 2 and all(e == 1 for p, e in factorization.items() if p not in primes_with_sq):
        p1 = primes_with_sq[0]
        p_rest = [p for p in primes if p != p1]

        graph = get_congruence_graph(primes)
        total_congs = sum(len(v) for v in graph.values())

        if total_congs == 1:
            p_cong = [p for p, targets in graph.items() if targets][0]
            q_cong = list(graph[p_cong])[0]

            # Case 5: n = p1^2 * q * ..., exactly one congruence q % p1 == 1
            if q_cong == p1:
                q_special = p_cong
                if q_special % (p1 * p1) != 1:
                    # Condition: "none of the other prime factors of g ... divides p1 + 1"
                    # "other" means primes in p_rest except q_special.
                    if not any((p1 + 1) % p == 0 for p in p_rest if p != q_special):
                        return True
            # Case 6: n = p1^2 * ..., exactly one congruence p%q=1, where q!=p1
            else: # q_cong != p1
                # Text says "none of these numbers divides p1^2-1". "these numbers" are p_rest = {p2, ..., p_lambda}.
                if not any((p1 * p1 - 1) % p == 0 for p in p_rest):
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

    # Helper to get congruence graph {p: {q | p%q==1}}
    def get_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and p1 % p2 == 1:
                    graph[p1].add(p2)
        return graph

    # Helper to get inverse congruence graph {q: {p | p%q==1}}
    def get_inverse_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and p1 % p2 == 1:
                    graph[p2].add(p1)
        return graph

    # Square-free cases (must be odd)
    if all(e == 1 for e in exponents):
        if 2 in primes: return False
        inv_graph = get_inverse_congruence_graph(primes)
        total_congs = sum(len(v) for v in inv_graph.values())

        # Case SF1: n=3qr..., exactly two primes are 1 mod 3, no other congruences
        if 3 in primes and len(inv_graph.get(3, [])) == 2 and total_congs == 2:
            return True

        # Case SF2/3: n=p1 p2 p3 p4..., p2%p1=1, p3%p1=1, and (p2%p4=1 XOR p3%p4=1), total 3 congs
        if len(primes) >= 4:
            graph = get_congruence_graph(primes)
            total_congs = sum(len(v) for v in graph.values())
            for p1 in primes:
                if len(inv_graph.get(p1, [])) == 2:
                    p2, p3 = tuple(inv_graph[p1])
                    p_rest = [p for p in primes if p not in [p1, p2, p3]]
                    for p4 in p_rest:
                        is_p2_mod_p4 = p2 % p4 == 1
                        is_p3_mod_p4 = p3 % p4 == 1
                        if is_p2_mod_p4 ^ is_p3_mod_p4:
                           # Check if these are the only three congruences
                           expected_congs = 3
                           current_congs = (p2 % p1 == 1) + (p3 % p1 == 1) + (is_p2_mod_p4 or is_p3_mod_p4)
                           if total_congs == expected_congs:
                               return True

    # Cube case: n = p1^3 * p2 * ...
    primes_with_cube = [p for p, e in factorization.items() if e >= 3]
    if len(primes_with_cube) == 1 and factorization[primes_with_cube[0]] == 3 and all(e == 1 for p_other, e in factorization.items() if p_other != primes_with_cube[0]):
        p1 = primes_with_cube[0]
        p_rest = [p for p in primes if p != p1]
        
        # Condition 1: "no one of these prime numbers is congruent to unity with respect to another one of them"
        graph = get_congruence_graph(primes)
        if not any(graph.values()):
            # Condition 2: "none of the orders of the groups of isomorphisms of the various groups of order p_1^3 
            # is divisible by one of the prime nunbers p_2,...,p_\lambda"
            if p1 == 2:
                # For p1=2, odd prime factors of Aut orders are 3, 7
                if 3 not in p_rest and 7 not in p_rest:
                    return True
            else:  # p1 is odd
                # For odd p1, prime factors of Aut orders are factors of p1^2-1 and p1^2+p1+1
                if not any((p1 * p1 - 1) % pr == 0 or (p1 * p1 + p1 + 1) % pr == 0 for pr in p_rest):
                    return True

    # Square cases
    primes_with_sq = [p for p, e in factorization.items() if e >= 2]
    if len(primes_with_sq) == 1 and factorization[primes_with_sq[0]] == 2 and all(e == 1 for p_other, e in factorization.items() if p_other != primes_with_sq[0]):
        p1 = primes_with_sq[0]
        p_rest = [p for p in primes if p != p1]
        
        graph = get_congruence_graph(primes)
        total_congs = sum(len(v) for v in graph.values())

        # Case SQ1: g=p1^2 q..., q%p1^2==1 is only congruence
        if total_congs == 1:
            p_cong = [p for p, targets in graph.items() if targets][0]
            q_cong_set = graph[p_cong]
            # There is only one congruence: p_cong % q_cong == 1.
            # Here the condition is that some prime is congruent to 1 mod p1^2
            # It must be that q_cong = p1, and p_cong % p1^2 == 1.
            # But p_cong % p1 == 1 is the congruence, not p_cong % p1^2.
            # The text says "congruent with respect to p1^2", which is a stronger condition.
            if p_cong % (p1*p1) == 1 and len(q_cong_set) == 1 and p1 in q_cong_set:
                # Text: "nor does any such prime factor divide p_1^2 - 1"
                if not any((p1 * p1 - 1) % pr == 0 for pr in p_rest if pr != p_cong):
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