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

    # Helper to count congruences (p-1) % q == 0
    def get_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and (p1 - 1) % p2 == 0:
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
                # implies no congruence between p1 and p2. E.g., if p1 > p2 and (p1-1) % p2 == 0,
                # then p2 would be in graph[p1], meaning p2 = q1.
                if q1 != q2 and p1 != q2 and p2 != q1:
                    return True

        # Case 3 (sq-free, even, n=2pq)
        if factorization.get(2, 0) == 1 and len(primes) == 3:
            odd_primes = [p for p in primes if p != 2]
            p, q = min(odd_primes), max(odd_primes)
            if (q - 1) % p != 0:
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
                # If e1=2, Aut orders can have factors of p1-1, p1+1.
                # Check for p2 dividing p1^2-1.
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

            # Case 5: n = p1^2 * q * ..., exactly one congruence q = 1 (mod p1)
            if q_cong == p1:
                q_special = p_cong
                if (q_special - 1) % (p1 * p1) != 0:
                    # Condition: "none of the other prime factors of g ... divides p1 + 1"
                    if not any((p1 + 1) % p == 0 for p in p_rest):
                        return True
            # Case 6: n = p1^2 * ..., exactly one congruence p=1(mod q), where q!=p1
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

    # Helper to get congruence graph {p: {q | (p-1)%q==0}}
    def get_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and (p1 - 1) % p2 == 0:
                    graph[p1].add(p2)
        return graph

    # Helper to get inverse congruence graph {q: {p | (p-1)%q==0}}
    def get_inverse_congruence_graph(p_set):
        graph = collections.defaultdict(set)
        for p1 in p_set:
            for p2 in p_set:
                if p1 != p2 and (p1 - 1) % p2 == 0:
                    graph[p2].add(p1)
        return graph

    # Square-free cases (must be odd)
    if all(e == 1 for e in exponents):
        if 2 in primes: return False
        
        graph = get_congruence_graph(primes)
        inv_graph = get_inverse_congruence_graph(primes)
        total_congs = sum(len(v) for v in graph.values())

        # Case SF1: Two congruences, both of the form p=1(mod 3)
        if 3 in primes and len(inv_graph.get(3, [])) == 2 and total_congs == 2:
            return True

        # Case SF2/3: Three congruences with a specific graph structure
        if total_congs == 3:
            # Condition 1: largest prime with an outgoing edge cannot have out-degree >= 2
            p_cong_list = list(graph.keys())
            if not p_cong_list: return False # Should be unreachable if total_congs > 0

            if len(graph[max(p_cong_list)]) >= 2:
                return False

            # Condition 2: Structure must have a prime with in-degree 2 AND a prime with out-degree 2.
            # The original code missed the out-degree check, causing the bug for n=7917.
            # For n=7917={3,7,13,29}, graph is {7->3, 13->3, 29->7}.
            # in_degrees are {3:2, 7:1}, out_degrees are {7:1, 13:1, 29:1}.
            # It has in_degree 2 but no out_degree 2, so it should be false.
            has_indegree_2 = any(len(p_set) == 2 for p_set in inv_graph.values())
            has_outdegree_2 = any(len(p_set) == 2 for p_set in graph.values())
            
            if has_indegree_2 and has_outdegree_2:
                return True
        
        # If square-free and none of the above cases match, it's not gnu=5
        return False

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
                # For odd p1, prime factors of Aut orders are factors of p1-1, p1+1, p1^2+p1+1
                if not any((p1 * p1 - 1) % pr == 0 or (p1 * p1 + p1 + 1) % pr == 0 for pr in p_rest):
                    return True

    # Square cases
    primes_with_sq = [p for p, e in factorization.items() if e >= 2]
    if len(primes_with_sq) == 1 and factorization[primes_with_sq[0]] == 2 and all(e == 1 for p_other, e in factorization.items() if p_other != primes_with_sq[0]):
        p1 = primes_with_sq[0]
        p_rest = [p for p in primes if p != p1]
        
        graph = get_congruence_graph(primes)
        total_congs = sum(len(v) for v in graph.values())

        # Case SQ1/SQ2: Exactly one congruence p=1(mod q) exists.
        if total_congs == 1:
            p_cong = list(graph.keys())[0]
            q_cong = list(graph[p_cong])[0]
            
            # For gnu=5, this congruence must be of the form p=1(mod p1).
            if q_cong == p1:
                # Case SQ1: p is congruent to 1 mod p1^2.
                if (p_cong - 1) % (p1 * p1) == 0:
                    # Text: "nor does any such prime factor divide p_1^2 - 1"
                    # "such prime factor" refers to p_rest \ {p_cong}.
                    if not any((p1 * p1 - 1) % pr == 0 for pr in p_rest if pr != p_cong):
                        return True
                # Case SQ2: gnu=4 case is p=1(mod p1) but not mod p1^2, and (p1+1) not div by other primes.
                # If the last condition fails (i.e. p1+1 is divisible by another prime), it becomes gnu=5.
                else: 
                    if any((p1 + 1) % p == 0 for p in p_rest):
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
                # The condition is that q1-1 is divisible by one and only one SUCH factor,
                # where SUCH factors are the other primes in p_rest.
                divisors_of_q1_minus_1 = [p for p in p_rest_no_q1 if (q1 - 1) % p == 0]
                if len(divisors_of_q1_minus_1) == 1:
                    return True

    return False