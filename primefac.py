#!/usr/bin/python3

from math import sqrt, log2, ceil, gcd, inf
from itertools import takewhile, compress, count
from multiprocessing import Process, Queue as mpQueue
from random import randrange#, seed; seed(0)
from time import time
from datetime import datetime, timedelta

version = "2.0.12"

def primegen(limit=inf):
    """
    Generates primes < limit almost lazily by a segmented sieve of Eratosthenes.
    Memory usage depends on the sequence of prime gaps.  Unconditionally, we can
    say that memory usage is O(p^0.7625), where p is the most-recently-yielded
    prime.  On the Riemann Hypothesis and Cramer's conjecture, we can bring this
    down to O(p^0.75 * log(p)) and O(sqrt(p) * log(p)^2), respectively.
    Input: limit -- a number (default = inf)
    Output: sequence of integers
    Examples:
    >>> list(islice(primegen(), 20))
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    >>> list(primegen(73))
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    """
    # We don't sieve 2, so we ought to be able to get sigificant savings by halving the length of the sieve.
    # But the tiny extra computation involved in that seems to exceed the savings.
    yield from takewhile(lambda x: x < limit, (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47))
    pl, pg = [3,5,7], primegen()
    for p in pl: next(pg)
    n = next(pg); nn = n*n
    while True:
        n = next(pg)
        ll, nn = nn, n*n
        sl = (nn - ll)
        sieve = bytearray([True]) * sl
        for p in pl:
            k = (-ll) % p
            sieve[k::p] = bytearray([False]) * ((sl-k)//p + 1)
        if nn > limit: break                                            # TODO bring this condition up to the while statement
        yield from compress(range(ll,ll+sl,2), sieve[::2])
        pl.append(n)
    yield from takewhile(lambda x: x < limit, compress(range(ll,ll+sl,2), sieve[::2]))

def ilog(x, b):                                 # TODO: investigate optimization starting from x.bin_length() * 2 // b
    """
    Greatest integer l such that b**l <= x
    Input: x, b -- integers
    Output: An integer
    Examples:
    >>> ilog(263789, 10)
    5
    >>> ilog(1023, 2)
    9
    """
    l = 0
    while x >= b:
        x //= b
        l += 1
    return l
    # TODO possible optimization route: x.bit_length() == ilog(x, 2) + 1; we can therefore use x.bit_length() * 2 // b as a
    #      1st approximation to ilog(x, b), then compute pow(b, x.bit_length() * 2 // b), then compare that to x and adjust.

try: from math import isqrt
except ImportError:
    def isqrt(n):
        """
        Greatest integer less than or equal to the square root of n.
        Shamelessly stolen from https://codegolf.stackexchange.com/a/9088.
        Input: n -- a whole number
        Output: An integer
        Examples:
        >>> list(map(isqrt, range(25)))
        [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        """
        if n < 0: return int(n)
        c = n*4//3
        d = c.bit_length()
        a = d>>1
        if d&1:
            x = 1 << a
            y = (x + (n >> a)) >> 1
        else:
            x = (3 << a) >> 2
            y = (x + (c >> a)) >> 1
        if x != y:
            x, y = y, (y + n//y) >> 1
            while y < x: x, y = y, (y + n//y) >> 1
        return x

try: from math import prod as iterprod
except ImportError:
    def iterprod(l):
        """
        Product of the elements of any iterable.  Product of empty iterable == 1.
        Input: l -- iterable
        Output: A number
        Examples:
        >>> iterprod(range(1, 8))
        5040
        """
        z = 1
        for x in l: z *= x
        return z

try:
    pow(2, -1, 3)
    def modinv(a, m):
        """
        Returns the inverse of a modulo m, normalized to lie between 0 and m-1.  If
        a is not coprime to m, return None.
        Input:
            a -- an integer coprime to m
            n -- a positive integer
        Output: None or an integer x between 0 and m-1 such that (a * x) % m == 1
        Examples:
        >>> [modinv(1,1), modinv(2,5), modinv(5,8), modinv(37,100), modinv(12,30)]
        [0, 3, 5, 73, None]
        """
        return pow(a, -1, m)
except ValueError:
    def modinv(a, m):
        """
        Returns the inverse of a modulo m, normalized to lie between 0 and m-1.  If
        a is not coprime to m, return None.
        Input:
            a -- an integer coprime to m
            n -- a positive integer
        Output: None or an integer x between 0 and m-1 such that (a * x) % m == 1
        Examples:
        >>> [modinv(1,1), modinv(2,5), modinv(5,8), modinv(37,100), modinv(12,30)]
        [0, 3, 5, 73, None]
        """
        if m <= 0: return None
        M, x, r = m, 1, 0
        while m != 0:
            q = a//m
            a, m, r, x = m, a%m, x-q*r, r
        return x % M if a == 1 else None

def introot(n, r=2):    # TODO Newton iteration?
    """
    For returns the rth root of n, rounded to the nearest integer in the
    direction of zero.  Returns None if r is even and n is negative.
    Input:
        n -- an integer
        r -- a natural number or None
    Output: An integer
    Examples:
    >>> [introot(-729, 3), introot(-728, 3), introot(1023, 2), introot(1024, 2)]
    [-9, -8, 31, 32]
    """
    if n < 0: return None if r%2 == 0 else -introot(-n, r)
    if n < 2: return n
    if r == 1: return n
    if r == 2: return isqrt(n)
    #if r % 2 == 0: return introot(isqrt(n), r//2)      # TODO Check validity of this line.
    lower = upper = 1 << (n.bit_length() // r)
    while lower ** r >  n: lower >>= 2
    while upper ** r <= n: upper <<= 2
    while lower != upper - 1:
        mid = (lower + upper) // 2
        m = mid**r
        if   m == n: return  mid
        elif m <  n: lower = mid
        elif m >  n: upper = mid
    return lower

def ispower(n, r=0):
    """
    If r == 0:
        If n is a perfect power, return a tuple containing largest integer (in
        terms of magnitude) that, when squared/cubed/etc, yields n as the first
        component and the relevant power as the second component.
        If n is not a perfect power, return None.
    If r > 0:
        We check whether n is a perfect rth power; we return its rth root if it
        is and None if it isn't.
    Input:
        n -- an integer
        r -- an integer
    Output: An integer, a 2-tuple of integers, or None
    Examples:
    >>> [ispower(n) for n in [64, 25, -729, 1729]]
    [(8, 2), (5, 2), (-9, 3), None]
    >>> [ispower(64, r) for r in range(7)]
    [(8, 2), 64, 8, 4, None, None, 2]
    """
    #if r == 0: return any(ispower(n, r) for r in primegen(n.bit_length()+1))
    #return n == introot(n, r) ** r
    if r == 0:
        if n in (0, 1, -1): return (n, 1)
        for r in primegen(n.bit_length()+1):
            x = ispower(n, r)
            if x is not None: return (x, r)
        return None
    # TODO tricks for special cases
    if (r == 2) and (n & 2): return None
    if (r == 3) and (n & 7) in (2,4,6): return None
    x = introot(n, r)
    return None if x is None else (x if x**r == n else None)

def jacobi(a, n):
    """
    The Jacobi symbol (a|n).
    Input:
        a -- any integer
        n -- odd integer
    Output: -1, 0, or 1
    Examples:
    >>> [jacobi(a, 15) for a in [-10, -7, -4, -2, -1, 0, 1, 2, 4, 7, 10]]
    [0, 1, -1, -1, -1, 0, 1, 1, 1, -1, 0]
    >>> [jacobi(a, 13) for a in [-10, -9, -4, -2, -1, 0, 1, 2, 4, 9, 10]]
    [1, 1, 1, -1, 1, 0, 1, -1, 1, 1, 1]
    >>> [jacobi(a, 11) for a in [-10, -9, -4, -2, -1, 0, 1, 2, 4, 9, 10]]
    [1, -1, -1, 1, -1, 0, 1, -1, 1, 1, -1]
    """
    if (n%2 == 0) or (n < 0): return None # n must be a positive odd number     TODO delete this check?
    if (a == 0) or (a == 1): return a
    a, t = a%n, 1
    while a != 0:
        while not a & 1:
            a //= 2
            if n & 7 in (3, 5): t *= -1
        a, n = n, a
        if (a & 3 == 3) and (n & 3) == 3: t *= -1
        a %= n
    return t if n == 1 else 0

def isprime(n, tb=(3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59)): # TODO optimize the basis, possibly by varying it with n
    """
    BPSW primality test variant: we use the strong Lucas PRP test and preface
    the computation with trial division for speed.  No composites are known to
    pass the test, though it is suspected that infinitely many will do so.
    There are definitely no such errors below 2^64.
    This function is mainly a streamlined version of bpsw().
    Input:
        n -- integer.  Number to be examined.
        tb -- iterable of primes.  Basis for trial division.
    Output: True if probably prime; False if definitely composite.
    Examples:
    >>> [n for n in range(91) if isprime(1000*n+1)]
    [3, 4, 7, 9, 13, 16, 19, 21, 24, 28, 51, 54, 55, 61, 69, 70, 76, 81, 88, 90]
    >>> [isprime(rpn("38 ! 1 +")), isprime(rpn("38 ! 1 -"))]
    [False, True]
    """
    # 1.  Do some trial division with tb as the basis.
    if n % 2 == 0 or n < 3: return n == 2
    for p in tb:
        if n % p == 0: return n == p
    
    # 2.  If sprp(n,2) fails, return False.  If it succeeds, continue.
    t, s = (n - 1) // 2, 1
    while t % 2 == 0: t //= 2; s += 1
    #assert 1 + 2**s * t == n
    x = pow(2, t, n)
    if x != 1 and x != n - 1:
        for j in range(1, s):
            x = pow(x, 2, n)
            if x == 1: return False
            elif x == n - 1: break
        else: return False
    
    # 3.  Select parameters for slprp.
    for D in count(5, 4):
        j = jacobi(D, n)
        if j == 0: return D == n
        if j == -1: break
        D = -2 - D
        j = jacobi(D, n)
        if j == 0: return -D == n
        if j == -1: break
        if D == -13 and ispower(n,2): return False      # If n is square, then this loop amounts to very slow trial division.
    
    # Now run slprp(n, 1, (1 - D) // 4) and return the result.
    b = (1 - D) // 4
    if 1 < gcd(n, b) < n: return False
    s, t = 1, (n + 1) // 2
    while t % 2 == 0: s += 1; t //= 2
    v, w, q, Q = 2, 1, 1, 1
    for k in bin(t)[2:]:
        q = (q*Q) % n
        if k == '1': Q, v, w = (q*b) % n, (w*v - q) % n, (w*w - 2*q*b) % n
        else:        Q, w, v =  q       , (w*v - q) % n, (v*v - 2*q  ) % n
    # assert ( (2*w-v) * modinv(D,n) ) % n, v == lucasmod(t, 1, b, n)
    if v == 0 or ( (2*w-v) * modinv(D,n) ) % n == 0: return True
    q = pow(b, t, n)
    for _ in range(1, s):
        v = (v*v - 2*q) % n
        if v == 0: return True
        q = (q*q) % n
    return False

def pollardrho_brent(n, verbose=False):
    """
    Factors integers using Brent's variation of Pollard's rho algorithm.
    If n is prime, we immediately return n; if not, we keep chugging until a
    nontrivial factor is found.  This function calls the randomizer; two
    successive calls may therefore return two different results.
    Input: n -- number to factor
    Output: A factor of n.
    Examples:
    >>> n = rpn("20 ! 1 +"); f = pollardrho_brent(n); n % f
    0
    """
    if isprime(n): return n
    g = n
    while g == n:
        y, c, m, g, r, q = randrange(1, n), randrange(1, n), randrange(1, n), 1, 1, 1
        while g==1:
            x, k = y, 0
            for i in range(r): y = (y**2 + c) % n
            while k < r and g == 1:
                ys = y
                for i in range(min(m, r-k)):
                    y = (y**2 + c) % n
                    q = q * abs(x-y) % n
                g, k = gcd(q, n), k+m
            r *= 2
        if g==n:
            while True:
                ys = (ys**2+c)%n
                g = gcd(x-ys, n)
                if g > 1: break
    return g

def pollard_pm1(n, B1=100, B2=1000, verbose=False):       # TODO: What are the best default bounds and way to increment them?
    """
    Integer factoring function.  Uses Pollard's p-1 algorithm.  Note that this
    is only efficient if the number to be factored has a prime factor p such
    that p-1's largest prime factor is "small".  In this implementation, that
    tends to mean less than 10,000,000 or so.
    Input:
        n -- number to factor
        B1 -- Natural number.  Bound for phase 1.  Default == 100.
        B2 -- Natural number > B1.  Bound for phase 2.  Default == 1000.
    Output: A factor of n.
    Examples:
    >>> pollard_pm1(rpn("28 ! 1 - 239 //"))
    1224040923709997
    """
    if isprime(n): return n
    m = ispower(n)
    if m: return m[0]
    while True:
        pg = primegen()
        q = 2           # TODO: what about other initial values of q?
        p = next(pg)
        while p <= B1: q, p = pow(q, p**ilog(B1, p), n), next(pg)
        g = gcd(q-1, n)
        if 1 < g < n: return g
        while p <= B2: q, p = pow(q, p, n), next(pg)
        g = gcd(q-1, n)
        if 1 < g < n: return g
        # These bounds failed.  Increase and try again.
        B1 *= 10
        B2 *= 10

def mlucas(v, a, n):
    # Helper for williams_pp1().  Multiplies along a Lucas sequence mod n.
    v1, v2 = v, (v**2 - 2) % n
    for bit in bin(a)[3:]: v1, v2 = ((v1**2 - 2) % n, (v1*v2 - v) % n) if bit == "0" else ((v1*v2 - v) % n, (v2**2 - 2) % n)
    return v1
def williams_pp1(n, verbose=False):      # TODO: experiment with different values of v0, and implement the two-phase version
    """
    Integer factoring function.  Uses Williams' p+1 algorithm, single-stage
    variant.  Note that this is only efficient when the number to be factored
    has a prime factor p such that p+1's largest prime factor is "small".
    Input: n -- integer to factor
    Output: Integer.  A nontrivial factor of n.
    Example:
    >>> williams_pp1(315951348188966255352482641444979927)
    12403590655726899403
    """
    if isprime(n): return n
    m = ispower(n)
    if m: return m[0]
    for v in count(3):
        for p in primegen():
            e = ilog(isqrt(n), p)
            if e == 0: break
            for _ in range(e): v = mlucas(v, p, n)
            g = gcd(v - 2, n)
            if 1 < g < n: return g
            if g == n: break

def ecadd(p1, p2, p0, n):
    # Helper for ecm().  Adds two points on a Montgomery curve mod n.
    x1,z1 = p1; x2,z2 = p2; x0,z0 = p0
    t1, t2 = (x1-z1)*(x2+z2), (x1+z1)*(x2-z2)
    return (z0*pow(t1+t2,2,n) % n, x0*pow(t1-t2,2,n) % n)
def ecdub(p, A, n):
    # Helper for ecm().  Doubles a point on a Montgomery curve mod n.
    x, z = p; An, Ad = A
    t1, t2 = pow(x+z,2,n), pow(x-z,2,n)
    t = t1 - t2
    return (t1*t2*4*Ad % n, (4*Ad*t2 + t*An)*t % n)
def ecmul(m, p, A, n):
    # Helper for ecm().  Multiplies a point on a Montgomery curve mod n.
    if m == 0: return (0, 0)
    elif m == 1: return p
    else:
        q = ecdub(p, A, n)
        if m == 2: return q
        b = 1
        while b < m: b *= 2
        b //= 4
        r = p
        while b:
            if m&b: q, r = ecdub(q, A, n), ecadd(q, r, p, n)
            else:   q, r = ecadd(r, q, p, n), ecdub(r, A, n)
            b //= 2
        return r
def secm(n, B1, B2, seed):
    """
    Seeded ECM.  Helper function for ecm().  Returns a possibly-trivial divisor
    of n given two bounds and a seed.  Uses the two-phase algorithm on
    Montgomery curves.  See https://wp.me/prTJ7-zI and https://wp.me/prTJ7-A7
    for more details.  Most of the code for this function's "helpers" were
    shamelessly copied from the first of those links.
    Input:
        n -- Integer to factor
        B1 -- Integer.  Number of iterations for the first phase.
        B2 -- Integer.  Number of iterations for the second phase.
        seed -- Integer.  Selects the specific curve we'll be working on.
    Output: Integer.  A possibly-trivial factor of n.
    Examples:
    >>> secm(rpn('24 ! 1 -'), 100, 1000, 22)
    991459181683
    """
    u, v = (seed**2 - 5) % n, 4*seed % n
    p = pow(u, 3, n)
    Q, C = (pow(v-u,3,n)*(3*u+v) % n, 4*p*v % n), (p, pow(v,3,n))
    pg = primegen()
    p = next(pg)
    while p <= B1: Q, p = ecmul(p**ilog(B1, p), Q, C, n), next(pg)
    g = gcd(Q[1], n)
    if 1 < g < n: return g
    while p <= B2:
        # There is a trick that can speed up the second stage.  Instead of multiplying each prime by Q, we iterate over i from
        # B1 + 1 to B2, adding 2Q at each step; when i is prime, the current Q can be accumulated into the running solution.
        # Again, we defer the calculation of the GCD until the end of the iteration.    TODO: Implement and compare performance.
        Q = ecmul(p, Q, C, n)
        g *= Q[1]
        g %= n
        p = next(pg)
    return gcd(g, n)
def ecmparams(n):   # TODO: Better parameters.
    counter = 0
    for i in count():
        for _ in range(2*i+1):
            yield (2**i, 10 * 2**i, randrange(6,n), counter)
            counter += 1
        for j in range(i+1):
            yield (2**j, 10 * 2**j, randrange(6,n), counter)
            counter += 1
def ecm(n, paramseq=ecmparams, nprocs=1, verbose=False):
    """
    "Modern" integer factoring via elliptic curves.  Uses Montgomery curves, the
    two-phase algorithm, and (optionally) multiple processes.  The hard work is
    done by secm(); this function just does the managerial work of pulling a
    sequence of parameters out of a generator and feeding them into secm().
    Input:
        n -- number to factor
        paramseq -- sequence of parameters to feed into secm().  It must be an
                    infinite generator of 4-tuples (a,b,c,d), where a is the
                    number of iterations for the first phase, b is the number of
                    iterations for the second phase, c is a seed to select the
                    curve to work on, and d is an auxiliary used to count the
                    parameters generated so far.  We need a < b and 6 <= c < n.
        nprocs -- number of processes to use.  Default == 1.  Setting this > 1
                  is discouraged on "small" inputs because managing multiple
                  processes incurs significant overhead.
    Output:
        A factor of n.
        Note that if the parameter sequence calls the randomizer (which is
        currently the default behavior), then two successive calls may therefore
        return two different results.
    Examples:
    >>> n = 625793187653 * 991459181683 # = 620448401733239439359999 = 24! - 1
    >>> f = ecm(n)
    >>> (n//f) * f
    620448401733239439359999
    """
    g = n % 6
    if g % 2 == 0: return 2
    if g % 3 == 0: return 3
    if isprime(n): return n
    m = ispower(n)
    if m: return m[0]
    if nprocs == 1:
        for (B1,B2,seed,i) in paramseq(n):
            f = secm(n, B1, B2, seed)
            if 1 < f < n: return f
    assert nprocs != 1
    def factory(params, output): output.put(secm(*params))
    ps, facs, procs = paramseq(n), mpQueue(), []
    procs = [Process(target=factory, args=((n,)+next(ps)[:3], facs)) for _ in range(nprocs)]
    for p in procs: p.start()
    while True:
        g = facs.get()
        if 1 < g < n:
            for p in procs: p.terminate()
            return g
        for p in range(nprocs):                                                      # TODO: Try doing this with a process pool.
            if not procs[p].is_alive():
                del procs[p]
                break
        procs.append(Process(target=factory, args=((n,)+next(ps)[:3], facs)))
        procs[-1].start()

def sqrtmod_prime(a, p):
    """
    Solves x**2 == a (mod p) for x.  We assume that p is a prime and a is a
    quadratic residue modulo p.  If either of these assumptions is false, the
    return value is meaningless.
    The Cipolla-Lehmer section is my own.  The rest appears to be derived from
    https://codegolf.stackexchange.com/a/9088.
    Input:
        a -- natural number
        p -- prime number
    Output: whole number less than p
    Examples:
    >>> sqrtmod_prime(4, 5)
    3
    >>> sqrtmod_prime(13, 23)
    6
    >>> sqrtmod_prime(997, 7304723089)
    761044645
    """
    a %= p
    if p%4 == 3: return pow(a, (p+1) >> 2, p)
    elif p%8 == 5:
        v = pow(a << 1, (p-5) >> 3, p)
        return (a*v*(((a*v*v<<1)%p)-1))%p
    elif p%8 == 1:
        # CranPom ex 2.31, pg 112 / 126.  Pretty sure this amounts to Cipolla-Lehmer.
        if a == 0: return 0     # Necessary to avoid an infinite loop in the legendre section
        h = 2
        # while legendre(h*h - 4*a, p) != -1:
        while pow(h*h - 4*a, (p-1) >> 1, p) != p - 1: h += 1                            # TODO compare speed vs random selection
        #return ( lucasmod((p+1)//2, h, a, p)[1] * modinv(2, p) ) % p
        k, v, w, q, Q = (p+1)//2, 2, h % p, 1, 1
        for kj in bin(k)[2:]:
            q = (q*Q) % p
            if kj == '1': Q, v, w = (q*a) % p, (w*v - h*q) % p, (w*w - 2*q*a) % p
            else:         Q, w, v =  q       , (w*v - h*q) % p, (v*v - 2*q  ) % p
        return (v*k) % p
    else: return a # p == 2

def siqs(n, verbose=False):
    """
    Uses the Self-Initializing Quadratic Sieve to extract a factor of n.
    We make some calls to randrange, so the output may change from call to call.
    This is derived from https://github.com/skollmann/PyFactorise.
    Input:
        n -- number to factor
        verbose -- if True, print progress reports.  Default == False.
    Output:
        If n is prime, we return n.
        If n is composte, we return a nontrivial factor of n.
        If n is too small, we raise a ValueError.
    Examples:
    >>> siqs(factorial(24) - 1) in (625793187653, 991459181683)
    True
    """
    
    if (not isinstance(n, int)) or n < 0: raise ValueError("Number must be a positive integer.")
    
    if n < 2**64:
        if verbose: print("Number is too small for SIQS.  Using Pollard Rho instead.")
        return pollardrho_brent(n)
    
    if verbose: print("Factorizing %d (%d digits)..." % (n, len(str(n))))
    
    if isprime(n): return n
    
    if verbose: print("Number is composite.")
    if verbose: print("Checking whether it is a perfect power...")
    perfect_power = ispower(n)
    if perfect_power:
        if verbose: print(n, "=", "%d^%d." % (perfect_power[0], perfect_power[1]))
        return perfect_power[0]
    
    if verbose: print("Not a perfect power.")
    if verbose: print("Using SIQS on %d (%d digits)..." % (n, len(str(n))))
    
    if verbose: starttime, sievetime, latime = time(), 0, 0
    
    # Choose parameters nf (sieve of factor base) and m (for sieving in [-m,m].
    # Using similar parameters as msieve-1.52
    dig = len(str(n))
    if   dig <= 34: nf, m = 200, 65536
    elif dig <= 36: nf, m = 300, 65536
    elif dig <= 38: nf, m = 400, 65536
    elif dig <= 40: nf, m = 500, 65536
    elif dig <= 42: nf, m = 600, 65536
    elif dig <= 44: nf, m = 700, 65536
    elif dig <= 48: nf, m = 1000, 65536
    elif dig <= 52: nf, m = 1200, 65536
    elif dig <= 56: nf, m = 2000, 65536 * 3
    elif dig <= 60: nf, m = 4000, 65536 * 3
    elif dig <= 66: nf, m = 6000, 65536 * 3
    elif dig <= 74: nf, m = 10000, 65536 * 3
    elif dig <= 80: nf, m = 30000, 65536 * 3
    elif dig <= 88: nf, m = 50000, 65536 * 3
    elif dig <= 94: nf, m = 60000, 65536 * 9
    else:           nf, m = 100000, 65536 * 9
    
    class FactorBasePrime:
        """A factor base prime for the SIQS"""
        def __init__(self, p, tmem, lp):
            self.p, self.soln1, self.soln2, self.tmem, self.lp, self.ainv = p, None, None, tmem, lp, None
    
    # Compute and return nf factor base primes suitable for a SIQS on n.
    factor_base = []
    for p in primegen():
        if pow(n, (p-1) >> 1, p) == 1:          # if n is a quadratic residue mod p
            t = sqrtmod_prime(n, p)
            lp = round(log2(p))                 # This gets rounded for the sake of speed.
            factor_base.append(FactorBasePrime(p, t, lp))
            if len(factor_base) >= nf: break
    
    if verbose: print("Factor base size: %d primes.\nLargest prime in base: %d." % (nf, factor_base[-1].p))
    
    npolys, i_poly, prev_cnt, relcount, smooth_relations, required_relations_ratio = 0, 0, 0, 0, [], 1.05
    
    p_min_i, p_max_i = None, None
    for (i,fb) in enumerate(factor_base):
        if p_min_i is None and fb.p >= 400:                                            # 400 is a tunable parameter.
            p_min_i = i
        if p_max_i is None and fb.p > 4000:                                           # 4000 is a tunable parameter.
            p_max_i = i - 1
            break
    
    # The following may happen if the factor base is small, make sure that we have enough primes.
    if p_max_i is None: p_max_i = len(factor_base) - 1
    if p_min_i is None or p_max_i - p_min_i < 20: p_min_i = min(p_min_i, 5)    # TODO This line is problematic for some small n.
    
    target = sqrt(2 * float(n)) / m
    target1 = target / sqrt((factor_base[p_min_i].p + factor_base[p_max_i].p) / 2)
    
    while True:
    
        if verbose: temptime = time()
    
        if verbose: print("*** Phase 1: Finding smooth relations ***")
        required_relations = round(nf * required_relations_ratio)
        if verbose: print("Target: %d relations" % required_relations)
        enough_relations = False
        while not enough_relations:
            
            if i_poly == 0: # Compute the first of a set of polynomials
                # find q such that the product of factor_base[q_i] is about sqrt(2n)/m; try a few sets to find a good one
                best_q, best_a, best_ratio = None, None, None
                for _ in range(30):
                    a, q = 1, []
                    
                    while a < target1:
                        p_i = 0
                        while p_i == 0 or p_i in q: p_i = randrange(p_min_i, p_max_i+1)
                        a *= factor_base[p_i].p
                        q.append(p_i)
                    
                    ratio = a / target
                    
                    # ratio too small seems to be not good
                    if (best_ratio is None or (ratio >= 0.9 and ratio < best_ratio) or best_ratio < 0.9 and ratio > best_ratio):
                        best_q, best_a, best_ratio = q, a, ratio
                a, q = best_a, best_q
                
                s, B = len(q), []
                for l in range(s):
                    fb_l = factor_base[q[l]]
                    q_l = fb_l.p
                    assert a % q_l == 0
                    gamma = (fb_l.tmem * modinv(a // q_l, q_l)) % q_l
                    if gamma > q_l // 2: gamma = q_l - gamma
                    B.append(a // q_l * gamma)
                
                b = sum(B) % a
                b_orig = b
                if (2 * b > a): b = a - b
                
                assert 0 < b and 2 * b <= a and (b * b - n) % a == 0
                
                g1, g2, g3, ga, gb = b * b - n, 2 * a * b, a * a, a, b_orig
                ha, hb = a, b
                for fb in factor_base:
                    if a % fb.p != 0:
                        fb.ainv = modinv(a, fb.p)
                        fb.soln1 = (fb.ainv * (fb.tmem - b)) % fb.p
                        fb.soln2 = (fb.ainv * (-fb.tmem - b)) % fb.p
            
            else: # Compute the (i+1)-th polynomial, given that g is the i-th polynomial.
                #v = lowest_set_bit(i) + 1
                #z = -1 if ceil(i / (2 ** v)) % 2 == 1 else 1
                z = -1 if ceil(i_poly / (1 + (i_poly ^ (i_poly-1)))) % 2 == 1 else 1
                #b = (g.b + 2 * z * B[v - 1]) % g.a
                b = (gb + 2 * z * B[(i_poly & (-i_poly)).bit_length() - 1]) % ga
                a = ga
                b_orig = b
                if (2 * b > a): b = a - b
                assert (b * b - n) % a == 0
                
                g1, g2, g3, ga, gb = b * b - n, 2 * a * b, a * a, a, b_orig
                ha, hb = a, b
                for fb in factor_base:
                    if a % fb.p != 0:
                        fb.soln1 = (fb.ainv * ( fb.tmem - b)) % fb.p
                        fb.soln2 = (fb.ainv * (-fb.tmem - b)) % fb.p
            
            
            i_poly += 1
            npolys += 1
            if i_poly >= 2 ** (len(B) - 1): i_poly = 0
            # BEGIN SIEVING.  Most of our time is spent between here and the "END SIEVING" comment.
            sieve_array = [0] * (2 * m + 1)
            for fb in factor_base:
                if fb.soln1 is None: continue
                p, lp = fb.p, fb.lp
                a_start_1 = fb.soln1 - ((m + fb.soln1) // p) * p
                if p > 20:
                    for a in range(a_start_1 + m, 2 * m + 1, p): sieve_array[a] += lp
                    a_start_2 = fb.soln2 - ((m + fb.soln2) // p) * p
                    for a in range(a_start_2 + m, 2 * m + 1, p): sieve_array[a] += lp
            # Perform the trial division step of the SIQS.
            limit = round(log2(m * sqrt(float(n)))) - 25                # 25 is a tunable parameter.  The rounding is for speed.
            for (i,sa) in enumerate(sieve_array):
                if sa >= limit:
                    x = i - m
                    gx = (g3 * x + g2) * x + g1
                    # Determine whether gx can be fully factorized into primes from the factor base.
                    # If so, store the indices of the factors from the factor base as divisors_idx. If not, set that to None.
                    a, divisors_idx = gx, []
                    for (fbi,fb) in enumerate(factor_base):
                        if a % fb.p == 0:
                            exp = 0
                            while a % fb.p == 0:
                                a //= fb.p
                                exp += 1
                            divisors_idx.append((fbi, exp))
                        if a == 1:
                            u = ha * x + hb
                            v = gx
                            assert (u * u) % n == v % n
                            smooth_relations.append((u, v, divisors_idx))
                            break
                    relcount = len(smooth_relations)
                    if relcount >= required_relations:
                        enough_relations = True
                        break
            # END SIEVING.  Most of our time is spent between here and the "BEGIN SIEVING" comment.
            
            if verbose and relcount > 0 and (relcount >= required_relations or i_poly % 8 == 0 or relcount > prev_cnt):
                frac = relcount / required_relations
                t = time() - starttime              # Time spent so far
                ett = t / frac                      # Estimated total time
                ettc = ett - t                      # Estimated time to completion
                eta = datetime.isoformat(datetime.now() + timedelta(seconds=ettc), sep=' ', timespec='seconds')
                print('\b'*256 + "Found %d rels (%.1f%%) on %d polys.  ETA %ds (%s). " % (relcount, 100*frac, npolys, ettc, eta), end='', flush=True)
                prev_cnt = relcount
        
        if verbose:
            sievetime += time() - temptime
            print("\nSieving time: %f seconds." % sievetime)
            temptime = time()
            
            print("*** Phase 2: Linear Algebra ***")
            print("Building matrix for linear algebra step...")
        M = [0] * nf
        mask = 1
        for sr in smooth_relations:
            for (j,exp) in sr[2]:
                if exp % 2: M[j] += mask
            mask <<= 1
        
        if verbose: print("Finding perfect squares using matrix and factors from perfect squares...")
        # Perform the linear algebra step of the SIQS: do fast Gaussian elimination to determine pairs of perfect squares mod n.
        # Use the optimisations described in [Çetin K. Koç and Sarath N. Arachchige. 'A Fast Algorithm for Gaussian Elimination
        # over GF(2) and its Implementation on the GAPP.' Journal of Parallel and Distributed Computing 13.1 (1991): 118-122].
        if verbose: gaussstart = time()
        row_is_marked = bytearray([False]) * relcount
        pivots = [-1] * nf
        for j in range(nf):
            M_j = M[j]
            i = (M_j & (-M_j)).bit_length() - 1         #i = -1 if M[j] == 0 else lowest_set_bit(M[j])
            # i is now the row of the first nonzero entry in column j, or -1 if no such row exists.
            if i > -1:
                pivots[j] = i
                row_is_marked[i] = True
                for k in range(nf):
                    if (M[k] >> i) & 1 and k != j:  # test M[i][k] == 1
                        M[k] ^= M_j         # add column j to column k mod 2
            if verbose and ((nf - j) % 100 == 0 or nf == j + 1):
                frac = (j+1) / nf
                t = time() - gaussstart             # Time spent so far
                ett = t / frac                      # Estimated total time
                ettc = ett - t                      # Estimated time to completion
                eta = datetime.isoformat(datetime.now() + timedelta(seconds=ettc), sep=' ', timespec='seconds')
                print('\b'*256 + "Gaussian elimination: %d/%d (%.1f%%).  ETA %ds (%s). " % (j+1, nf, 100*frac, ettc, eta), end='', flush=True)
        if verbose: print("\nGaussian elimination time: %f seconds" % (time() - gaussstart))
        attempts = 0
        for i in range(relcount):
            if not row_is_marked[i]:
                square_indices = [i]
                for j in range(nf):
                    if (M[j] >> i) & 1:  # test M[i][j] == 1
                        square_indices.append(pivots[j])
                # Given the solution encoded by square_indices, try to find a factor of n, and return it.
                attempts += 1
                if verbose: print('\b'*42 + "Attempt #%d" % attempts, end='', flush=True)
                # Given on of the solutions returned by siqs_solve_matrix and the corresponding smooth relations,
                # calculate the pair (sqrt1, sqrt2) such that sqrt1^2 = sqrt2^2 (mod n).
                sqrt1, sqrt2 = 1, 1
                for idx in square_indices:
                    sqrt1 *= smooth_relations[idx][0]
                    sqrt2 *= smooth_relations[idx][1]
                sqrt2 = isqrt(sqrt2)
                assert (sqrt1 * sqrt1) % n == (sqrt2 * sqrt2) % n
                factor = gcd(sqrt1 - sqrt2, n)
                if 1 != factor != n:
                    if verbose:
                        print(" succeeded.")
                        latime += time() - temptime
                        print("Linear algebra time: %f seconds." % latime)
                        totaltime = time() - starttime
                        print("Total time: %f seconds." % (totaltime))
                    return factor
        
        if verbose:
            print('\b'*42 + "All %d attempts failed." % attempts)
            latime += time() - temptime
            print("Linear algebra time: %f seconds." % latime)
            
            print("Failed to find a solution. Finding more relations...")
        required_relations_ratio += 0.05

def multifactor(n, methods=(pollardrho_brent, pollard_pm1, williams_pp1, ecm, siqs), verbose=False):
    """
    Integer factoring function.  Uses several methods in parallel.  Waits for a
    function to return, kills the rest, and reports.  Note that two successive
    calls may return different results depending on which method finishes first
    and whether any methods call the randomizer.
    Input:
        n -- number to factor
        methods -- list of functions to run.
    Output: A factor of n.
    Examples:
    >>> methods = [pollardrho_brent, pollard_pm1, williams_pp1, ecm, siqs]
    >>> n = rpn("24 ! 1 -"); f = multifactor(n, methods); n%f
    0
    """
    # Note that the multiprocing incurs relatively significant overhead.  Only call this if n is proving difficult to factor.
    def factory(method, n, verbose, output): output.put((method(n, verbose=verbose), str(method).split()[1]))
    factors = mpQueue()
    procs = [Process(target=factory, args=(m, n, verbose, factors)) for m in methods]
    for p in procs: p.start()
    (f, g) = factors.get()
    for p in procs: p.terminate()
    names = {"pollardrho_brent":"prb", "pollard_pm1":"p-1", "williams_pp1":"p+1"}
    return (f, names.get(g, g))

def primefac(n, trial=1000, rho=42000, verbose=False, methods=(pollardrho_brent,)):
    """
    Generates the prime factors of the input.  Factors that appear x times are
    yielded x times.
    Input:
        n -- the number to be factored
        trial -- Trial division is performed up to this limit and no further.
                 We trial divide by 2 and 3 whether the user wants to or not.
                 Default == 1000.
        rho -- How long to have Pollard's rho chug before declaring a cofactor
               difficult.  Default == 42,000 iterations.  Floating-point inf is
               an acceptable value.
        methods -- Use these methods on difficult cofactors.  If the tuple has
                   more than 1 element, we have multifactor handle it.  Calling
                   multifactor has a high overhead, so when the tuple has a
                   single element, we call that function directly.  The default
                   is (pollardrho_brent,).  Each function f in methods must
                   accept a single number n as its argument.  If n is prime,
                   f(n) must return n. If n is composite, f(n) must return a
                   number strictly between 1 and n that evenly divides n.
                   Giving up is not allowed.
    Output: Prime factors of n
    Examples:
    >>> list(primefac(1729))
    [7, 13, 19]
    >>> list(sorted(primefac(rpn("24 ! 1 -"))))
    [625793187653, 991459181683]
    """
    # Obtains a complete factorization of n, yielding the prime factors as they are obtained.
    # If the user explicitly specifies a splitting method, use that method.  Otherwise,
    # 1.  Pull out small factors with trial division.
    # 2.  Do a few rounds of Pollard's Rho algorithm.
    # 3.  Launch multifactor on the remainder.  Note that multifactor's multiprocessing incurs relatively significant overhead.
    
    if verbose: print("\nFactoring %d (%d digits):" % (n, len(str(n))))
    if n < 0:
        if verbose: print("Trial division: -1")
        yield -1; n *= -1
    if n < 2: return
    if isprime(n):
        if verbose: print("Number is prime.")
        yield n
        return
    if verbose: print("Number is composite.")
    
    # Trial division
    factors, nroot = [], isqrt(n)
    for p in primegen(max(4,trial)+1): # Note that we check for 2 and 3 whether the user wants to or not.
        if verbose: print('\b'*80 + "Trial division:", p, end='', flush=False)
        if n%p == 0:
            while n%p == 0:
                if verbose: print('\b'*80 + "Trial division:", p)
                yield p
                n //= p
            nroot = isqrt(n)
        if p > nroot:
            if n != 1:
                if verbose:
                    laststr = "Trial division: " + str(p)
                    print('\b'*80 + str(n) + " " * len(laststr))
                yield n
            return
    if isprime(n):
        if verbose:
            laststr = "Trial division: " + str(p)
            print('\b'*80 + str(n) + " " * len(laststr))
        yield n
        return
    
    if verbose: print("\nTrial division finished.")
    
    # TODO: Fermat's method?
    
    # Pollard's rho
    factors, difficult = [n], []
    while len(factors) != 0:
        rhocount = 0
        n = factors.pop()
        if verbose: print("Beginning Pollard Rho on %d." % n)
        try:
            g = n
            while g == n:
                x, c, g = randrange(1, n), randrange(1, n), 1
                y = x
                while g==1:
                    if rhocount >= rho: raise Exception
                    rhocount += 1
                    x = (x**2 + c) % n
                    y = (y**2 + c) % n
                    y = (y**2 + c) % n
                    g = gcd(x-y, n)
            # We now have a nontrivial factor g of n.  If we took too long to get here, we're actually at the except statement.
            f = n // g
            if verbose: print("Pollard Rho split %d into %d and %d." % (n, g, f))
            if isprime(g):
                if verbose: print(g, "is prime.")
                yield g
            else:
                if verbose: print(g, "is composite.")
                factors.append(g)
            if isprime(f):
                if verbose: print(f, "is prime.")
                yield f
            else:
                if verbose: print(f, "is composite.")
                factors.append(f)
        except Exception:
            if verbose: print("Declaring", n, "difficult.")
            difficult.append(n) # Factoring n took too long.  We'll have multifactor chug on it.
    
    # TODO: P-1 by itself?
    # TODO: ECM by itself?
    
    names = {"pollardrho_brent":"prb", "pollard_pm1":"p-1", "williams_pp1":"p+1", "ecm":"ecm", "siqs":"siqs"}
    factors = difficult
    if len(methods) == 1:
        method = methods[0]
        name = names[str(method).split()[1]]
        while len(factors) != 0:
            n = factors.pop()
            if verbose: print("Beginning %s on %d." % (name, n))
            f = method(n, verbose=verbose)
            g = n // f
            if verbose: print("%s split %d into %d and %d." % (name, n, f, g))
            if isprime(f):
                if verbose: print(f, "is prime.")
                yield f
            else:
                if verbose: print(f, "is composite.")
                factors.append(f)
            if isprime(g):
                if verbose: print(g, "is prime.")
                yield g
            else:
                if verbose: print(g, "is composite.")
                factors.append(g)
    else:
        while len(factors) != 0:
            n = factors.pop()
            if verbose:
                #methodstring = ", ".join(names[str(m).split()[1]] for m in methods)
                if len(methods) == 2: methodstring = " and ".join(names[str(m).split()[1]] for m in methods)
                else:
                    assert len(methods) > 2
                    methodstring = ", ".join(names[str(m).split()[1]] for m in methods[:-1])
                    methodstring += ", and " + names[str(methods[-1]).split()[1]]
                if verbose: print("Beginning %s on %d." % (methodstring, n))
            f, name = multifactor(n, methods=methods, verbose=verbose)
            g = n // f
            if verbose: print("%s split %d into %d and %d." % (name, n, f, g))
            if isprime(f):
                if verbose: print(f, "is prime.")
                yield f
            else:
                if verbose: print(f, "is composite.")
                factors.append(f)
            if isprime(g):
                if verbose: print(g, "is prime.")
                yield g
            else:
                if verbose: print(g, "is composite.")
                factors.append(g)

def rpn(instr):
    stack = []
    for token in instr.split():
        if set(token).issubset("1234567890"): stack.append(int(token))
        elif len(token) > 1 and token[0] == '-' and set(token[1:]).issubset("1234567890"): stack.append(int(token))
        elif token in ('+', '-', '*', '/', '//', '%', '**', 'x', 'xx'):   # binary operators
            b = stack.pop()
            a = stack.pop()
            if   token == '+' : res = a  + b
            elif token == '-' : res = a  - b
            elif token == '*' : res = a  * b
            elif token == 'x' : res = a  * b
            elif token == '/' : res = a // b
            elif token == '//': res = a // b
            elif token == '%' : res = a  % b
            elif token == '**': res = a ** b
            elif token == 'xx': res = a ** b
            stack.append(res)
        elif token in ('!', 'f', '#', 'p'):                             # unary operators
            a = stack.pop()
            if   token == '!' : res = iterprod(range(1, a+1))
            elif token == 'f' : res = iterprod(range(1, a+1))
            elif token == '#' : res = iterprod(primegen(a+1))
            elif token == 'p' : res = iterprod(primegen(a+1))
            stack.append(res)
        else: raise Exception
    for x in stack: yield int(x)

usage = """
This is primefac version 2.0.11.

USAGE:
    primefac [-vs|-sv] [-v|--verbose] [-s|--summary|--summarize] [-t=NUMBER]
             [-r=NUMBER] [-m=[prb][,p-1][,p+1][,ecm][,siqs]] rpn
    
    The value "rpn" is an expression in reverse Polish notation, also called
    postfix notation, and is evaluated using integer arithmetic.  The values
    that remain on the stack after evaluation are then factored in sequence.
    
    "-t" specifies the largest prime to use for trial division.  The default
    value for this parameter is 1000.  Using "-t=inf" will make primefac use
    trial division exclusively.
    
    "-r" is the number of iterations of Pollard's rho algorithm to do before
    calling a cofactor "difficult".  Default == 42,000.  Use "-r=inf" to use
    Pollard's rho algorithm exclusively once trial division is completed.
    
    If verbosity is invoked, then we provide progress reports and also state
    which algorithms produced which factors during the multifactor phase.
    
    If the summary and verbosity flags are absent, then the output should be
    identical to the output of the GNU factor command, modulo permutation of
    the factors.  If the verbosity flag is invoked, then we provide progress
    reports, turn on the summary flag, and state which methods yielded which
    factors during the multifactor phase.
    
    If the summary flag is present, then the output is modified by including
    a single newline between each item's output, before the first, and after
    the last.  Each item's output is also modified by printing a second line
    of data summarizing the results by indicating the number of digits (base
    10) in the input, the number of digits (base 10) in each factor, and the
    factors' multiplicities.  For example:
    
        >>> user@computer:~$ primefac  -s   24 ! 1 -   7 f
        >>> 
        >>> 620448401733239439359999: 991459181683 625793187653
        >>> Z24  =  P12 x P12  =  625793187653 x 991459181683
        >>> 
        >>> 5040: 2 2 2 2 3 3 5 7
        >>> Z4  =  P1^4 x P1^2 x P1 x P1  =  2^4 x 3^2 x 5 x 7
        >>> 
        >>> user@computer:~$
    
    Note that primes in the ordinary output lines are listed in the order in
    which they were found, while primes in the summary lines are reported in
    strictly-increasing order.
    
    The -v and -s flags may be collapsed into a single flag, -vs or -sv, but
    recall that -v implies -s.
    
    The "-m" flag controls what methods are used during the difficult-factor
    phase.  The "prb" and "ecm" options can be provided several times to use
    multiple concurrent instances of these methods.  Concurrent applications
    of the p-1, p+1, or SIQS methods confers no benefit over a single usage,
    so repeated listings of those methods are ignored.
    
    This program can be imported into Python scripts as a module; however, I
    recommend importing them from the module "labmath" instead.  This module
    is available via pip (https://pypi.org/project/labmath/).


INTERNAL DETAILS:
    Factoring:
        We use a three-stage factoring algorithm.
        1.  Trial divide with all primes less than or equal to the specified
            limit.  We trial divide by 2 and 3 regardless of the limit.
        2.  Run Pollard's rho algorithm on whatever remains.  This algorithm
            may split a number into two composite cofactors.  Such cofactors
            remain here until they survive the specified number of rounds of
            the rho algorithm.
        3.  Subject each remaining cofactor to a set of up to five factoring
            methods in parallel:
                Pollard's rho algorithm with Brent's improvement,
                Pollard's p-1 method,
                Williams' p+1 method,
                the elliptic curve method,
                and the self-initializing quadratic sieve.
    
    RPN:
        The available binary operators are +, -, *, //, %, and **, which all
        indicate the same operations here as they indicate in Python3 source
        code; i.e., they denote addition, subtraction, multiplication, floor
        division, remaindering, and powering.  The available unary operators
        are ! and #, which denote the factorial and primorial, respectively.
        For terminal syntax compatibility reasons, the RPN expression may be
        enclosed in quotes, and five aliases are allowed: x for *, / for //,
        xx for **, f for !, and p for #.


CREDITS:
    Significant parts of this code are derived or outright copied from other
    people's work.  In particular, the SIQS code was derived mostly verbatim
    from https://github.com/skollmann/PyFactorise by Stephan Kollmann, while
    the functions to manipulate points on elliptic curves were copied from a
    reply to the blog post at http://programmingpraxis.com/2010/04/23/.  The
    rest, I believe, is my own work, but I may have forgotten something.

"""

# TODO graceful exit on keyboard interrupt
# TODO timeout?
# TODO niceness in multifactor?
# TODO can mlucas be replaced with lucasu_mod/lucasv_mod from gmpy2?
if __name__ == "__main__":
    from sys import exit, argv
    if len(argv) == 1: exit(usage)
    if any("h" in arg for arg in argv[1:]) or any("?" in arg for arg in argv[1:]): exit(usage)
    start, rpx, tr, rr, verbose, summarize = 1, [], 1000, 42000, False, False
    ms = {"prb":pollardrho_brent, "p-1":pollard_pm1, "p+1":williams_pp1, "ecm":ecm, "siqs":siqs}
    methods = (pollardrho_brent, pollard_pm1, williams_pp1, ecm, siqs)
    try:
        for arg in argv[1:]:
            if arg in ("-v", "--verbose"): verbose, summarize = True, True
            elif arg in ("-s", "--summary", "--summarize"): summarize = True
            elif arg in ("-vs", "-sv"): verbose, summarize = True, True
            elif arg[:3] == "-t=": tr = inf if arg[3:] == "inf" else int(arg[3:])    # Maximum number for trial division
            elif arg[:3] == "-r=": rr = inf if arg[3:] == "inf" else int(arg[3:])    # Number of rho rounds before multifactor
            elif arg[:3] == "-m=": #methods = tuple(ms[x] for x in arg[3:].split(',') if x in ms)
                methods = []
                for x in arg[3:].split(','):
                    if x in ms:
                        if x in ("p-1", "p+1", "siqs") and ms[x] in methods: continue
                        methods.append(ms[x])
            else: rpx.append(arg)
        nums = rpn(' '.join(rpx))
    except: exit("Error while parsing arguments.  To view usage, invoke with -h, --help, or -?.")
    if summarize: print()
    for n in nums:
        assert isinstance(n, int)
        print("%d:" % n, end='', flush=True)
        f = {}
        flist = []
        for p in primefac(n, trial=(n if tr == "inf" else tr), rho=rr, verbose=verbose, methods=methods):
            f[p] = f.get(p, 0) + 1
            flist.append(p)
            if not verbose: print(" %d" % p, end='', flush=True)
            assert isprime(p) and n%p == 0, (n, p)
        if verbose:
            print()
            print("%d:" % n, end='')
            for p in flist: print(" %d" % p, end='')
        print()
        if summarize:
            print("Z%d  = " % len(str(n)), end=' ')
            outstr = ""
            for p in sorted(f):
                if f[p] == 1: outstr += "P%d x " % len(str(p))
                else: outstr += "P%d^%d x " % (len(str(p)), f[p])
            outstr = outstr[:-2] + " = "
            for p in sorted(f):
                if f[p] == 1: outstr += " %d x" % p
                else: outstr += " %d^%d x" % (p, f[p])
            print(outstr[:-2])
            print()
    exit()

# testdata = [factorial(24) - 1, factorial(38) + 1, factorial(40) - 1, factorial(44) + 1, factorial(54) + 1]

