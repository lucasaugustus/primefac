primefac version 2.0.12
=======================

This is a module and command-line utility for factoring integers.  As a module, we provide a primality test, several functions for extracting a non-trivial factor of an integer, a generator that yields all of a number's prime factors (with multiplicity), and ancillary functions used in calculating these things.  As a command-line utility, this project aims to replace GNU's ``factor`` command with a more versatile utility --- in particular, this utility can operate on arbitrarily large numbers, uses multiple cores in parallel, uses better algorithms, handles input in reverse Polish notation, and can be tweaked via command-line flags.  Specifically,

 - GNU's ``factor`` command won't factor anything greater than 2\ :sup:`127`\ -1.  primefac handles arbitrarily large integers.
 - GNU's ``factor`` command uses Pollard's rho algorithm.  While this extracts small factors quickly, large factors take a while to find.  primefac uses, among other things, the elliptic curve method and the self-initializing quadratic sieve, which are far more efficient at extracting large factors.
 - GNU's ``factor`` command is a single-threaded application.  primefac uses by default five threads to take advantage of the multiple cores typically available on modern machines.  Each of these threads uses a different algorithm to factor the number:

   - One thread runs Brent's variation on Pollard's rho algorithm.  This is good for extracting smallish factors quickly.
   - One thread runs the two-stage version of Pollard's *p*\ -1 method.  This is good at finding factors *p* for which *p*\ -1 is a product of small primes.
   - One thread runs Williams' *p*\ +1 method.  This is good at finding factors *p* for which *p*\ +1 is a product of small primes.
   - One thread runs the elliptic curve method.  This is a bit slower than Pollard's rho algorithm when the factors extracted are small, but it has significantly better performance on difficult factors.
   - One thread runs the self-initializing quadratic sieve.  This is the best algorithm for factoring "hard" numbers short of the general number field sieve.  However, it's (relatively speaking) more than a little slow when the numbers are small, and the time it takes depends only on the size of the number being factored rather than the size of the factors being extracted as with Pollard's rho algorithm and the elliptic curve method, so we use the preceding algorithms to handle those.

 - We also extend the utility by interpreting the command-line arguments as an expression in reverse Polish notation and factoring the numbers remaining on the evaluation stack when interpretation is complete.  For example, the command::

    python3 -m primefac 24 ! 1 - 38 ! 1 +

  will factor the numbers 24! - 1 = 620448401733239439359999 and 38! + 1 = 523022617466601111760007224100074291200000001.


Module Usage
============
The primary functions are ``isprime`` and ``primefac``, but we define a number of helper functions along the way.  A few of these functions are already available via Python3's built-in math module, but get defined anyway for the sake of PyPy3 compatibility.

.. code:: python

    isqrt(n)

Computes the greatest integer whose square does not exceed the non-negative integer ``n``.

.. code:: python

    introot(n, r=2)

For non-negative ``n``, returns the greatest integer less than or equal to the ``r``\ :sup:`th`\  root of ``n``.

For negative ``n``, returns the least integer greater than or equal to the ``r``\ :sup:`th`\  root of ``n``, or ``None`` if ``r`` is even.

.. code:: python

    primegen(limit=inf)

Non-terminating generator.  Yields the prime numbers.  It amounts to a segmented Sieve of Eratosthenes.

.. code:: python

    iterprod(l)

Returns the product of the elements of ``l``, which can be any iterable (but should obviously terminate; e.g., ``iterprod(primegen())`` would be a bad idea).

.. code:: python

    jacobi(a, p)

Computes the Jacobi symbol ``(a|p)``, where ``p`` is a positive odd number.  This is used in ``isprime``.

.. code:: python

    isprime(n, tb=(3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59))

The main primality test.  It iss an implementation of the BPSW test (Baillie-Pomerance-Selfridge-Wagstaff) and is deterministic for all numbers less than 2\ :sup:`64` --- in fact, while infinitely many false positives are conjectured to exist, no false positives are currently known.  We preface the test with some trial division for the sake of speed.  We trial divide by 2 whether the caller wants to or not.

.. code:: python

    ilog(x, b)

Returns the greatest integer ``l`` such that  ``b**l <= x``.

.. code:: python

    ispower(n, r=0)

If ``r == 0``, then we return the largest integer that, when squared/cubed/etc, yields ``n``, or ``None`` if no such integer exists.  If ``r != 0``, then we check whether ``n`` is a perfect ``r``th power, returning its ``r``th root if it is and ``None`` if it is not.

.. code:: python

    pollardrho_brent(n)

Brent's improvement on Pollard's rho algorithm.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    pollard_pm1(n, B1=100, B2=1000)

Pollard's *p*\ +1 algorithm, two-phase version.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    mlucas(v, a, n)

Helper function for ``williams_pp1``.  Multiplies along a Lucas sequence modulo ``n``.

.. code:: python

    williams_pp1(n)

Williams' *p*\ +1 algorithm.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    ecadd(p1, p2, p0, n)

Helper function for ``ecm``.  Adds points on a Montgomery curve modulo ``n``.

.. code:: python

    ecdub(p, A, n)

Helper function for ``ecm``.  Doubles a point on a Montgomery curve modulo ``n``.

.. code:: python

    ecmul(m, p, A, n)

Helper function for ``ecm``.  Multiplies a point on a Montgomery curve modulo ``n``.

.. code:: python

    secm(n, B1, B2, seed)

Seeded ECM.  Helper function for ``ecm``.  Returns a possibly-trivial divisor of n given two bounds and a seed.  Uses the two-phase algorithm on Montgomery curves.  See https://wp.me/prTJ7-zI and https://wp.me/prTJ7-A7 for more details.  Most of the code for this function's "helpers" were shamelessly copied from the first of those links.

.. code:: python

    ecmparams(n)

Generates a sequence of parameters to be applied to ``secm``.

.. code:: python

    ecm(n, paramseq=ecmparams, nprocs=1)

"Modern" integer factoring via elliptic curves.  Uses Montgomery curves, the two-phase algorithm, and (optionally) multiple processes.  The hard work is done by secm(); this function just does the managerial work of pulling a sequence of parameters out of a generator and feeding them into secm().  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.  For more details see `these`_ `two`_ Programming Praxis posts.

.. _these: http://programmingpraxis.com/2010/04/23/modern-elliptic-curve-factorization-part-1/
.. _two: http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/

.. code:: python

    sqrtmod_prime(n, p)

Computes a square root of ``n`` modulo the prime number ``p``.  The return value is not meaningful if ``n`` has no square root modulo ``p`` or if ``p`` is composite.

.. code:: python

    modinv(a, m)

Computes a multiplicative inverse of ``a`` modulo ``m``.  The return value is not meaningful if ``gcd(a,m) != 1``.

.. code:: python

    siqs(n)

Factors ``n`` using the self-initializing quadratic sieve.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.  This function was copied mostly verbatim from `this GitHub repository`__.

__ https://github.com/skollmann/PyFactorise

.. code:: python

    multifactor(n, methods=(pollardrho_brent, pollard_pm1, williams_pp1, ecm, siqs))

Runs several factoring algorithms on ``n`` simultaneously by loading them into their own threads via the ``multiprocessing`` module.  When one function returns, everything is killed off and that value gets returned.

.. code:: python

    primefac(n, trial_limit=1000, rho_rounds=42000, verbose=False,
             methods=(pollardrho_brent))

Generator.  Yields the prime factors of ``n``, with multiplicity.

.. code:: python

    rpn(instr)

Evaluates the string ``instr`` as an arithmetical expression in reverse Polish notation.
The available binary operators are +, -, \*, //, %, and \*\*, which all indicate the same operations here as they indicate in Python3 source code; i.e., they denote addition, subtraction, multiplication, floor division, remaindering, and powering.  The available unary operators are ! and #, which denote the factorial and primorial, respectively.  For terminal syntax compatibility reasons, the RPN expression may be enclosed in quotes, and five aliases are allowed: x for \*, / for //, xx for \*\*, f for !, and p for #.


Dependencies
------------

This package imports items from ``multiprocessing``, ``random``, ``itertools``, ``math``, ``time``, and ``datetime``.  These are all in the Python standard library.


Command-Line Usage
==================

.. code:: sh

    python3 -m primefac [-vs|-sv] [-v|--verbose] [-s|--summary|--summarize] [-t=NUMBER]
                        [-r=NUMBER] [-m=[prb][,p-1][,p+1][,ecm][,siqs]] rpn

``rpn`` is an expression in revese Polish notation and is evaluated using integer arithmetic.  Each number that remains on the stack after evaluation is then factored.

``-t`` sets the trial division limit; the default value is 1000.  Use ``-t=inf`` to use trial division exclusively.

``-r`` sets the number of rounds of Pollard's rho algorithm to try before calling a factor "difficult".  The default value is 42,000.  Use ``-r=inf`` to use Pollard's rho exclusively once the trial division is completed.

If verbosity is invoked, then we provide progress reports and also state which algorithms produced which factors during the multifactor phase.

If the summary and verbosity flags are absent, then the output should be identical to the output of the GNU factor command, modulo permutation of the factors.  If the verbosity flag is invoked, then we provide progress reports, turn on the summary flag, and state which methods yielded which factors during the multifactor phase.

If the summary flag is present, then the output is modified by including a single newline between each item's output, before the first, and after the last.  Each item's output is also modified by printing a second line of data summarizing the results by indicating the number of digits (base 10) in the input, the number of digits (base 10) in each factor, and the factors' multiplicities.  For example:

        >>> user@computer:~$ python3 -m primefac  -s   24 ! 1 -   7 f
        >>> 
        >>> 620448401733239439359999: 991459181683 625793187653
        >>> Z24  =  P12 x P12  =  625793187653 x 991459181683
        >>> 
        >>> 5040: 2 2 2 2 3 3 5 7
        >>> Z4  =  P1^4 x P1^2 x P1 x P1  =  2^4 x 3^2 x 5 x 7
        >>> 
        >>> user@computer:~$

Note that primes in the ordinary output lines are listed in the order in which they were found, while primes in the summary lines are reported in strictly-increasing order.

The ``-v`` and ``-s`` flags may be combined into a single flag in either order --- i.e., into ``-vs`` or ``-sv``.

The ``-m=`` flag controls the functions used during the ``multifactor`` phase.  The options are ``prb``, ``p-1``, ``p+1``, ``ecm``, and ``siqs``, representing Pollard's rho, Pollard's *p*\ -1, Williams' *p*\ +1, the elliptic curve method, and the self-initializing quadratic sieve, respectively.  The options must be separated by commas.  The options can be repeated: if ``prb`` is listed twice, for example, then ``multifactor`` will run two instances of ``pollardrho_brent`` simultaneously.  In the case of ``prb`` and ``ecm``, this decreases the expectation value of the time to find a factor, whereas the other three algorithms (*p*\ -1, *p*\ +1, and MPQS) have no randomized component so that running duplicate instances of these three algorithms confers no benefit.  We therefore ignore repeated listings of the latter three methods: for example, calling

.. code:: sh

    python3 -m primefac -m=prb,prb,ecm,ecm,ecm,mpqs,mpqs 38 ! 1 +

will run during the multifactor phase two instances of Pollard's rho, three instances of the elliptic curve method, and one instance of the MQPS.  Invoking more methods than you have cores available is unlikely to confer any benefit.
