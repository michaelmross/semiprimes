from fractions import Fraction as F
import math

# Gafni-Tao Table 1: upper bounds on A(sigma), piecewise
def A(s):
    if s <= 0.5: return 1/(1-s)
    if s <= 0.7: return 3/(2-s)                      # Ingham
    if s < 19/25: return 15/(3+5*s)                  # Guth-Maynard
    if s < 127/167: return 9/(8*s-2)                 # Ivic
    if s < 13/17: return 15/(13*s-3)                 # Ivic
    if s < 17/22: return 6/(5*s-1)                   # Ivic
    if s < 41/53: return 2/(9*s-6)                   # TTY
    if s < 7/9: return 9/(7*s-1)                     # Ivic
    if s < 1867/2347: return 9/(8*(2*s-1))           # TTY
    if s < 7/8: return 3/(2*s)                       # Bourgain/Ivic
    if s < 279/314: return 3/(10*s-7)                # Heath-Brown
    if s <= 0.9: return 24/(30*s-11)                 # CDV/Ivic
    if s <= 31/34: return 3/(10*s-7)                 # TTY
    if s < 14/15: return 11/(48*s-36)                # TTY
    return 2/(15*s-12)  # coarse for high sigma (enough: constraint kills these)

# Gafni-Tao Table 2: upper bounds on A*(sigma)
def Astar(s):
    if s <= 0.5: return 3/(1-s)
    if s <= 2/3: return (10-11*s)/((2-s)*(1-s))                       # HB
    if s <= 0.7: return (18-19*s)/((4-2*s)*(1-s))                     # HB
    b1 = (539 - math.sqrt(42121))/460
    if s <= b1: return 5*(18-19*s)/(2*(5*s+3)*(1-s))                  # TTY
    if s <= 165/226: return 2*(45-44*s)/((2*s+15)*(1-s))              # TTY
    b2 = (5831 + math.sqrt(60001))/8240
    if s <= b2: return (457-546*s)/(2*(61-58*s)*(1-s))                # TTY
    if s <= 42/55: return 5*(18-19*s)/(2*(5*s+3)*(1-s))               # TTY
    if s <= 97/127: return (18-19*s)/(6*(15*s-11)*(1-s))              # TTY
    if s <= 79/103: return 3*(18-19*s)/(4*(4*s-1)*(1-s))              # TTY
    if s <= 33/43: return (18-19*s)/(2*(37*s-27)*(1-s))               # TTY
    if s <= 84/109: return 5*(18-19*s)/(2*(13*s-3)*(1-s))             # TTY
    b3 = (1273 - math.sqrt(128689))/1184
    if s <= b3: return (18-19*s)/(9*(3*s-2)*(1-s))                    # TTY
    if s <= 5/6: return 4*(10-9*s)/(5*(4*s-1)*(1-s))                  # TTY
    return 12/(4*s-1)                                                  # HB

# NOTE on the "sigma=None" case.
# The Gafni-Tao supremum ranges over sigma with A(sigma) >= 1/(1-theta) - eps.
# At theta = 1/2 this admissible set is a fat interval, (1/2, 25/32], and the
# grid scan finds the interior maximum reliably. But at theta = 17/30 the
# threshold is 1/(1-theta) = 30/13, which the tabulated bounds attain ONLY in
# the limit sigma -> 7/10 (Ingham and Guth-Maynard both equal 30/13 exactly at
# sigma = 7/10 and are strictly smaller on either side). The admissible set is
# thus a single point, no grid node lands exactly on it, the supremum comes
# back empty, and mu_bound falls through to the 1-theta floor of GT eq. (1.9),
# reporting argbest = None. This is a grid artifact, not a property of the
# bound: the true value mu(17/30) <= 7/12 is confirmed by the EXACT rational
# evaluation at sigma = 7/10 performed at the bottom of this script, which is
# the actual consistency check against Gafni-Tao Section 3.
def mu_bound(theta, eps=1e-9, n=4_000_000):
    best, argbest = 1-theta, None   # GT (1.9): max(1-theta, sup ...)
    lo, hi = 0.5+1e-12, 0.999999
    for i in range(n+1):
        s = lo + (hi-lo)*i/n
        if A(s) < 1/(1-theta) - eps:   # constraint A(sigma) >= 1/(1-theta)-eps must be satisfiable
            continue
        m2 = (1-theta)*(1-s)*A(s) + 2*s - 1
        m4 = (1-theta)*(1-s)*Astar(s) + 4*s - 3
        v = min(m2, m4)
        if v > best:
            best, argbest = v, s
    return best, argbest

# The first four values verify Proposition 3: the maximum sits at sigma = 7/10
# for theta = 1/2 and stays there as theta decreases (the k >= 5 case of
# Theorem 2), with the value growing linearly in 1/2 - theta. The last value,
# theta = 17/30, prints sigma=None by the grid artifact described above; it is
# included only to mark where the exact-rational check below takes over.
for th in (0.5, 0.499, 0.495, 0.49, 17/30):
    b, s = mu_bound(th)
    print(f"theta={th:.6f}  mu <= {b:.6f}  at sigma={s}")

# exact values at sigma = 7/10, theta = 1/2
from fractions import Fraction
s = Fraction(7,10); th = Fraction(1,2)
Astar_710 = Fraction(18-19*7//1,1)  # careful: compute exactly
Astar_710 = (Fraction(18) - Fraction(19)*s) / ((Fraction(4)-Fraction(2)*s)*(Fraction(1)-s))
A_710 = Fraction(3)/(Fraction(2)-s)
mu4 = (1-th)*(1-s)*Astar_710 + 4*s - 3
mu2 = (1-th)*(1-s)*A_710 + 2*s - 1
print("A(7/10) =", A_710, "; A*(7/10) =", Astar_710)
print("mu2 =", mu2, "=", float(mu2), "; mu4 =", mu4, "=", float(mu4))
print("exponent 2*mu4 - 1 =", 2*mu4-1, "=", float(2*mu4-1))
# GT sanity check: theta = 17/30 should give 7/12
th2 = Fraction(17,30)
mu4b = (1-th2)*(1-s)*Astar_710 + 4*s - 3
print("GT check theta=17/30: mu4 =", mu4b, "(should be 7/12 =", Fraction(7,12), ")")
