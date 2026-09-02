# Implicit-CNT axial VACF: discrete longitudinal-mode reconstruction protocol

**Status:** current C88 low-k reconstruction protocol, validated by the
synthetic-N800 to N1600 blind back-test on 2026-09-01.  This is a finite,
discrete-mode construction; it is not a continuum identity for a complete
microscopic VACF and does not supply an unmeasured high-k/local term.

## 1. Scope and observable

For a periodic axial length `Lz`, use positive discrete axial wavevectors

`k_n = 2*pi*n/Lz`, `n=1,...,M`.

For each velocity seed `r`, measure at the same physical `k_n`:

- `W_n^(r) = C_JJ(k_n,0) / [N_O <v_z^2>]`;
- `Fs_n^(r)(t)`, the self intermediate scattering factor;
- `Phi_n^(r)(t) = C_JJ(k_n,t)/C_JJ(k_n,0)`, the normalized longitudinal
  current carrier.

The physical positive/negative-k paired discrete sum is

`P^(r)(t) = 2 sum_(n=1)^M W_n^(r) Fs_n^(r)(t) Phi_n^(r)(t)`.

Perform this sum for each seed first, then calculate the ensemble mean and
conditional velocity-seed SEM.  Do not average fields across seeds before
multiplication.  The factor 2 is the real-field `+/-k` degeneracy factor.

## 2. Current carrier parameterization

At every measured source mode and seed fit the *real* carrier to

`Phi_n(t) = exp[-Gamma_n t] [ a_n cos(omega_n t) + b_n sin(omega_n t) ]`.

`omega_n`, `Gamma_n`, `a_n`, and `b_n` are separate fitted fields.  The
fit is a representation of the measured real autocorrelation, not a claim
that a single unshifted DHO is exact.  In particular, do not impose `b=0`,
and do not replace an observed matched-k carrier by its fitted curve.

Because the carrier is normalized, fitted `a` should be audited against its
expected near-unity value; it is not a global VACF amplitude knob.  `b` is a
phase/initial-slope parameter, not a substitute for a local baseline.

## 3. Cross-box mapping: source box L to target box 2L

Let the source be N800 and target be N1600 in the diagnostic, or N1600 and
target N3200 in a prospective application.  At a shared physical wavevector,
target even `n=2m` maps exactly to source `m`.

### A. Exact matched-mode rule (non-negotiable)

For every even target mode, retain the full source field at every time lag:

`Phi_target,2m(t) = Phi_source,m(t)`;

`Fs_target,2m(t) = Fs_source,m(t)`;

`W_target,2m = W_source,m / 2`.

The last relation is the discrete per-mode scaling associated with the
doubling of axial length, consistent with approximately constant `N W`.
Never regenerate a matched carrier from interpolated DHO parameters.

### B. Newly inserted odd modes

For target odd `n=3,5,...`, independently interpolate in physical `k`:

`omega(k), Gamma(k), a(k), b(k)`.

Construct only the missing carrier using the phase-DHO form above.  At every
lag independently interpolate `log Fs(k,t)` (positive field).  Interpolate
the static `W(k)` separately.  The product `W Fs Phi` must not be interpolated
as a single opaque kernel because it hides the source of error.

### C. New lowest mode n=1

Treat `k_1` as a distinct low-k extrapolation, not merely an odd interpolation.
For the current diagnostic, independently extrapolate `omega,a,b` from the
first four source modes without imposing `omega=ck`.  The default damping
continuation is the no-zero-wave-number-damping form

`Gamma(k) = A k^p`, with `A>0` and `p>0`,

obtained by a low-k log-log fit.  Do not use a linear fit of `log Gamma`
against `k` as the default: it has a finite `Gamma(k=0)` intercept and, in the
N800-to-N1600 held-out n=1 test, over-damped the new lowest mode by about a
factor of three.  A finite `Gamma(0)` is an exception that must be established
by an explicit independent control, not an unconstrained extrapolation
artifact.  Use `log Fs` versus `k^2` for each lag and extrapolate `W` versus
`k^2` as a separate static operation.

This rule is provisional and must be blind-back-tested whenever a smaller
physical k becomes measurable.  No claim of a nonzero `Gamma(k=0)` is made
unless direct low-k data establish it.

## 4. Required validation hierarchy

1. **Carrier check:** source-parameter fit quality; show `omega,Gamma,a,b`
   versus k and representative carrier traces.
2. **Exact-anchor check:** all matched-even target modes must be reproduced
   exactly from serialized source fields, before testing any missing modes.
3. **Odd-mode check:** compare inserted modes with held-out target data,
   mode by mode, including first-zero timing and phase morphology.
4. **New-k1 check:** report it separately; do not conceal its error in a
   global statistic.
5. **Field-isolation check:** compare Phi-only, Phi+Fs, and Phi+Fs+W sums
   against held-out truth, while retaining the provenance of each field.
6. **Residual localization:** decompose `Delta K_n(t) = 2 W_n Fs_n
   (Phi_n^pred-Phi_n^truth)` to identify the mode(s) responsible for each
   VACF lobe residual.

## 5. Prohibited shortcuts

- No local/residual/taper/direct-VACF baseline term.
- No fitted global amplitude multiplying the final sum.
- No time-window-dependent reassignment of a residual to high-k modes.
- No `omega=ck` forced through the origin merely because finite-k data appear
  nearly linear.
- No averaging across velocity seeds before forming the nonlinear field product.
- No treating a synthetic coarse box derived from a fine box as an independent
  physical validation.

## 6. Evidence and current limits

The N800-to-N1600 synthetic blind back-test has shown that retaining exact
anchors plus phase-DHO interpolation materially repairs the prior failure.
It is an algorithmic closure test only because the synthetic N800 source was
made from N1600 even modes.  The remaining sensitivity is dominated by the
new lowest k carrier.  Independent N800/N1600/N3200 trajectories remain
necessary for physical validation and for any asserted low-k dispersion law.

Authoritative asset packages:

- `results/collective_mode_response/C88_N1600_syntheticN800_roundtrip_phaseDHO_ab/2026-09-01/`
- `results/collective_mode_response/C88_N1600_syntheticN800_roundtrip_phaseDHO_Fs_W/2026-09-01/`
