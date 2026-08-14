"""A Skov-type two-state HMM over outgroup-private variant counts.

The Skov et al. (2018) caller bins a genome into fixed windows, counts the
variants an ingroup individual carries that are absent from an unadmixed
outgroup, and decodes a two-state hidden Markov model with Poisson emissions:
a low-rate modern-human state and a high-rate archaic state.

This module reimplements that model so the caller can be run on simulated
genotypes with known truth. That is the only way to measure how posterior
decoding distorts run lengths without the circularity of rescaling a decoded
distribution by the parameter it is supposed to test.

Two quantities come out of a fit and they are not the same number:

``fitted_generations``
    Derived from the archaic-to-human transition probability. This is the
    model's own admixture-time parameter, comparable to the Skov S4 column.
``decoded_generations``
    The exponential decay of the archaic runs that posterior decoding actually
    returns, comparable to what a downstream tract analysis would measure.

Their ratio is the decoder inflation. On the real 89 Papuan individuals it is
1.56x; reproducing that on simulated data is what validates this
reimplementation.

The forward-backward recursion is a scalar loop over millions of windows, so it
is JIT-compiled. Without that the calibration is not laptop-feasible.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

import numpy as np

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only in stripped envs
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def wrap(func):
            return func

        if args and callable(args[0]):
            return args[0]
        return wrap


DEFAULT_WINDOW_BP = 1000
DEFAULT_RECOMBINATION_RATE = 1.2e-8


# --------------------------------------------------------------------------
# emissions
# --------------------------------------------------------------------------


def poisson_emissions(obs: np.ndarray, rates: np.ndarray) -> np.ndarray:
    """Poisson emission likelihoods for every window and state.

    Counts are small non-negative integers, so the pmf is evaluated once per
    distinct count and gathered, which is far cheaper than evaluating it per
    window.
    """
    from scipy.special import gammaln

    obs = np.asarray(obs, dtype=np.int64)
    rates = np.asarray(rates, dtype=np.float64)
    kmax = int(obs.max()) if obs.size else 0
    ks = np.arange(kmax + 1, dtype=np.float64)
    log_table = (
        ks[:, None] * np.log(rates)[None, :]
        - rates[None, :]
        - gammaln(ks + 1.0)[:, None]
    )
    return np.exp(log_table)[obs]


# --------------------------------------------------------------------------
# scaled forward-backward
# --------------------------------------------------------------------------


@njit(cache=True, fastmath=True)
def _forward_backward(e, A, pi):
    n = e.shape[0]
    k = e.shape[1]

    alpha = np.empty((n, k))
    scale = np.empty(n)

    total = 0.0
    for s in range(k):
        alpha[0, s] = pi[s] * e[0, s]
        total += alpha[0, s]
    scale[0] = total
    for s in range(k):
        alpha[0, s] /= total

    for t in range(1, n):
        total = 0.0
        for s in range(k):
            acc = 0.0
            for r in range(k):
                acc += alpha[t - 1, r] * A[r, s]
            alpha[t, s] = acc * e[t, s]
            total += alpha[t, s]
        scale[t] = total
        for s in range(k):
            alpha[t, s] /= total

    beta = np.empty((n, k))
    for s in range(k):
        beta[n - 1, s] = 1.0
    for t in range(n - 2, -1, -1):
        inv = 1.0 / scale[t + 1]
        for s in range(k):
            acc = 0.0
            for r in range(k):
                acc += A[s, r] * e[t + 1, r] * beta[t + 1, r]
            beta[t, s] = acc * inv

    gamma = np.empty((n, k))
    for t in range(n):
        total = 0.0
        for s in range(k):
            gamma[t, s] = alpha[t, s] * beta[t, s]
            total += gamma[t, s]
        for s in range(k):
            gamma[t, s] /= total

    xi = np.zeros((k, k))
    for t in range(n - 1):
        inv = 1.0 / scale[t + 1]
        for s in range(k):
            for r in range(k):
                xi[s, r] += alpha[t, s] * A[s, r] * e[t + 1, r] * beta[t + 1, r] * inv

    loglik = 0.0
    for t in range(n):
        loglik += np.log(scale[t])

    return gamma, xi, loglik


@njit(cache=True, fastmath=True)
def _weighted_rate_update(obs, gamma):
    k = gamma.shape[1]
    num = np.zeros(k)
    den = np.zeros(k)
    for t in range(obs.shape[0]):
        for s in range(k):
            num[s] += gamma[t, s] * obs[t]
            den[s] += gamma[t, s]
    return num / den


# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------


@dataclass
class HMMFit:
    rates: np.ndarray
    transitions: np.ndarray
    initial: np.ndarray
    loglik: float
    iterations: int
    converged: bool
    archaic_fraction: float
    window_bp: int
    recombination_rate: float

    @property
    def fitted_generations(self) -> float:
        """Admixture time implied by the archaic-to-human transition rate.

        A tract laid down ``t`` generations ago is broken by recombination at
        rate ``t * r`` per base pair, so the per-window probability of leaving
        the archaic state is ``t * r * window_bp``.
        """
        p_leave = float(self.transitions[1, 0])
        return p_leave / (self.recombination_rate * self.window_bp)

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["rates"] = self.rates.tolist()
        out["transitions"] = self.transitions.tolist()
        out["initial"] = self.initial.tolist()
        out["fitted_generations"] = self.fitted_generations
        return out


def _initial_parameters(
    obs: np.ndarray,
    archaic_fraction: float,
    admixture_generations: float,
    window_bp: int,
    recombination_rate: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = float(obs.mean()) if obs.size else 1.0
    # Start the two states well separated; EM pulls them to the data.
    rates = np.array([mean * 0.5, mean * 4.0], dtype=np.float64)
    rates = np.maximum(rates, 1e-6)

    p_leave = admixture_generations * recombination_rate * window_bp
    p_leave = min(max(p_leave, 1e-6), 0.5)
    p_enter = p_leave * archaic_fraction / (1.0 - archaic_fraction)
    p_enter = min(max(p_enter, 1e-9), 0.5)

    A = np.array([[1.0 - p_enter, p_enter], [p_leave, 1.0 - p_leave]])
    pi = np.array([1.0 - archaic_fraction, archaic_fraction])
    return rates, A, pi


def fit_hmm(
    obs: np.ndarray,
    *,
    archaic_fraction: float = 0.05,
    admixture_generations: float = 1500.0,
    window_bp: int = DEFAULT_WINDOW_BP,
    recombination_rate: float = DEFAULT_RECOMBINATION_RATE,
    max_iter: int = 50,
    tol: float = 1e-4,
) -> HMMFit:
    """Baum-Welch fit of the two-state Poisson HMM.

    ``archaic_fraction`` and ``admixture_generations`` only seed the search;
    both are free parameters of the fit.
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be at least 1, got {max_iter}")
    obs = np.ascontiguousarray(np.asarray(obs, dtype=np.int64))
    if obs.size < 2:
        raise ValueError(f"need at least 2 windows to fit, got {obs.size}")
    rates, A, pi = _initial_parameters(
        obs, archaic_fraction, admixture_generations, window_bp, recombination_rate
    )

    previous = -np.inf
    converged = False
    iterations = 0
    gamma = None

    for iterations in range(1, max_iter + 1):
        e = poisson_emissions(obs, rates)
        gamma, xi, loglik = _forward_backward(e, A, pi)

        A = xi / xi.sum(axis=1, keepdims=True)
        rates = _weighted_rate_update(obs, gamma)
        rates = np.maximum(rates, 1e-9)
        pi = gamma[0].copy()

        # Keep state 1 as the archaic (high-rate) state so downstream code can
        # rely on the ordering rather than re-deriving it.
        if rates[0] > rates[1]:
            rates = rates[::-1].copy()
            A = A[::-1, ::-1].copy()
            pi = pi[::-1].copy()
            gamma = gamma[:, ::-1].copy()

        if abs(loglik - previous) < tol * max(1.0, abs(previous)):
            converged = True
            previous = loglik
            break
        previous = loglik

    fraction = float(gamma[:, 1].mean()) if gamma is not None else float("nan")
    return HMMFit(
        rates=rates,
        transitions=A,
        initial=pi,
        loglik=float(previous),
        iterations=iterations,
        converged=converged,
        archaic_fraction=fraction,
        window_bp=window_bp,
        recombination_rate=recombination_rate,
    )


def posterior_decode(obs: np.ndarray, fit: HMMFit) -> np.ndarray:
    """Posterior state probabilities for the archaic state, per window.

    Returns probabilities rather than calls; thresholding belongs to
    :func:`extract_runs`, so that one decode can be re-thresholded cheaply.
    """
    obs = np.ascontiguousarray(np.asarray(obs, dtype=np.int64))
    e = poisson_emissions(obs, fit.rates)
    gamma, _, _ = _forward_backward(e, fit.transitions, fit.initial)
    return gamma[:, 1]


def extract_runs(
    posterior: np.ndarray, threshold: float = 0.5
) -> tuple[np.ndarray, np.ndarray]:
    """Maximal runs of archaic-state windows.

    Returns half-open ``[start, end)`` window indices.
    """
    called = posterior >= threshold
    if not called.any():
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    padded = np.concatenate(([False], called, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return edges[0::2].astype(np.int64), edges[1::2].astype(np.int64)


# --------------------------------------------------------------------------
# decay estimation
# --------------------------------------------------------------------------


def decay_generations(
    lengths_morgans: np.ndarray, min_length_morgans: float = 0.0
) -> float:
    """MLE exponential decay rate of left-truncated tract lengths.

    For lengths drawn from an exponential with rate ``lam`` and observed only
    above ``T``, the MLE is ``1 / mean(x - T)``. The rate is in generations
    because tract length in Morgans decays at ``t`` per generation.
    """
    lengths = np.asarray(lengths_morgans, dtype=np.float64)
    kept = lengths[lengths >= min_length_morgans]
    if kept.size == 0:
        return float("nan")
    excess = kept.mean() - min_length_morgans
    if excess <= 0:
        return float("nan")
    return 1.0 / excess


def run_lengths_morgans(
    starts: np.ndarray,
    ends: np.ndarray,
    window_bp: int = DEFAULT_WINDOW_BP,
    recombination_rate: float = DEFAULT_RECOMBINATION_RATE,
) -> np.ndarray:
    return (ends - starts).astype(np.float64) * window_bp * recombination_rate


def call_individual(
    obs: np.ndarray,
    *,
    window_bp: int = DEFAULT_WINDOW_BP,
    recombination_rate: float = DEFAULT_RECOMBINATION_RATE,
    min_length_morgans: float = 0.0,
    threshold: float = 0.5,
    **fit_kwargs: Any,
) -> dict[str, Any]:
    """Fit, decode, and summarise one individual's window counts."""
    fit = fit_hmm(
        obs,
        window_bp=window_bp,
        recombination_rate=recombination_rate,
        **fit_kwargs,
    )
    posterior = posterior_decode(obs, fit)
    starts, ends = extract_runs(posterior, threshold=threshold)
    lengths = run_lengths_morgans(starts, ends, window_bp, recombination_rate)
    decoded = decay_generations(lengths, min_length_morgans)
    fitted = fit.fitted_generations
    return {
        "fit": fit,
        "starts": starts,
        "ends": ends,
        "lengths_morgans": lengths,
        "n_runs": int(starts.size),
        "decoded_generations": decoded,
        "fitted_generations": fitted,
        "decoded_over_fitted": decoded / fitted if fitted > 0 else float("nan"),
        "decoded_archaic_fraction": float((ends - starts).sum()) / float(obs.size),
    }
