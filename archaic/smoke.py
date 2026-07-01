"""
smoke.py — synthetic-data smoke test for the shared engine every phase builds
on: PackedGeno (file format) -> Panel (frequencies) -> stats (f4-ratio/D).

This is deliberately SECONDARY to the real AADR-based pipeline: it never reads
or writes anything under the configured aadr_dir, and it does not validate
scientific accuracy (see validate_simulation.py / VALIDATION.md for that). Its
only job is to catch plumbing regressions (a broken reader, a sign flip, a
reference mix-up) in seconds, without needing the ~4 GB non-redistributable
AADR files — useful in CI and for anyone who clones the repo without a copy of
the data.

    python -m archaic.smoke          # or: archaic-pipeline smoke-test
"""
from __future__ import annotations
import shutil
import tempfile

from .synthetic import write_synthetic_panel
from .panel import Panel
from . import stats as st
from .log_utils import get_logger

log = get_logger("archaic.smoke")

N_SNP = 30_000
N_BLOCKS = 50
ALPHA_TRUE = 0.05


def run_smoke_test(seed: int = 1, verbose: bool = True) -> bool:
    """Build a synthetic panel, run it through Panel + stats, and check the
    estimator behaves sanely. Returns True iff every check passes."""
    tmp = tempfile.mkdtemp(prefix="archaic_smoke_")
    checks = []
    try:
        prefix, info = write_synthetic_panel(
            tmp, n_snp=N_SNP, alpha_true=ALPHA_TRUE, seed=seed)
        panel = Panel(prefix, autosomes_only=True)
        block = st.assign_blocks(panel.n_snp, N_BLOCKS)

        freq, finfo = panel.frequencies({
            "Altai": dict(ids=["Altai"]), "Vindija": dict(ids=["Vindija"]),
            "Denisova": dict(ids=["Denisova"]), "Chimp": dict(ids=["Chimp"]),
            "Mbuti": dict(pops=["Mbuti"]), "Test": dict(pops=["Test"]),
            "AltaiEcho": dict(ids=["AltaiEcho"]),
        })

        checks.append(("panel loads with the right shape",
                        f"n_snp={panel.n_snp:,} n_ind={len(panel.ind)}",
                        panel.n_snp == N_SNP and len(panel.ind) == sum(info["n_per_pop"].values())))

        echo = st.f4_ratio(freq, "Altai", "Chimp", "AltaiEcho", "Mbuti", "Vindija",
                           block, N_BLOCKS)
        checks.append(("archaic-as-its-own-test reads ~100% Neanderthal",
                        f"alpha={echo['theta']*100:.0f}% z={echo['z']:.1f}",
                        0.6 <= echo["theta"] <= 1.4 and echo["z"] > 5))

        test = st.f4_ratio(freq, "Altai", "Chimp", "Test", "Mbuti", "Vindija",
                           block, N_BLOCKS)
        checks.append(("synthetic Test population recovers a positive, "
                        "significant Neanderthal signal",
                        f"alpha={test['theta']*100:.1f}% (true {ALPHA_TRUE*100:.0f}%) "
                        f"z={test['z']:.1f}",
                        test["z"] > 1.5 and 0 < test["theta"] < 0.3))

        dnea = st.dstat(freq, "Test", "Mbuti", "Altai", "Chimp", block, N_BLOCKS)
        checks.append(("D_Nea(Test,Mbuti;Altai,Chimp) is positive (Test carries "
                        "more Neanderthal than the non-introgressed baseline)",
                        f"D={dnea['theta']:.4f} z={dnea['z']:.1f}",
                        dnea["theta"] > 0))

        dden = st.dstat(freq, "Test", "Mbuti", "Denisova", "Chimp", block, N_BLOCKS)
        checks.append(("D_Den(Test,Mbuti;Denisova,Chimp) is small (Test has no "
                        "simulated Denisovan admixture)",
                        f"D={dden['theta']:.4f} z={dden['z']:.1f}",
                        abs(dden["z"]) < 4))

        del panel
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    allpass = all(ok for _, _, ok in checks)
    if verbose:
        log.info("SYNTHETIC SMOKE TEST (secondary to the AADR-based pipeline)")
        for label, detail, ok in checks:
            log.info(f"  [{'PASS' if ok else 'FAIL'}] {label:60s} {detail}")
        log.info(f"OVERALL: {'PASS' if allpass else 'FAIL'}")
    return allpass


if __name__ == "__main__":
    import sys
    sys.exit(0 if run_smoke_test() else 1)
