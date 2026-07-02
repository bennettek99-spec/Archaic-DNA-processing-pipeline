"""
ancestry.py — general-purpose West-Eurasian ancestry decomposition.

Where the archaic side of this pipeline asks "how much Neanderthal/Denisovan?",
this module asks the complementary question: "of what *human* source populations
is a genome a mixture?" It models any target cohort as a mixture of the canonical
ancient source populations that structure West-Eurasian variation —

    WHG   Western Hunter-Gatherers      (Loschbour / Villabruna cline)
    EHG   Eastern Hunter-Gatherers      (Karelia, Samara)
    CHG   Caucasus Hunter-Gatherers     (Kotias, Satsurblia)
    Iran_N Iranian Neolithic            (Ganj Dareh)
    Anatolia_N  Anatolian Neolithic     (early European farmer / EEF proxy)
    Levant_N / Natufian                 (basal-Eurasian-rich Near East)
    Steppe_Yamnaya  steppe pastoralists (EHG + CHG; the Bronze-Age expansion)
    ANE   Ancient North Eurasians       (Mal'ta MA-1, Afontova Gora)

— using the pipeline's own f4-statistic machinery (archaic.stats) and qpAdm
(archaic.qpadm), with block-jackknife SEs. Both the classic unconstrained qpAdm
(Haak et al. 2015) and a simplex-constrained "supervised admixture" fit are run,
so every target gets an interpretable, valid set of proportions even where the
unconstrained model strays off the simplex.

All source/outgroup definitions below were verified against the local AADR v66.1
1240K release; missing labels are dropped gracefully at run time.
"""
from __future__ import annotations
from collections import OrderedDict
import numpy as np

from . import qpadm as qp
from . import kinship as kin
from . import profiles as pf


# --------------------------------------------------------------- source library
# Each source is a predicate on the *lower-cased* AADR group_id. Order is the
# canonical West->East / HG->farmer->steppe ordering used for stacked bars.
SOURCES = OrderedDict([
    ("WHG", dict(
        color="#3b82c4",
        desc="Western Hunter-Gatherers (Villabruna cluster)",
        pred=lambda g: (("mesolithic" in g or g.endswith("_hg") or "_hg-" in g
                         or "_hg_" in g or "loschbour" in g)
                        and any(c in g for c in (
                            "france", "spain", "england", "belgium", "netherlands",
                            "germany", "luxembourg", "switzerland", "italy",
                            "iberia", "croatia", "poland"))
                        and "irongates" not in g))),
    ("EHG", dict(
        color="#7ab648",
        desc="Eastern Hunter-Gatherers (Karelia, Samara)",
        pred=lambda g: any(s in g for s in (
            "karelia_mesolithic", "russia_samara_hg", "russia_samara_en_mesolithic",
            "russia_vologda_mesolithic", "sidelkino", "russia_veretye")))),
    ("CHG", dict(
        color="#c0392b",
        desc="Caucasus Hunter-Gatherers (Kotias, Satsurblia)",
        pred=lambda g: "kotiasklde" in g or "satsurblia" in g)),
    ("Iran_N", dict(
        color="#e08e0b",
        desc="Iranian Neolithic (Ganj Dareh)",
        pred=lambda g: "iran_ganjdareh_n" in g)),
    ("Anatolia_N", dict(
        color="#f1c40f",
        desc="Anatolian Neolithic farmers (EEF proxy)",
        pred=lambda g: "anatolia_n" in g or "turkey_n" in g)),
    ("Levant_N", dict(
        color="#d98880",
        desc="Levantine Neolithic (PPNB)",
        pred=lambda g: "jordan_ppnb" in g or ("israel" in g and "ppn" in g))),
    ("Natufian", dict(
        color="#a04000",
        desc="Epipalaeolithic Natufians (Levant)",
        pred=lambda g: "natufian" in g)),
    ("Steppe_Yamnaya", dict(
        color="#8e44ad",
        desc="Steppe pastoralists (Yamnaya)",
        pred=lambda g: "yamnaya" in g)),
    ("ANE", dict(
        color="#16a085",
        desc="Ancient North Eurasians (Mal'ta MA-1, Afontova Gora)",
        pred=lambda g: "russia_malta_up" in g or "afontovagora" in g)),
])


def source_color(name):
    return SOURCES.get(name, {}).get("color", "#888888")


# ---------------------------------------------------------- outgroup ("right") set
# Distal outgroups that are always safe (never a West-Eurasian source), verified
# present in the 1240K panel.  name -> ("pop"|"id", value).
BASE_RIGHT = OrderedDict([
    ("Mbuti", ("pop", "Mbuti")),
    ("Han", ("pop", "Han")),
    ("Papuan", ("pop", "Papuan")),
    ("Karitiana", ("pop", "Karitiana")),
    ("Ust_Ishim", ("id", "Ust_Ishim.DG")),
    ("Kostenki14", ("id", "Kostenki14.SG")),
    ("MA1", ("id", "MA1.SG")),
])
# Near-source outgroups that sharpen basal-Eurasian / HG resolution. They ARE
# entries of the source library (same frequency key), so they cost no extra read
# and are simply dropped from the right set when they double as a model source.
EXTRA_RIGHT = ["Natufian", "Iran_N", "CHG", "EHG"]


# --------------------------------------------------------------- model presets
MODELS = OrderedDict([
    ("west3", ["Anatolia_N", "Steppe_Yamnaya", "WHG"]),
    ("west4", ["Anatolia_N", "Steppe_Yamnaya", "WHG", "Iran_N"]),
    ("deep4", ["WHG", "EHG", "CHG", "Anatolia_N"]),
    ("deep5", ["WHG", "EHG", "CHG", "Anatolia_N", "Iran_N"]),
    ("hg3",   ["WHG", "EHG", "CHG"]),
])


def default_outgroups(sources):
    """Ordered list of outgroup *names* (frequency keys) for a model: the always-
    safe distal set plus any near-source outgroup not itself a model source."""
    srcs = set(sources)
    names = list(BASE_RIGHT.keys())
    names += [x for x in EXTRA_RIGHT if x not in srcs]
    return names


# ---------------------------------------------------------------- cohort resolve
def _mask_for(meta_gl, pred):
    return meta_gl.map(pred).to_numpy()


def resolve_cohorts(panel, meta, specs, maxn=60, kin_prune=True, kin_snp=40000,
                    seed=0, verbose=False):
    """Turn cohort specifications into arrays of .geno column indices.

    specs : name -> spec, where spec is one of
        callable(g)             predicate on lower-cased group_id  (ancient set)
        ("grp", substring)      ancient group_id substring match
        ("pop", label)          present-day / labelled population from the panel
        ("id",  genetic_id)     a single named individual
    Cohorts are capped to `maxn` individuals, then optionally kinship-pruned
    (archaic.kinship) so relatives do not bias the mean frequencies.
    Returns name -> np.int64 column array (empty if nothing matched).
    """
    gl = meta["group_id"].str.lower()
    id_col = panel._id_to_col
    rng = np.random.default_rng(seed)
    out = {}
    for name, spec in specs.items():
        if callable(spec):
            ids = meta.loc[_mask_for(gl, spec), "genetic_id"].tolist()
            cols = np.array([id_col[i] for i in ids if i in id_col], dtype=np.int64)
        elif isinstance(spec, tuple) and spec[0] == "grp":
            m = gl.str.contains(spec[1], regex=False, na=False).to_numpy()
            ids = meta.loc[m, "genetic_id"].tolist()
            cols = np.array([id_col[i] for i in ids if i in id_col], dtype=np.int64)
        elif isinstance(spec, tuple) and spec[0] == "pop":
            cols = panel.cols_for(pops=[spec[1]])
        elif isinstance(spec, tuple) and spec[0] == "id":
            cols = panel.cols_for(ids=[spec[1]])
        else:
            cols = np.array([], dtype=np.int64)

        if len(cols) > maxn:
            cols = np.sort(rng.choice(cols, maxn, replace=False))
        if kin_prune and len(cols) >= 4:
            try:
                keep, _dropped, _ = kin.prune(panel, cols, n_snp=kin_snp)
                cols = keep
            except Exception:
                pass
        out[name] = cols
        if verbose:
            print(f"  {name:20s} n={len(cols)}")
    return out


def cohort_freqs(panel, cohort_cols):
    """mean-genome allele frequencies for cohorts with >=1 individual (one read
    each, memory-safe). Returns (freq, info)."""
    nz = {k: v for k, v in cohort_cols.items() if len(v)}
    return pf.cohort_frequencies(panel, nz)


# ------------------------------------------------------------------- decompose
def decompose(freq, target, sources, outgroups, block, n_blocks=50):
    """Both qpAdm fits (unconstrained + simplex-constrained) for one target.
    Returns dict with weights/se/p for each, plus the sources actually used
    (those present in `freq`) and the outgroups used."""
    srcs = [s for s in sources if s in freq and np.isfinite(freq[s]).any()]
    outs = [o for o in outgroups if o in freq and np.isfinite(freq[o]).any()]
    if len(srcs) < 2 or len(outs) < len(srcs):
        return dict(target=target, sources=srcs, outgroups=outs, ok=False)
    r_free = qp.qpadm(freq, target, srcs, outs, block, n_blocks)
    r_con = qp.qpadm_constrained(freq, target, srcs, outs, block, n_blocks)
    return dict(target=target, sources=srcs, outgroups=outs, ok=True,
                free=r_free, constrained=r_con,
                n_snp=r_con["n_snp"])


def weights_dict(result, which="constrained"):
    """Flatten a decompose() result into {source: (weight, se)} for the chosen fit."""
    r = result.get(which)
    if not result.get("ok") or r is None:
        return {}
    return {s: (w, se) for s, w, se in zip(r["sources"], r["weights"], r["se"])}


def decompose_best(freq, target, models, block, n_blocks=50):
    """Run every candidate model in `models` (name -> source-name list) against
    one target, each with its own methodologically-valid outgroup set (a source
    is never also used as its own outgroup — see default_outgroups), and rank
    them best-first. Ranking: an unconstrained fit that already lands on the
    simplex (no negative/>1 weights) is preferred, then higher GLS fit p-value
    (more plausible), then fewer sources (Occam's razor).

    Returns a list of decompose()-shaped dicts, each tagged with 'model'."""
    results = []
    for name, srcs in models.items():
        outs = default_outgroups(srcs)
        r = decompose(freq, target, srcs, outs, block, n_blocks)
        r["model"] = name
        results.append(r)

    def rank_key(r):
        if not r["ok"]:
            return (1, 0.0, 99)
        p = r["free"]["p"]
        return (0 if r["free"]["feasible"] else 0.5,
                -(p if np.isfinite(p) else -1.0),
                len(r["sources"]))
    results.sort(key=rank_key)
    return results
