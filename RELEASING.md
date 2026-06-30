# Releasing & archiving (Zenodo DOI)

To give the pipeline a citable DOI for the bioRxiv submission:

1. **Connect the repo to Zenodo** (one time): sign in at <https://zenodo.org> with
   your GitHub account, go to *Settings → GitHub*, and flip the toggle **on** for
   `bennettek99-spec/Archaic-DNA-processing-pipeline`.
2. **Tag a release on GitHub**: Releases → *Draft a new release* → tag e.g.
   `v0.2.0` → Publish. Zenodo automatically archives that tarball and mints a DOI.
   (`.zenodo.json` in this repo supplies the title, authors, license and keywords.)
3. **Cite it**: Zenodo gives both a version DOI and a permanent "all-versions" DOI.
   Put the all-versions DOI in the paper's *Code availability* section and in
   `CITATION.cff`.

The repo already ships `CITATION.cff` (GitHub shows a "Cite this repository"
button) and `.zenodo.json` (Zenodo metadata). CI (`.github/workflows/ci.yml`) runs
the unit tests on every push so the archived release is known-green.

## Versioning
Bump `version` in `pyproject.toml`, `CITATION.cff`, and `.zenodo.json` together,
then tag.
