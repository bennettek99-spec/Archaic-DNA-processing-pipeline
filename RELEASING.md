# Releasing and archiving

Releases are created from a green `main` branch and archived through Zenodo when
the repository integration is enabled.

## Release checklist

1. Move completed changelog entries from `Unreleased` into a dated version.
2. Set the same version in:
   - `pyproject.toml`
   - `archaic/__init__.py`
   - `CITATION.cff`
   - `.zenodo.json`
3. Update `date-released` in `CITATION.cff` and `publication_date` in
   `.zenodo.json`.
4. Run:

   ```bash
   python tools/check_repo_docs.py
   pyflakes archaic archaic_admixture_dating scripts tests tools oase1_bam_pipeline
   python -m pytest -q
   archaic-pipeline smoke-test
   python -m build
   ```

5. Merge the release pull request after every required CI check passes.
6. Create an annotated `vX.Y.Z` tag on the merge commit and push it.
7. Create the matching GitHub release from the changelog entry.
8. If Zenodo is connected, confirm the archive and add its version DOI and
   all-versions DOI to `CITATION.cff` and the paper's code-availability section.

## Zenodo setup

Sign in to Zenodo with GitHub, open **Settings -> GitHub**, and enable
`bennettek99-spec/Archaic-DNA-processing-pipeline`. `.zenodo.json` supplies the
software metadata used when a GitHub release is archived.

## Versioning

Use semantic versioning:

- patch: compatible bug or documentation correction;
- minor: new analysis module, report contract, or substantial reproducibility
  improvement;
- major: incompatible estimator, configuration, or output-schema change.

Never tag a release containing partial heavy-run outputs or unreviewed personal
genomic data.
