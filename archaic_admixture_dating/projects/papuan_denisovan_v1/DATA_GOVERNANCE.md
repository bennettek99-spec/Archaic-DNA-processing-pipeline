# Data governance

Papuan and Papua New Guinean genomic data are not ordinary interchangeable
software inputs. Technical availability does not imply ethical permission.

## Required checks

Before analysis, record:

- the official data source and exact release;
- consent and permitted secondary-use scope;
- whether Indigenous governance or community-specific conditions apply;
- controlled-access approval and expiration where applicable;
- restrictions on population labels, linkage, redistribution, and publication;
- whether tract-level derivatives may be retained or shared;
- the approved people and systems that may access the files.

## Operational rules

- Never bypass authentication, access committees, application procedures, or
  download restrictions.
- Never commit raw genomes, sample-level restricted metadata, credentials,
  signed URLs, or controlled-access derivatives.
- Use approved sample identifiers and respectful population names; retain a
  private mapping only when explicitly permitted.
- Store inputs outside Git and reference them through a machine-local manifest.
- Configure downloads as `controlled` or `manual` unless automated public reuse
  is clearly permitted. `--force` cannot override this classification.
- Publish only aggregate results allowed by the relevant terms.
- Remove or archive data according to access agreements, not merely project
  convenience.

## Indigenous data sovereignty

Analyses involving Indigenous genomes should consider governance frameworks
such as CARE in addition to technical FAIR principles. Researchers must seek
appropriate community and institutional guidance rather than treating this
document as ethics approval.

## Why raw genomes are excluded from Git

Git is designed for broad replication and permanent history. It is unsuitable
for large files, revocable access, controlled redistribution, sensitive
metadata, and participant/community governance. The repository therefore
contains only synthetic fixtures, code, aggregate machine-readable outputs
that are permitted for sharing, and provenance that does not expose protected
paths or credentials.
