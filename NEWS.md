# Changes in IPY

This page describes the most important changes in `IPY`. The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic
Versioning](https://semver.org/spec).

## Unreleased

### Added

- New hyperbolic loss function.

- New Huber loss function.

- New L1-L2 white regularization.

### Removed

- Dependency on `yfitsio` has been suppressed.

- Roughness regularization has been removed as it is redundant with other
  regularizations.

### Fixed

- Define `M_PI` if not in `<math.h>`.
