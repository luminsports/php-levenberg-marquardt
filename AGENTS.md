# AGENTS.md

Repository-specific guidance for agents working on this standalone PHP library.

## Project Snapshot

- PHP 8.1+ library; Composer package `luminsports/levenberg-marquardt`
- Fits parameterized nonlinear functions with a Levenberg–Marquardt implementation
- Runtime code: `src/`; PHPUnit tests: `tests/`
- Matrix operations: `markbaker/matrix`; PHPStan and PHP CS Fixer for quality checks

## Source of Truth

- Dependencies, autoloading, and scripts: `composer.json`
- Solver options, iteration, error calculation, and prediction: `src/LevenbergMarquardt.php`
- Fitted result contract: `src/Curve.php`; predicted point contract: `src/Point.php`
- Input-length exception: `src/SeriesCountMismatch.php`
- Numerical and validation contracts: `tests/LevenbergMarquardtTest.php`, `tests/PointTest.php`
- PHPUnit, PHPStan, and formatting: `phpunit.xml`, `phpstan.neon.dist`, `.php-cs-fixer.dist.php`

## Commands

```bash
composer install
composer check-style
composer fix-style
composer phpstan
composer test
composer test -- --coverage-text
./vendor/bin/phpunit tests/PointTest.php
./vendor/bin/phpunit --filter test_it_can_guess_the_next_point_on_a_linear_model
```

## Fast Test Runs

- Pest 5's test impact analysis runs only the tests your changes affect: `./vendor/bin/pest --parallel --tia`. The
  first run records a dependency graph from code coverage and takes about as long as a full parallel run. Later runs
  re-run the affected tests and replay cached results for the rest.
- TIA needs PCOV or Xdebug, and every test must be a Pest test: one PHPUnit test class stops the run with "Tia mode
  requires Pest tests". Convert class-based tests first. A `defineEnvironment()` override can move into a trait that
  the test file applies with `uses()`.
- The graph lives in `~/.pest/tia/<project-key>/` (`./vendor/bin/pest --baseline` prints the path) and is shared across
  worktrees. It re-records itself when `composer.lock` or `phpunit.xml` changes; `--tia --fresh` re-records it by hand.
  Edits to comments, docblocks, or formatting re-run nothing.
- `--parallel` alone also cuts a full run several-fold.
- TIA is for local iteration. Do not add `--tia` to CI, and do not count a replayed run as a full-suite run where this
  file asks for one.
- This repository is on Pest 4, which has no TIA. Upgrade to Pest 5 before using `--tia`.

## Package Rules

- Keep the library framework-independent; do not add Laravel-specific runtime behavior.
- Parameter order must remain consistent across the model closure, initial values, bounds, weights, and gradients.
- Reset the cached curve whenever a setter changes solver inputs or options.
- Preserve named fitted parameters for non-variadic model closures and positional values for variadic closures.
- Validate paired series lengths and minimum input size before matrix operations.
- Treat damping, finite-difference, convergence, weighting, and bounds behavior as numerical contracts.
- Use tolerances for floating-point assertions; cover both exact and nonlinear fits.
- Preserve the ml.js attribution and MIT license when changing documentation or packaging.
- Do not weaken PHPStan or reformat unrelated files.

## Handoff

- Run Composer validation, style, PHPStan, focused PHPUnit tests, the full suite, and coverage.
- Report inherited numerical, dependency, or PHP-version warnings separately.
