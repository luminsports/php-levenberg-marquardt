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
