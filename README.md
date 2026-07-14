# Levenberg–Marquardt Curve Fitting for PHP

[![Code Style](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/php-cs-fixer.yml/badge.svg?branch=main)](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/php-cs-fixer.yml)
[![PHPStan](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/phpstan.yml/badge.svg?branch=main)](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/phpstan.yml)
[![Tests](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/tests.yml)

Fit nonlinear parameterized functions to paired numeric data with a fluent, dependency-light PHP API. The package supports bounds, weights, forward or central finite differences, damping controls, error and iteration reporting, and predictions from the fitted curve.

## Requirements

- PHP 8.1 or later
- `markbaker/matrix` 2.1 or later

## Installation

```bash
composer require luminsports/levenberg-marquardt
```

## Example

```php
use LuminSports\LevenbergMarquardt\LevenbergMarquardt;

$xValues = [0, 1, 2, 3, 4];
$yValues = [1, 3, 5, 7, 9];

$model = (new LevenbergMarquardt)
    ->setParameterizedFunction(
        fn (float $slope, float $intercept) =>
            fn (float $x): float => $slope * $x + $intercept,
    )
    ->setInitialValues([1, 0])
    ->setXValues($xValues)
    ->setYValues($yValues);

$curve = $model->getCurve();

$curve->getParameters(); // approximately ['slope' => 2.0, 'intercept' => 1.0]
$curve->getError();
$curve->getIterations();
$model->predict([5, 6]); // array of Point objects
```

The parameterized function receives the values to fit and must return a function of the independent variable. Initial, minimum, and maximum values follow the same parameter order. `predict()` returns `Point` objects; call `toArray()` when an `['x' => ..., 'y' => ...]` representation is needed.

## Development

```bash
composer install
composer check-style
composer phpstan
composer test
```

Generate a text coverage report with:

```bash
composer test -- --coverage-text
```

## Background and Credit

The implementation is based on the derivative-free nonlinear least-squares methods described by [Brown and Dennis](https://doi.org/10.1007/BF01404679) and is a PHP port of [mljs/levenberg-marquardt](https://github.com/mljs/levenberg-marquardt). The upstream MIT license and attribution are retained in [LICENSE.md](LICENSE.md).

## License

Released under the MIT License. See [LICENSE.md](LICENSE.md) and [COPYRIGHT](COPYRIGHT).
