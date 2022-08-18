# Levenberg Marquardt Curve Fitting

[![Code Style](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/php-cs-fixer.yml/badge.svg?branch=main)](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/php-cs-fixer.yml)
[![Tests](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/run-tests.yml/badge.svg?branch=main)](https://github.com/luminsports/php-levenberg-marquardt/actions/workflows/run-tests.yml)

This algorithm is based on the article [Brown, Kenneth M., and J. E. Dennis. "Derivative free analogues of the Levenberg-Marquardt and Gauss algorithms for nonlinear least squares approximation." Numerische Mathematik 18.4 (1971): 289-297.](https://doi.org/10.1007/BF01404679) and [http://people.duke.edu/~hpgavin/ce281/lm.pdf](http://people.duke.edu/~hpgavin/ce281/lm.pdf)

In order to get a general idea of the problem you could also check the [Wikipedia article](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm).

## Credit

This package is a PHP port of [mljs/levenberg-marquardt](https://github.com/mljs/levenberg-marquardt). As such, we have retained their MIT license.

## Installation

```composer require luminsports/levenberg-marquardt```

## Example

```php
use LuminSports\LevenbergMarquardt\LevenbergMarquardt;

// Set up the curve-fitting model
$model = (new LevenbergMarquardt)
    // the parameters and returns a function with the independent variable as a parameter
    ->setParameterizedFunction(function ($criticalPower, $pMax, $tau) {
        return fn (float $t) => ($pMax - $criticalPower) * exp(-$t / $tau) + $criticalPower;
    })
    
    // array of initial parameter values
    ->setInitialValues([
        'criticalPower' => $relative ? 4 : 300,
        'wPrime'        => $relative ? 285 : 20000,
        'tau'           => $relative ? 4 : 300,
    ])
    
    // minimum allowed values for parameters
    ->setMinValues([0, 0, 0])
    
    // maximum allowed values for parameters
    ->setMaxValues([1000, 600, 300])
    
    // Levenberg-Marquardt parameter, small values of the damping parameter λ result in a Gauss-Newton update and large values of λ result in a gradient descent update (default 1E-2)
    ->setDamping(1.5)
    
    // factor to reduce the damping (Levenberg-Marquardt parameter) when there is not an improvement when updating parameters (default: 9)
    ->setDampingStepDown(9)
    
    // factor to increase the damping (Levenberg-Marquardt parameter) when there is an improvement when updating parameters (default: 11)
    ->setDampingStepUp(11)
    
    // the threshold to define an improvement through an update of parameters (default: 1E-3)
    ->setImprovementThreshold(1E-3)
    
    // the step size to approximate the jacobian matrix (default: 10E-2)
    ->setGradientDifference(10E-2)
    
    // if true the jacobian matrix is approximated by central differences otherwise by forward differences (default: false)
    ->setCentralDifference(true)
    
    // maximum of allowed iterations (default: 100)
    ->setMaxIterations(1000)
    
    // minimum uncertainty allowed for each point (default: 10E-3)
    ->setErrorTolerance(10E-3);
    
$curve = $model->setXValues($knownXValues)
    ->setYValues($knownYValues)
    ->getCurve();
    
// outputs the parameters that formed the best fitting curve
$curve->getParameters(); 

// outputs the number of iterations it took to reach this solution
$curve->getIterations(); 

// the sum of the weighted squares of the errors (or weighted residuals) between the known y-coordinates and the curve-fit function.
$curve->getError(); 

// returns an array of points for the provided x-values, with the predicted y-values
$model->predict($xValuesToPredictFor);
```
