<?php

/*
 * Copyright (c) 2021 Lumin Sports Technology - All Rights Reserved
 *  Unauthorized copying of this file, via any medium is strictly prohibited
 *  Proprietary and confidential.
 *
 */

namespace LuminSports\LevenbergMarquardt;

use Closure;
use InvalidArgumentException;
use Matrix\Builder;
use Matrix\Functions;
use Matrix\Matrix;

use function Matrix\multiply;

use ReflectionFunction;

/**
 * Adapted from Javascript implementation.
 *
 * @see https://github.com/mljs/levenberg-marquardt
 */
class LevenbergMarquardt
{
    protected float $damping = 1E-2;

    protected int $dampingStepUp = 11;

    protected int $dampingStepDown = 9;

    protected int $maxIterations = 100;

    protected float $errorTolerance = 1E-7;

    protected bool $centralDifference = false;

    protected float|array $gradientDifference = 10E-2;

    protected float|array $weights = 1;

    protected float $improvementThreshold = 1E-3;

    protected array $minValues = [];

    protected array $maxValues = [];

    protected array $xValues = [];

    protected array $yValues = [];

    protected Closure $parameterizedFunction;

    protected array $initialValues = [];

    protected ?Curve $curve = null;

    public function predict(array $xValues): array
    {
        $curve = $this->getCurve();

        $points = [];

        foreach ($xValues as $x) {
            $points[] = new Point($x, $curve->predictY($x));
        }

        return $points;
    }

    public function getCurve(): Curve
    {
        if (isset($this->curve)) {
            return $this->curve;
        }

        $xCount = count($this->getXValues());
        $yCount = count($this->getYValues());

        if ($xCount !== $yCount) {
            throw new SeriesCountMismatch("Number of elements in arrays do not match {$xCount}:{$yCount}");
        }

        $reflection = new ReflectionFunction($function = $this->getParameterizedFunction());

        $initialValues = $this->getInitialValues();

        if (count($initialValues) && count($initialValues) !== $reflection->getNumberOfParameters() && ! $reflection->isVariadic()) {
            throw new InvalidArgumentException('$initialValues must be an array with a length matching the number of arguments for $parameterizedFunction');
        }

        $parameters = count($initialValues) ? $initialValues : array_fill(0, $reflection->getNumberOfParameters(), 1);
        $damping = $this->damping;

        $error = $this->calculateError($parameters);

        $optimalError = $error;
        $optimalParameters = $parameters;

        $converged = $error <= $this->getErrorTolerance();

        $minValues = $this->getMinValues();
        $maxValues = $this->getMaxValues();

        for ($iteration = 0; $iteration < $this->getMaxIterations() && ! $converged; $iteration++) {
            $previousError = $error;

            [$pertubations, $jacobianWeigthResidualError] = $this->step($parameters, $damping);

            foreach ($parameters as $parameterIndex => $parameterValue) {
                $newValue = $parameterValue - $pertubations->getValue($parameterIndex + 1, 1);

                if (isset($minValues[$parameterIndex])) {
                    $newValue = max($minValues[$parameterIndex], $newValue);
                }

                if (isset($maxValues[$parameterIndex])) {
                    $newValue = min($maxValues[$parameterIndex], $newValue);
                }

                $parameters[$parameterIndex] = (float) $newValue;
            }

            $error = $this->calculateError($parameters);

            if (is_nan($error)) {
                break;
            }

            if ($error < $optimalError - $this->getErrorTolerance()) {
                $optimalError = $error;
                $optimalParameters = $parameters;
            }

            $improveMetric = ($previousError - $error) / $pertubations
                    ->transpose()
                    ->multiply(multiply($pertubations, $damping)->add($jacobianWeigthResidualError))
                    ->getValue(1, 1);

            if ($improveMetric > $this->getImprovementThreshold()) {
                $damping = max($damping / $this->getDampingStepDown(), 1e-7);
            } else {
                $damping = min($damping * $this->getDampingStepUp(), 1e7);
            }

            $converged = $error <= $this->getErrorTolerance();
        }

        if (! $reflection->isVariadic()) {
            $keys = array_map(function ($param) {
                return $param->getName();
            }, $reflection->getParameters());

            $optimalParameters = array_combine($keys, $optimalParameters);
        }

        $this->curve = new Curve($optimalParameters, $optimalError, $iteration, call_user_func_array($function, array_values($optimalParameters)));

        return $this->curve;
    }

    public function getDamping(): float
    {
        return $this->damping;
    }

    public function getDampingStepUp(): int
    {
        return $this->dampingStepUp;
    }

    public function getDampingStepDown(): int
    {
        return $this->dampingStepDown;
    }

    public function getMaxIterations(): int
    {
        return $this->maxIterations;
    }

    public function getErrorTolerance(): float
    {
        return $this->errorTolerance;
    }

    public function isCentralDifference(): bool
    {
        return $this->centralDifference;
    }

    public function getGradientDifference(): array|float
    {
        return $this->gradientDifference;
    }

    public function getImprovementThreshold(): float
    {
        return $this->improvementThreshold;
    }

    public function getMinValues(): array
    {
        return $this->minValues;
    }

    public function getMaxValues(): array
    {
        return $this->maxValues;
    }

    public function getXValues(): array
    {
        return $this->xValues;
    }

    public function getYValues(): array
    {
        return $this->yValues;
    }

    public function getParameterizedFunction(): Closure
    {
        return $this->parameterizedFunction;
    }

    public function getInitialValues(): array
    {
        return $this->initialValues;
    }

    public function setInitialValues(array $initialValues): static
    {
        $this->initialValues = array_values($initialValues);
        $this->resetCurve();

        return $this;
    }

    public function setWeights(int|array $weights): static
    {
        $this->weights = $weights;
        $this->resetCurve();

        return $this;
    }

    public function getWeights(): int|array
    {
        return $this->weights;
    }

    public function getWeightSquare(): array
    {
        $weights = $this->getWeights();
        $weightSquare = array_fill(0, count($this->getXValues()), 0);

        foreach ($this->getYValues() as $i => $yValue) {
            $weightSquare[$i] = 1 / pow(is_array($weights) ? $weights[$i] : $weights ?? 1, 2);
        }

        return $weightSquare;
    }

    public function setDamping(float $damping): static
    {
        if ($damping <= 0) {
            throw new InvalidArgumentException('The damping option must be a positive number');
        }

        $this->damping = $damping;
        $this->resetCurve();

        return $this;
    }

    public function setDampingStepUp(int $dampingStepUp): static
    {
        $this->dampingStepUp = $dampingStepUp;
        $this->resetCurve();

        return $this;
    }

    public function setDampingStepDown(int $dampingStepDown): static
    {
        $this->dampingStepDown = $dampingStepDown;
        $this->resetCurve();

        return $this;
    }

    public function setMaxIterations(int $maxIterations): static
    {
        $this->maxIterations = $maxIterations;
        $this->resetCurve();

        return $this;
    }

    public function setErrorTolerance(float $errorTolerance): static
    {
        $this->errorTolerance = $errorTolerance;
        $this->resetCurve();

        return $this;
    }

    public function setCentralDifference(bool $centralDifference): static
    {
        $this->centralDifference = $centralDifference;
        $this->resetCurve();

        return $this;
    }

    public function setGradientDifference(float|array $gradientDifference): static
    {
        $this->gradientDifference = $gradientDifference;
        $this->resetCurve();

        return $this;
    }

    public function setImprovementThreshold(float $improvementThreshold): static
    {
        $this->improvementThreshold = $improvementThreshold;
        $this->resetCurve();

        return $this;
    }

    public function setMinValues(array $minValues): static
    {
        $this->minValues = $minValues;
        $this->resetCurve();

        return $this;
    }

    public function setMaxValues(array $maxValues): static
    {
        $this->maxValues = $maxValues;
        $this->resetCurve();

        return $this;
    }

    public function setXValues(array $xValues): static
    {
        if (count($xValues) < 2) {
            throw new InvalidArgumentException('$xValues must be an array with more than 2 points');
        }

        $this->xValues = $xValues;
        $this->resetCurve();

        return $this;
    }

    public function setYValues(array $yValues): static
    {
        if (count($yValues) < 2) {
            throw new InvalidArgumentException('$yValues must be an array with more than 2 points');
        }

        $this->yValues = $yValues;
        $this->resetCurve();

        return $this;
    }

    public function setParameterizedFunction(Closure $parameterizedFunction): static
    {
        $this->parameterizedFunction = $parameterizedFunction;
        $this->resetCurve();

        return $this;
    }

    public function calculateError(array $parameters)
    {
        $error = 0;
        $func = call_user_func_array($this->getParameterizedFunction(), $parameters);

        $weightSquare = $this->getWeightSquare();

        foreach ($this->getXValues() as $i => $x) {
            $error += pow($this->yValues[$i] - $func($x), 2) / $weightSquare[$i];
        }

        return $error;
    }

    /**
     * @return Matrix[]
     */
    protected function step(array $parameters, float $damping): array
    {
        $identity = Builder::createFilledMatrix(0, count($parameters))->toArray();

        for ($x = 0; $x < count($parameters); $x++) {
            $identity[$x][$x] = $damping;
        }

        $identityMatrix = new Matrix($identity);

        $func = call_user_func_array($this->getParameterizedFunction(), $parameters);

        $evaluatedData = [];

        foreach ($this->getXValues() as $x) {
            $evaluatedData[] = $func($x);
        }

        $gradientMatrix = $this->gradientFunction($parameters, $evaluatedData);
        $residualErrorMatrix = $this->matrixFunction($evaluatedData);

        $inverseMatrix = Functions::inverse(
            $identityMatrix->add(
                $gradientMatrix->multiply(
                    $gradientMatrix->transpose()
                )
            )
        );

        $jacobianWeigthResidualError = $gradientMatrix->multiply($residualErrorMatrix);
        $pertubations = $inverseMatrix->multiply($jacobianWeigthResidualError);

        return [$pertubations, $jacobianWeigthResidualError];
    }

    protected function gradientFunction(array $parameters, array $evaluatedData): Matrix
    {
        $matrix = Builder::createFilledMatrix(0, count($parameters), count($this->getXValues()))->toArray();

        $rowIndex = 0;

        $gradientDifference = is_array($this->getGradientDifference()) ? $this->getGradientDifference() : array_fill(0, count($parameters), $this->getGradientDifference());

        for ($param = 0; $param < count($parameters); $param++) {
            $delta = $gradientDifference[$param] ?? 0;

            if ($delta === 0) {
                continue;
            }

            $auxParams = $parameters;
            $auxParams[$param] += $delta;
            $func = call_user_func_array($this->getParameterizedFunction(), $auxParams);

            if (! $this->isCentralDifference()) {
                foreach ($this->getXValues() as $xIndex => $x) {
                    $matrix[$rowIndex][$xIndex] = ($evaluatedData[$xIndex] - $func($x)) / $delta;
                }
            } else {
                $auxParams = $parameters;
                $auxParams[$param] -= $delta;
                $delta *= 2;
                $func2 = call_user_func_array($this->getParameterizedFunction(), $auxParams);
                foreach ($this->getXValues() as $xIndex => $x) {
                    $matrix[$rowIndex][$xIndex] = ($func2($x) - $func($x)) / $delta;
                }
            }
            $rowIndex++;
        }

        return new Matrix($matrix);
    }

    protected function matrixFunction(array $evaluatedData): Matrix
    {
        $matrix = Builder::createFilledMatrix(null, count($this->getXValues()), 1)->toArray();

        foreach ($this->getYValues() as $yIndex => $y) {
            $matrix[$yIndex][0] = $y - $evaluatedData[$yIndex];
        }

        return new Matrix($matrix);
    }

    protected function resetCurve(): void
    {
        $this->curve = null;
    }
}
