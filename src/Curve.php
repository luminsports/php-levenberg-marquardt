<?php

/*
 * Copyright (c) 2021 Lumin Sports Technology - All Rights Reserved
 *  Unauthorized copying of this file, via any medium is strictly prohibited
 *  Proprietary and confidential.
 *
 */

namespace LuminSports\LevenbergMarquardt;

class Curve
{
    public function __construct(protected array $parameters, protected float $error, protected int $iterations, protected \Closure $curveFunction)
    {
    }

    public function getParameters(): array
    {
        return $this->parameters;
    }

    public function getError(): float
    {
        return $this->error;
    }

    public function getIterations(): int
    {
        return $this->iterations;
    }

    public function predictY(float $x): float
    {
        return call_user_func($this->curveFunction, $x);
    }
}
