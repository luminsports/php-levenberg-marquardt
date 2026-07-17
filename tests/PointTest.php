<?php

namespace LuminSports\LevenbergMarquardt\Test;

use LuminSports\LevenbergMarquardt\Point;

it('gets and sets point coordinates', function () {
    $point = new Point(1, 2);

    expect($point->getX())->toEqual(1)
        ->and($point->getY())->toEqual(2);

    $point->setX(10.2);
    $point->setY(5.82E-10);

    expect($point->getX())->toEqual(10.2)
        ->and($point->getY())->toEqual(5.82E-10)
        ->and($point->toArray())->toEqual([
            'x' => 10.2,
            'y' => 5.82E-10,
        ]);
});
