<?php

namespace LuminSports\LevenbergMarquardt\Test;

use InvalidArgumentException;
use LuminSports\LevenbergMarquardt\LevenbergMarquardt;
use LuminSports\LevenbergMarquardt\SeriesCountMismatch;
use PHPUnit\Framework\TestCase;

/**
 * These tests (and this entire package) are basically a carbon-copy of mljs/levenberd-marquadt.
 *
 * https://github.com/mljs/levenberg-marquardt/blob/c8317bc28753976e1182cf9b381080a361d56aaa/src/__tests__/curve.test.js
 */
class LevenbergMarquardtTest extends TestCase
{
    public function contrivedProblems(): array
    {
        return [
            [
                'name'              => 'bennet5([2, 3, 5])',
                'fn'                => fn ($b1, $b2, $b3) => fn (float $t) => $b1 * pow($t + $b2, -1 / $b3),
                'n'                 => 154,
                'xStart'            => -2.6581,
                'xEnd'              => 49.6526,
                'problemParameters' => [2, 3, 5],
                'options'           => [
                    'damping'       => 0.00001,
                    'maxValues'     => [11, 11, 11],
                    'minValues'     => [1, 2.7, 1],
                    'initialValues' => [3.5, 3.8, 4],
                ],
            ],
            [
                'name'              => '4*sin(2*t)',
                'fn'                => fn ($a, $b) => fn (float $x) => $a * sin($b * $x),
                'n'                 => 20,
                'xStart'            => 0,
                'xEnd'              => 19,
                'problemParameters' => [4, 2],
                'options'           => [
                    'damping'         => 0.1,
                    'dampingStepUp'   => 1,
                    'dampingStepDown' => 1,
                    'initialValues'   => [5.8, 2.4],
                ],
            ],
            [
                'name'                       => 'Sigmoid',
                'fn'                         => fn ($a, $b, $c) => fn (float $t) => $a / ($b + exp(-$t * $c)),
                'n'                          => 20,
                'xStart'                     => 0,
                'xEnd'                       => 19,
                'problemParameters'          => [2, 2, 2],
                'options'                    => [
                    'damping'       => 0.1,
                    'initialValues' => [3, 3, 3],
                ],
                'decimalsForErrorValues'     => 2,
                'decimalsForParameterValues' => 1,
            ],
            [
                'name'                       => 'Sum of lorentzians',
                'fn'                         => fn (...$p) => function (float $t) use ($p) {
                    $result = 0;

                    for ($i = 0; $i < count($p); $i += 3) {
                        $p2 = pow($p[$i + 2] / 2, 2);
                        $factor = $p[$i + 1] * $p2;
                        $result += $factor / (pow($t - $p[$i], 2) + $p2);
                    }

                    return $result;
                },
                'n'                          => 100,
                'xStart'                     => 0,
                'xEnd'                       => 99,
                'problemParameters'          => [1.05, 0.1, 0.3, 4, 0.15, 0.3],
                'options'                    => [
                    'damping'            => 0.1,
                    'gradientDifference' => [0.01, 0.0001, 0.0001, 0.01, 0.0001, 0],
                    'initialValues'      => [1.1, 0.15, 0.29, 4.05, 0.17, 0.3],
                ],
                'decimalsForErrorValues'     => 2,
                'decimalsForParameterValues' => 1,
            ],
            [
                'name'                       => 'Sum of lorentzians, central differences',
                'fn'                         => fn (...$p) => function (float $t) use ($p) {
                    $result = 0;

                    for ($i = 0; $i < count($p); $i += 3) {
                        $p2 = pow($p[$i + 2] / 2, 2);
                        $factor = $p[$i + 1] * $p2;
                        $result += $factor / (pow($t - $p[$i], 2) + $p2);
                    }

                    return $result;
                },
                'n'                          => 100,
                'xStart'                     => 0,
                'xEnd'                       => 99,
                'problemParameters'          => [1, 0.1, 0.3, 4, 0.15, 0.3],
                'options'                    => [
                    'damping'            => 0.1,
                    'gradientDifference' => [0.01, 0.0001, 0.0001, 0.01, 0.0001],
                    'centralDifference'  => true,
                    'initialValues'      => [1.1, 0.15, 0.29, 4.05, 0.17, 0.28],
                    'errorTolerance'     => 10e-8,
                ],
                'decimalsForErrorValues'     => 2,
                'decimalsForParameterValues' => 1,
            ],
        ];
    }

    public function realWorldProblems()
    {
        return [
            [
                'name'     => 'fourParamEq',
                'fn'       => fn ($a, $b, $c, $d) => fn ($t) => $a + ($b - $a) / (1 + pow($c, $d) * pow($t, -$d)),
                'x'        => [9.22e-12, 5.53e-11, 3.32e-10, 1.99e-9, 1.19e-8, 7.17e-8, 4.3e-7, 0.00000258, 0.0000155, 0.0000929],
                'y'        => [7.807, -3.74, 21.119, 2.382, 4.269, 41.57, 73.401, 98.535, 97.059, 92.147],
                'expected' => [
                    'iterations'      => 200,
                    'parameterError'  => 16398.0009709,
                    'parameterValues' => [-16.7697, 43.4549, 1018.8938, -4.3514],
                ],
                'options'  => [
                    'damping'       => 0.00001,
                    'maxIterations' => 200,
                    'weights'       => 1,
                    'initialValues' => [0, 100, 1, 0.1],
                ],
            ],
        ];
    }

    /**
     * @dataProvider contrivedProblems
     *
     * @group        contrived
     */
    public function test_it_fits_each_contrived_problem($name, $curveFunction, $n, $xStart, $xEnd, $problemParameters, $options, $decimalsForErrorValues = 2, $decimalsForParameterValues = 3)
    {
        $xCoords = array_map(function ($i) use ($xStart, $xEnd, $n) {
            return $xStart + ($i * ($xEnd - $xStart)) / ($n - 1);
        }, range(0, $n));

        $yCoords = array_map($curveFunction(...$problemParameters), $xCoords);

        $curve = (new LevenbergMarquardt)->setParameterizedFunction($curveFunction)
            ->setDamping($options['damping'] ?? 1E-2)
            ->setDampingStepDown($options['dampingStepDown'] ?? 9)
            ->setDampingStepUp($options['dampingStepUp'] ?? 11)
            ->setMaxIterations($options['maxIterations'] ?? 100)
            ->setGradientDifference($options['gradientDifference'] ?? 10E-2)
            ->setCentralDifference($options['centralDifference'] ?? false)
            ->setErrorTolerance($options['errorTolerance'] ?? 1E-7)
            ->setMaxValues($options['maxValues'] ?? [])
            ->setMinValues($options['minValues'] ?? [])
            ->setInitialValues($options['initialValues'] ?? [])
            ->setXValues($xCoords)
            ->setYValues($yCoords)
            ->getCurve();

        $this->assertEqualsWithDelta($problemParameters, array_values($curve->getParameters()), 1 / pow(10, $decimalsForParameterValues ?? 3));
        $this->assertEqualsWithDelta(0, $curve->getError(), 1 / pow(10, $decimalsForErrorValues ?? 2));
    }

    /**
     * @dataProvider realWorldProblems
     *
     * @group        realWorld
     */
    public function test_it_fits_each_real_world_problem($name, $curveFunction, $xCoords, $yCoords, $expected, $options)
    {
        $curve = (new LevenbergMarquardt)->setParameterizedFunction($curveFunction)
            ->setDamping($options['damping'] ?? 1E-2)
            ->setDampingStepDown($options['dampingStepDown'] ?? 9)
            ->setDampingStepUp($options['dampingStepUp'] ?? 11)
            ->setMaxIterations($options['maxIterations'] ?? 100)
            ->setGradientDifference($options['gradientDifference'] ?? 10E-2)
            ->setCentralDifference($options['centralDifference'] ?? false)
            ->setErrorTolerance($options['errorTolerance'] ?? 1E-7)
            ->setMaxValues($options['maxValues'] ?? [])
            ->setMinValues($options['minValues'] ?? [])
            ->setInitialValues($options['initialValues'] ?? [])
            ->setXValues($xCoords)
            ->setYValues($yCoords)
            ->getCurve();

        $this->assertEqualsWithDelta($expected['parameterError'], $curve->getError(), 0.001);
        $this->assertEqualsWithDelta($expected['parameterValues'], array_values($curve->getParameters()), 0.001);
        $this->assertEquals($expected['iterations'], $curve->getIterations());
    }

    /**
     * @group other
     */
    public function test_it_should_return_solution_with_lowest_error()
    {
        $xCoords = [
            0, 0.6283185307179586, 1.2566370614359172, 1.8849555921538759,
            2.5132741228718345, 3.141592653589793, 3.7699111843077517,
            4.39822971502571, 5.026548245743669, 5.654866776461628,
        ];
        $yCoords = [
            0, 1.902113032590307, 1.1755705045849465, -1.175570504584946,
            -1.9021130325903073, -4.898587196589413e-16, 1.902113032590307,
            1.1755705045849467, -1.1755705045849456, -1.9021130325903075,
        ];

        $sinFunction = fn ($a, $b) => fn (float $t) => $a * sin($b * $t);

        $model = (new LevenbergMarquardt)
            ->setParameterizedFunction($sinFunction)
            ->setDamping(1.5)
            ->setMaxIterations(100)
            ->setGradientDifference(1E-2)
            ->setErrorTolerance(1E-2)
            ->setInitialValues([0.594398586701882, 0.3506424963635226])
            ->setXValues($xCoords)
            ->setYValues($yCoords);

        $curve = $model->getCurve();

        $manualCalculatedError = array_sum(array_map(function ($xCoord, $yCoord) use ($curve, $sinFunction) {
            return pow($yCoord - $sinFunction(...$curve->getParameters())($xCoord), 2);
        }, $xCoords, $yCoords));

        $this->assertEqualsWithDelta($curve->getError(), $manualCalculatedError, $model->getErrorTolerance());
        $this->assertEqualsWithDelta($curve->getError(), 15.52, $model->getErrorTolerance());
    }

    /**
     * @group other
     */
    public function test_it_can_get_model_settings()
    {
        $nullFunction = fn () => null;

        $model = (new LevenbergMarquardt)
            ->setParameterizedFunction($nullFunction)
            ->setWeights([1, 2, 3])
            ->setMinValues([4, 5, 6])
            ->setMaxValues([7, 8, 9])
            ->setInitialValues([5, 6, 7])
            ->setImprovementThreshold(1E-2)
            ->setXValues([1, 2])
            ->setYValues([3, 4]);

        $this->assertEquals($model->getDamping(), 1E-2);
        $this->assertEquals($model->getDampingStepUp(), 11);
        $this->assertEquals($model->getDampingStepDown(), 9);
        $this->assertEquals($model->getMaxIterations(), 100);
        $this->assertEquals($model->getErrorTolerance(), 1E-7);
        $this->assertEquals($model->getGradientDifference(), 10E-2);
        $this->assertEquals($model->getImprovementThreshold(), 1E-2);
        $this->assertEquals($model->getMinValues(), [4, 5, 6]);
        $this->assertEquals($model->getMaxValues(), [7, 8, 9]);
        $this->assertEquals($model->getInitialValues(), [5, 6, 7]);
        $this->assertEquals($model->getWeights(), [1, 2, 3]);
        $this->assertEquals($model->getXValues(), [1, 2]);
        $this->assertEquals($model->getYValues(), [3, 4]);
    }

    /**
     * @group other
     */
    public function test_error_should_be_zero_for_an_exact_fit()
    {
        $model = $this->createLinearModel(1, 1, 10);

        $this->assertEquals($model->calculateError([1, 1]), 0, 0.001);
    }

    /**
     * @group other
     */
    public function test_it_can_guess_the_next_point_on_a_linear_model()
    {
        $model = $this->createLinearModel(1, 1, 10);

        $curve = $model->getCurve();
        $points = array_map(fn ($p) => $p->toArray(), $model->predict([11, 12]));

        $this->assertEquals($curve->getIterations(), 0);
        $this->assertEquals($curve->getError(), 0);
        $this->assertEquals($curve->getParameters(), ['slope' => 1, 'intercept' => 1]);
        $this->assertEquals($points, [
            [
                'x' => 11,
                'y' => 12,
            ],
            [
                'x' => 12,
                'y' => 13,
            ],
        ]);
    }

    /**
     * @group other
     */
    public function test_error_should_match__the_sum_of_absolute_difference_between_the_model_and_the_data()
    {
        $model = $this->createLinearModel(1, 1, 10);

        $this->assertEquals($model->calculateError([1, 1]), 0, 0.001);
    }

    /**
     * @group other
     */
    public function test_curve_is_not_recalculated_when_parameters_are_unchanged()
    {
        $model = $this->createLinearModel(1, 1, 10);

        $this->assertSame($model->getCurve(), $model->getCurve());
    }

    /**
     * @group other
     */
    public function test_curve_is_recalculated_when_parameters_change()
    {
        $model = $this->createLinearModel(1, 1, 10);
        $curve = $model->getCurve();

        $this->assertEquals($curve->getParameters(), ['slope' => 1, 'intercept' => 1]);

        $model->setYValues(array_map(fn ($x) => 2 * $x + 1, range(0, 9)));
        $curve = $model->getCurve();

        $this->assertEqualsWithDelta($curve->getParameters(), ['slope' => 2, 'intercept' => 1], 0.1);
    }

    /**
     * @group exceptions
     */
    public function test_it_should_throw_an_error_when_negative_damping_provided()
    {
        $this->expectException(InvalidArgumentException::class);

        (new LevenbergMarquardt())->setDamping(-1);
    }

    /**
     * @group exceptions
     */
    public function test_it_should_throw_an_error_when_too_few_x_coords_are_set()
    {
        $this->expectException(InvalidArgumentException::class);

        (new LevenbergMarquardt())
            ->setXValues([1])
            ->getCurve();
    }

    /**
     * @group exceptions
     */
    public function test_it_should_throw_an_error_when_too_few_y_coords_are_set()
    {
        $this->expectException(InvalidArgumentException::class);

        (new LevenbergMarquardt())
            ->setYValues([1])
            ->getCurve();
    }

    /**
     * @group exceptions
     */
    public function test_it_should_throw_an_error_when_x_coords_and_y_coords_are_not_the_same_length()
    {
        $this->expectException(SeriesCountMismatch::class);

        (new LevenbergMarquardt())
            ->setXValues([1, 2, 3])
            ->setYValues([1, 2, 3, 4])
            ->getCurve();
    }

    /**
     * @group exceptions
     */
    public function test_it_should_throw_an_error_when_initial_values_do_not_match_function_parameter_count()
    {
        $this->expectException(InvalidArgumentException::class);

        (new LevenbergMarquardt())
            ->setParameterizedFunction(fn ($a, $b) => fn ($t) => $a * $b * $t)
            ->setInitialValues([1, 2, 3])
            ->setXValues([1, 2, 3])
            ->setYValues([1, 2, 3])
            ->getCurve();
    }

    /**
     * @group exceptions
     */
    public function test_it_should_return_initial_values_if_function_evaluates_to_nan_after_starting()
    {
        // Note: This test is identical to the other fourParamEq test, except for the
        // damping parameter. The increased damping option leads to a case where
        // c < 0 && d is not an integer so Math.pow(c, d) is NaN
        $model = (new LevenbergMarquardt())
            ->setParameterizedFunction(fn ($a, $b, $c, $d) => fn ($t) => $a + ($b - $a) / (1 + pow($c, $d) * pow($t, -$d)))
            ->setDamping(0.01)
            ->setMaxIterations(200)
            ->setInitialValues([0, 100, 1, 0.1])
            ->setXValues([9.22e-12, 5.53e-11, 3.32e-10, 1.99e-9, 1.19e-8, 7.17e-8, 4.3e-7, 0.00000258, 0.0000155, 0.0000929])
            ->setYValues([7.807, -3.74, 21.119, 2.382, 4.269, 41.57, 73.401, 98.535, 97.059, 92.147]);

        $curve = $model->getCurve();

        $this->assertEquals($curve->getIterations(), 0);
        $this->assertEqualsWithDelta($curve->getError(), 19289.706, 0.001);
        $this->assertEqualsWithDelta(array_values($curve->getParameters()), [0, 100, 1, 0.1], $model->getErrorTolerance());
    }

    protected function createLinearModel(float $slope, float $intercept, int $numberOfPoints): LevenbergMarquardt
    {
        $linearFunction = fn ($slope, $intercept) => fn ($x) => $slope * $x + $intercept;

        $xCoords = range(0, $numberOfPoints - 1);
        $yCoords = array_map($linearFunction($slope, $intercept), $xCoords);

        return (new LevenbergMarquardt)
            ->setParameterizedFunction($linearFunction)
            ->setWeights(1)
            ->setXValues($xCoords)
            ->setYValues($yCoords);
    }
}
