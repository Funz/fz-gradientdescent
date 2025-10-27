# fz-gradientdescent

Gradient Descent optimization algorithm plugin for fz (Funz parametric computing framework).

## Description

This repository contains a port of the classic GradientDescent algorithm from the old Funz plugin system to the new fz format. The Gradient Descent algorithm is a first-order local optimization method that uses finite differences to estimate gradients and iteratively move toward local optima.

**Algorithm reference:** [Gradient Descent on Wikipedia](http://en.wikipedia.org/wiki/Gradient_descent)

## Installation

This algorithm can be used with the fz framework by placing the `.fz/algorithms/gradientdescent.R` file in your project's `.fz/algorithms/` directory, or by referencing it directly.

## Usage

### Basic Example (Minimization)

```python
import fz

results = fz.fzd(
    input_file="input.txt",
    input_variables={"x1": "[0;1]", "x2": "[0;1]"},
    model="mymodel",
    output_expression="result",
    algorithm=".fz/algorithms/gradientdescent.R",
    calculators=["sh://bash calc.sh"],
    algorithm_options={
        "yminimization": True,
        "max_iterations": 50,
        "ytol": 0.01,
        "delta": 0.1,
        "epsilon": 0.01
    }
)
```

### Basic Example (Maximization)

```python
results = fz.fzd(
    input_file="input.txt",
    input_variables={"temperature": "[0;100]", "pressure": "[1;10]"},
    model="mymodel",
    output_expression="efficiency",
    algorithm=".fz/algorithms/gradientdescent.R",
    calculators=["sh://bash calc.sh"],
    algorithm_options={
        "yminimization": False,  # Maximize
        "max_iterations": 100,
        "ytol": 0.001,
        "target": 0.95  # Stop early if efficiency > 0.95
    }
)
```

## Algorithm Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `yminimization` | boolean | `true` | Minimize output value? Set to `false` for maximization |
| `ytol` | numeric | `0.1` | Convergence precision on output value |
| `max_iterations` | integer | `100` | Maximum number of iterations |
| `delta` | numeric | `1` | Gradient step factor (initial value, auto-adjusted) |
| `epsilon` | numeric | `0.01` | Relative finite difference step for gradient estimation |
| `target` | numeric | `Inf` | Output target limit for early convergence (use `-Inf` for maximization) |
| `x0` | string | `""` | Starting input values (comma-separated), e.g., "0.5,0.5" |

## Algorithm Behavior

1. **Initialization**: Starts from a random point (or provided `x0`) and evaluates function values at the point plus finite differences to estimate the gradient.

2. **Iteration**: At each iteration:
   - Estimates the gradient using finite differences
   - Takes a step in the direction of steepest descent (minimization) or ascent (maximization)
   - Auto-adjusts step size (`delta`) if the objective function doesn't improve
   - Keeps points within bounds using reflection

3. **Convergence**: Stops when:
   - Maximum iterations reached
   - Change in objective value is below `ytol`
   - Target value is reached (if specified)

## Output

The algorithm provides:

- **Final analysis**: 
  - Optimum value found
  - Location of optimum (input values)
  - Number of iterations and evaluations
  - Visualization plot (pairs plot for multi-dimensional, scatter plot for 1D)

- **Intermediate progress**:
  - Current iteration number
  - Current best value
  - Number of evaluations so far

## Technical Details

### Algorithm Type
- **Type**: Local optimization
- **Order**: First-order (uses gradient information)
- **Method**: Steepest descent/ascent with finite difference gradients

### Dependencies
- Base R (no special packages required)
- Optional: `base64enc` for HTML visualization output

### Porting Notes

This algorithm has been ported from the original Funz plugin format to the new fz format:
- Original: [algorithm-GradientDescent](https://github.com/Funz/algorithm-GradientDescent)
- New format reference: [fz examples/algorithms/montecarlo_uniform.R](https://github.com/Funz/fz/blob/implement-algorithms/examples/algorithms/montecarlo_uniform.R)

Key changes:
- Constructor pattern using S3 classes
- Methods renamed: `getInitialDesign` → `get_initial_design`, `getNextDesign` → `get_next_design`, etc.
- Input/output format adapted to new fz expectations (list of points instead of matrices)
- Return empty `list()` to signal completion instead of `NULL`
- State management using environments for mutable state

## License

BSD 3-Clause License (same as original Funz project)

## Author

Yann Richet (ported to new fz format)