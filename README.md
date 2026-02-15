# fz-gradientdescent

A [Funz](https://github.com/Funz/fz) algorithm plugin implementing **Gradient Descent** first-order local optimization in R.

This repository provides the Gradient Descent algorithm for the fz framework, ported from the original [algorithm-GradientDescent](https://github.com/Funz/algorithm-GradientDescent) Funz plugin.

**Algorithm reference:** [Gradient Descent on Wikipedia](http://en.wikipedia.org/wiki/Gradient_descent)

## Features

### Algorithm Interface (R S3 class)

The algorithm implements the fz R algorithm interface:

- `GradientDescent(...)`: S3 constructor accepting algorithm-specific options
- `get_initial_design.GradientDescent(obj, input_variables, output_variables)`: Return initial point + finite differences
- `get_next_design.GradientDescent(obj, X, Y)`: Compute gradient, take a step, return new point + finite differences, or `list()` when converged
- `get_analysis.GradientDescent(obj, X, Y)`: Return optimum value, location, and optional visualization
- `get_analysis_tmp.GradientDescent(obj, X, Y)`: Return intermediate progress (current iteration, best value)

### Algorithm Behavior

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

## Requirements

- **R** must be installed on your system
- **rpy2** Python package: `pip install rpy2`
- **fz** framework: `pip install git+https://github.com/Funz/fz.git`

## Installation

```bash
pip install git+https://github.com/Funz/fz.git
pip install rpy2
```

### Install the Algorithm Plugin

```python
import fz
fz.install_algorithm("gradientdescent")
```

Or from a URL:
```python
fz.install_algorithm("https://github.com/Funz/fz-gradientdescent")
```

Or using the CLI:
```bash
fz install gradientdescent
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
| `x0` | string | `""` | Starting input values (comma-separated), e.g., `"0.5,0.5"` |

## Usage

### Without fzd (standalone algorithm testing)

You can test the algorithm without any simulation code, using rpy2 directly:

```python
from rpy2 import robjects

# Source the R algorithm
robjects.r.source(".fz/algorithms/gradientdescent.R")
r_globals = robjects.globalenv

# Create an instance
r_algo = robjects.r["GradientDescent"](
    yminimization=True, max_iterations=20, ytol=0.01,
    delta=0.1, epsilon=0.01
)

# Define input variable ranges as R list
r_input_vars = robjects.r('list(x1 = c(0.0, 1.0), x2 = c(0.0, 1.0))')
r_output_vars = robjects.StrVector(["output"])

# Get initial design
r_design = r_globals['get_initial_design'](r_algo, r_input_vars, r_output_vars)
print(f"Initial design: {len(r_design)} points")
```

Or via fz's automatic wrapper:

```python
from fz.algorithms import load_algorithm

# Load R algorithm (fz handles rpy2 wrapping automatically)
algo = load_algorithm("gradientdescent",
                      yminimization=True, max_iterations=20, ytol=0.01)

# Same Python interface as Python algorithms
design = algo.get_initial_design(
    {"x1": (0.0, 1.0), "x2": (0.0, 1.0)}, ["output"]
)
print(f"Initial design: {len(design)} points")
```

### With fzd (coupled with a model)

Use `fz.fzd()` to run the algorithm coupled with a model and calculators:

```python
import fz

# Install model and algorithm plugins
fz.install("Model")
fz.install_algorithm("gradientdescent")

# Run optimization
analysis = fz.fzd(
    input_path="examples/Model/input.txt",
    input_variables={"x": "[0;10]", "y": "[-5;5]"},
    model="Model",
    output_expression="result",
    algorithm="gradientdescent",
    algorithm_options={
        "yminimization": True,
        "max_iterations": 50,
        "ytol": 0.01,
        "delta": 0.1,
        "epsilon": 0.01
    },
    calculators="localhost_Model",
    analysis_dir="analysis_results_gd"
)

print(analysis)
```

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

## Directory Structure

```
fz-gradientdescent/
├── .fz/
│   └── algorithms/
│       └── gradientdescent.R         # R algorithm implementation (S3 class)
├── .github/
│   └── workflows/
│       └── test.yml                  # CI workflow (includes R setup)
├── tests/
│   └── test_plugin.py                # Test suite (uses rpy2)
├── examples/
│   └── test_gradientdescent.R        # R-only test script
├── example_standalone.ipynb          # Notebook: algorithm without fzd
├── example_with_fzd.ipynb            # Notebook: algorithm with fzd
├── LICENSE
└── README.md
```

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

Key changes from old format:
- Constructor pattern using S3 classes
- Methods renamed: `getInitialDesign` → `get_initial_design`, `getNextDesign` → `get_next_design`, etc.
- Input/output format adapted to new fz expectations (list of points instead of matrices)
- Return empty `list()` to signal completion instead of `NULL`
- State management using environments for mutable state

## Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Or directly
python tests/test_plugin.py
```

## License

BSD 3-Clause License (same as original Funz project)

## Author

Yann Richet (ported to new fz format)
