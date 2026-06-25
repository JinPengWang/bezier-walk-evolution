# Bézier Walk Evolution (BWE)

A novel metaheuristic optimization algorithm that leverages the geometric properties of Bézier curves to guide the search process in continuous global optimization problems.

> **Paper**: Jinpeng Wang, Xingguo Xu, Yujing Sun, Jiguang Yu, Kaichen Ouyang, and Yuansheng Gao. *Random Walk on Bézier Curves for Global Optimization*.

## Overview

BWE introduces three evolution strategies inspired by the orders of Bézier curves:

| Strategy | Bézier Order | Behavior |
|----------|-------------|----------|
| Cubic Bézier Evolution | 3rd-order | High exploration |
| Quadratic Bézier Evolution | 2nd-order | Balanced |
| Linear Bézier Evolution | 1st-order | High exploitation |

The algorithm adaptively selects among these strategies using Bernstein polynomial-based weights that shift from exploration to exploitation over the course of the optimization.

## Repository Structure

```
├── main.m                  % Entry point — run this file
├── BWE.m                   % Core algorithm implementation
├── GetFunctionsDetails.m   % CEC2017 benchmark function definitions
├── cec17_func.cpp          % CEC2017 C++ source (MEX)
├── cec17_func.mexw64       % Pre-compiled MEX binary (Windows x64)
├── input_data/             % CEC2017 test data (rotation, shift, shuffle matrices)
├── license.txt             % BSD 3-Clause License
└── BWE.pdf                 % Full paper (preprint)
```

## Getting Started

### Prerequisites

- MATLAB R2023b or later
- A C compiler supported by MATLAB (if you need to recompile the MEX file)

### Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/JinPengWang/bezier-walk-evolution.git
   cd bezier-walk-evolution
   ```

2. Open MATLAB and navigate to the project directory.

3. Run the demo:
   ```matlab
   main
   ```

This will optimize CEC2017 benchmark function F1 in 50 dimensions and display the convergence curve.

### Custom Usage

To use BWE on your own objective function, replace the CEC2017 setup in `main.m` with:

```matlab
% Define your objective function (must accept a row vector, return a scalar)
fun = @(x) sum(x.^2);  % Example: Sphere function

% Set problem dimensions and bounds
nvars = 30;             % Number of decision variables
lb = -100 * ones(1, nvars);  % Lower bounds
ub =  100 * ones(1, nvars);  % Upper bounds

% Algorithm parameters
N = 50;   % Population size
T = 500;  % Maximum iterations

% Run BWE
[x, fval, ConvergenceCurve] = BWE(fun, nvars, lb, ub, N, T);
```

## API Reference

### `BWE`

```matlab
[TargetX, TargetF, ConvergenceCurve] = BWE(fun, nvars, lb, ub, N, T)
```

| Parameter | Description |
|-----------|-------------|
| `fun` | Function handle to the objective function |
| `nvars` | Number of decision variables |
| `lb` | Lower bound (scalar or `1 x nvars` vector) |
| `ub` | Upper bound (scalar or `1 x nvars` vector) |
| `N` | Population size |
| `T` | Maximum number of iterations |

| Output | Description |
|--------|-------------|
| `TargetX` | Best solution found |
| `TargetF` | Best objective function value |
| `ConvergenceCurve` | Historical best fitness at each iteration |

## Citation

If you use BWE in your research, please cite:

```bibtex
@article{wang2025bezier,
  title   = {Random Walk on Bézier Curves for Global Optimization},
  author  = {Wang, Jinpeng and Xu, Xingguo and Sun, Yujing and Yu, Jiguang and Ouyang, Kaichen and Gao, Yuansheng},
  year    = {2025}
}
```

## License

This project is licensed under the BSD 3-Clause License — see [license.txt](license.txt) for details.

## Contact

Jinpeng Wang — wangjinpengchunuo@163.com | 2211060117@stu.lntu.edu.cn

Northwestern Polytechnical University
