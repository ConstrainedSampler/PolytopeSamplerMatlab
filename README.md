# PolytopeSampler

PolytopeSampler is a `Matlab` implementation of constrained Riemannian Hamiltonian Monte Carlo for sampling from high dimensional disributions on polytopes. It is able to sample efficiently from sets and distributions with more than 100K dimensions.

## Quick Tutorial

PolytopeSampler samples from distributions of the form `exp(-f(x))`, for a convex function `f`, subject to constraints `Aineq * x <= bineq`, `Aeq * x = beq` and `lb <= x <= ub`. 

The function `f` can be specified by arrays containing its first and second derivative or function handles. Only the first derivative is required. By default, `f` is empty, which represents a uniform distribution. If the first derivative is a function handle, then the function and its second derivatives must also be provided.

To sample `N` points from a polytope `P`, you can call `sample(P, N)`. The function `sample` will 
1. Find an initial feasible point 
2. Run constrained Hamiltonian Monte Carlo
3. Test convergence of the sampling algorithm by computing Effective Sample Size (ESS) and terminate when `ESS >= N`. If the target distribution is uniform, a uniformity test will also be performed.

Extra parameters can be set up using `opts`. Some useful parameters include `maxTime` and `maxStep`. By default, they are set to 
```
                        maxTime: 86400 (max sampling time in seconds)
                        maxStep: 300000 (maximum number of steps)
```
The output is a struct `o`, which stores samples generated in `o.samples` and a summary of the sample in `o.summary`. `o.samples` is an array of size `dim x #steps`.

                
### Example

We demonstrate PolytopeSampler using a simple example, sampling uniformly from a simplex.
The polytope is defined by

```
>> P = struct;
>> d = 10;
>> P.Aeq = ones(1, d);
>> P.beq = 1;
>> P.lb = zeros(d, 1);
```
The polytope has dimension `d = 10` with constraint `sum_i x_i = 1` and `x >= 0`. This is a simplex.
To generate 200 samples uniformly from the polytope `P`, we call the function `sample()`. 
```
>> o = sample(P, 200);
  Time spent |  Time reamin |                  Progress | Samples |  AccProb | StepSize |  MixTime
00d:00:00:01 | 00d:00:00:00 | ######################### | 211/200 | 0.989903 | 0.200000 |     11.2
Done!
```
We can access the samples generated using
```
>> o.samples
```
We can print a summary of the samples: 
```
>> o.summary

ans =

  10×7 table

                     mean        std         25%         50%         75%      n_ess      r_hat 
                   ________    ________    ________    ________    _______    ______    _______

    samples[1]     0.093187    0.091207    0.026222    0.064326    0.13375    221.51    0.99954
    samples[2]     0.092815    0.086905    0.027018    0.066017    0.13221    234.59     1.0301
    samples[3]      0.10034    0.090834    0.030968    0.075631    0.13788    216.56     1.0159
    samples[4]      0.10531    0.092285    0.035363    0.077519     0.1481    235.25     1.0062
    samples[5]      0.10437    0.087634    0.034946    0.080095     0.1533    212.54    0.99841
    samples[6]       0.1029    0.093724    0.028774    0.074354    0.15135     227.6     1.0052
    samples[7]       0.1042    0.083084    0.038431    0.081964    0.15352    231.54     1.0008
    samples[8]     0.088778    0.086902    0.025565    0.062473    0.11837    229.69     1.0469
    samples[9]      0.10627     0.09074    0.036962    0.084294    0.15125    211.64    0.99856
    samples[10]     0.10184    0.084699    0.035981    0.074923    0.14578    230.63     1.0277
```
`n_ess` shows the effective sample size of the samples generated. 
`r_hat` tests the convergence of the sampling algorithm. 
A value of `r_hat` close to 1 indicates that the algorithm has converged properly. 

See `demo.m` for more examples, including examples of sampling from non-uniform distributions.  
