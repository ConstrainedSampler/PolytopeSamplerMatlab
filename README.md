# PolytopeSampler

PolytopeSampler is a `Matlab` implementation of constrained Hamiltonian Monte Carlo for sampling from high dimensional disributions on polytopes. It is able to sample efficiently from sets and distributions of more than 100K dimensions.

## Quick Tutorial

PolytopeSampler samples from distributions of the form `exp(-f(x))`, for a convex function `f`, subject to constraints `Aineq * x <= bineq`, `Aeq * x = beq` and `lb <= x <= ub`. 

The function `f` can be specified by arrays containing its first and second derivative or function handles. Only the first derivative is required. By default, `f` is empty, which represents a uniform distribution. If the first derivative is a function handle, then the function and its second derivatives must also be provided.

To sample `N` points from a polytope `P`, you can call `sample(P, N)`. The function `sample` will 
1. Find an initial feasible point 
2. Run constrained Hmailton Monte Carlo
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
The polytope has dimension `d = 10` with constraint `sum_i x_i = 1` and `x >= 0`. 
To generate 100 samples uniformly from the polytope `P`, we call the function `sample()`. 
```
>> o = sample(P, 200);
  Time spent |  Time reamin |                  Progress | Samples |  AccRate | StepSize |  MixTime
00d:00:00:01 | 00d:00:00:00 | ######################### | 214/200 | 0.984784 | 0.200000 |     11.0
Done!
```
We can access the samples generated using
```
>> o.samples
```
We can print the the summary of the samples. 
```
>> o.summary

ans =

  10×7 table

                     mean        std         25%         50%         75%      n_ess      r_hat 
                   ________    ________    ________    ________    _______    ______    _______

    samples[1]      0.10023    0.086023    0.034496    0.073419    0.14485    229.98     1.0013
    samples[2]      0.10641    0.085198    0.036705    0.083089    0.15976    242.72     1.0186
    samples[3]     0.096408    0.087197     0.03118     0.07296    0.13157    220.33     1.0153
    samples[4]      0.10016    0.088525    0.036215    0.075896    0.13846    221.46     1.0067
    samples[5]      0.10046    0.087502    0.035919     0.07427    0.14041    226.46    0.99859
    samples[6]      0.09742    0.088638    0.031311    0.074424    0.13731    215.08     1.0042
    samples[7]      0.10408    0.095544    0.034284    0.076434    0.14492    236.23     1.0617
    samples[8]     0.098973    0.080342    0.041192    0.081112     0.1321    226.55     1.0049
    samples[9]      0.10089    0.090578    0.029221     0.07661    0.14485     228.1     1.0023
    samples[10]    0.094967    0.085318    0.030154    0.069197    0.13642    224.57      1.001
```
`n_ess` shows the effective sample size of the samples generated. `r_hat` tests the convergence of the sampling algorithm. A value of `r_hat` close to 1 indicates that the algorithm has converged properly. 

See `demo.m` for more examples, including examples of sampling from non-uniform distributions.  
