Quick readme: 

Start with demo_template.m, uniform_sampler.m or gaussian_sampler.m and modify the problem to the one you want to sample; each file has one or more examples. 

In more detail:

The main component of this package is the function 
sample in sample.m
to sample according to a density proportional to a given function of the form 

f(x) = exp(-sum f_i(x_i)) 

restricted to a polyhedron defined by
{Aineq x <= bineq, Aeq x = beq, lb <= x <= ub}

The function f is specified by a vector function of its 1st, 2nd and 3rd derivatives.
Only the first derivative is required. For uniform sampling this is empty ([]).
If 2nd derivative is provided, 3rd derivative must also be provided, else it is assumed to be zero.
If the first derivative is a function handle, then the function and its second and third derivatives must also be provided.

This core function sample.m is supplemented by functions to: 
1. find an initial feasible point 
2. test convergence of the sampling algorithm with Effective Sample Size, and a uniformity test (for uniform sampling).

Before using sample, we set up the parameters for sampling via a struct called opts with the following properties:

			   seed: 'shuffle'
                        nSketch: JL dimension for fast leverage score computation (0 means no JL).
               adaptiveStepSize: True/False
                weightedBarrier: True/False
                   dynamicBound: True/False
                      odeMethod: default: @implicit_midpoint
                     outputFunc: [function_handle]
                      debugMode: True/False
                        maxTime: 86400 (max sampling time in seconds)
                        maxStep: 300000 (maximum number of steps)
                    minStepSize: 1.0000e-08  (minimum step size of the ODE solver)

 
The output of sample is a struct "o" with fields including:

                samples: dim x #steps
	    prepareTime: pre-processing time 
                    dim: dimension of instance
               nSamples: Number of nearly independent samples
                    ess: Effective Sample Size in each coordinate
             mixingTime: estimated #steps to mix
             sampleTime: total sampling time


