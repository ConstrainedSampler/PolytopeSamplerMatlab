Quick readme: 

See demo.m on how to generate random samples from polytope

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
2. test convergence of the sampling algorithm with Effective Sample Size.

You can set up the parameters for sampling via a struct called opts. Here are some parameters you may want to set:

                     outputFunc: [function_handle]
                        maxTime: 86400 (max sampling time in seconds)
                        maxStep: 300000 (maximum number of steps)

 
The output of sample is a struct "o" with fields including:

                samples: dim x #steps
            prepareTime: pre-processing time
               nSamples: Number of nearly independent samples
                    ess: Effective Sample Size in each coordinate
             sampleTime: total sampling time


