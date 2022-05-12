function sModel = optGpSampler(model, warmupPts, nSamples, nSteps, nThreads, solverName, verbose)
% optGpSampler; a fast sampling tool for sampling (genome-scale) constraint-based metabolic networks
%
% Author: WL Megchelenbrink, Radboud University Nijmegen, Nijmegen, the Netherlands
% Version 1.1
% Date: 30/04/2014
% 
% INPUTS:
% - model: a metabolic model in COBRA format. Mandatory fields are {S, lb, ub, c, b}
% - warmupPts: previously collected warmupPts. These can be obtained from sModel.warmupPts
% - nSamples: nr of samples to collect
% - nThreads: nr of parallel threads to use 
% - nSteps: nr of steps between samples (for larger models, you need to increase the step size to obtain better samples)
% - solverName: the LP to use for the warmup points; choose from {'cplex', 'gurobi', 'glpk'}
% - verbose: output logging to the screen, 0=no output, 1=output
%
% OUTPUTS:
% - sModel: the model with attached sample and warmup points
%
% EXAMPLE:
%  sModel = optGpSampler(model, [], 1e5, 500, 4, 'cplex', 1)
%  samples 500.000 sample points using 500 steps between each point. Uses the cplex LP solver and 4 parallel threads.


if nargin < 7
    verbose = 1;
end

% Scale very large fluxes bounds to a maximum of -1000, 1000 for better
% robustness
scaleFactor = max( max(abs(model.lb)), max(abs(model.ub))) / 1000;
if scaleFactor > 1
    model.lb = model.lb / scaleFactor;
    model.ub = model.ub / scaleFactor;
end

[warmups, samples] = cOptGpSampler(model, warmupPts,  nSamples, nSteps, nThreads, solverName, verbose);

% Undo the scaling effect
if scaleFactor > 1
    warmups = warmups * scaleFactor;
    samples = samples * scaleFactor;
end

sModel = model;
sModel.warmups = warmups;
sModel.points = samples;


end