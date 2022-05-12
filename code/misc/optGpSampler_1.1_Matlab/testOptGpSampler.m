function sModel = testOptGpSampler(solverName)
    
    if nargin < 1
        solverName = 'glpk';
    end
   
    load Ecoli_MFA;
    nSamples = 1e3;
    nSteps = 50;
    nCores = 1;
    verbose = 1;
    
    sModel = optGpSampler(model, [], nSamples, nSteps, nCores, solverName, verbose);
end