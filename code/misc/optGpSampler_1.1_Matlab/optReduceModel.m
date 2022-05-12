function [reducedModel, fixedRxnIds, fixedFluxValues, removedMetIDs,removedRxnIDs] = optReduceModel(model, solverName, stoichiometryTolerance, fluxTol)
% optReduceModel; reduces the constraint-based model by removing all reactions and associated metabolites that can not carry flux.
% 
% optReduceModel is very similar to the COBRA Toolbox function "reduceModel". Main differences are that optReduceModel allows
% to omit very small stoichiometric coefficients for a more stable model. This is motivated by the fact that in some models, 
% the objective value can not be  reached after propagating the flux lower- and upper bounds. By removing these small 
% coefficients (usually removing coefficients < 1e-4 from the objective function), the model can also produce biomass 
% after flux propagation (FVA without objective function).
%
% Author: WL Megchelenbrink, Radboud University Nijmegen, Nijmegen, the Netherlands
% Version 1.1
% Date: 30/04/2014
% 
% INPUTS:
% - model: a metabolic model in COBRA format. Mandatory fields are {S, lb, ub, c, b, mets, rxns}
% - solverName: the LP to use for model reduction; choose from {'cplex', 'gurobi', 'glpk'}
% - stoichiometryTolerance: remove stoichiometric coefficients with smaller absolute values than 'stoichiometryTolerance'. Default = 0 (no removal)
% - fluxTol: reactions with ub-lb < fluxTol are considered fixed. Moreover, if abs(ub) < tol, the reaction 
%   is considered "dead" and will be removed from the model. Default = 1e-6.
%
% OUTPUTS:
% - reducedModel: the model with all "dead" reactions and orphanized metabolites removed
% - fixedRxnIds: reactionIds that have a fixed flux (lb=ub)
% - fixedFluxValues: flux of the fixed reaction ids
% - removedMetIDs: row indices of the metabolites removed in the pre-processing step (by setting stoichiometryTol > 0)
% - removedRxnIDs: column indices of the reactions removed in the pre-processing step (by setting stoichiometryTol > 0)
%
% EXAMPLE:
%  sModel = optReduceModel(model, 'cplex',0)
%  Reduces the model by removing reactions that have ub-lb < 1e-6. The stoichiometric matrix is not changed prior to reduction.

verbose = 1;

% Set small stoichiometric coefficients (i.e. S_i,j < stoichiometricTolerance) to zero
if nargin < 3 || stoichiometryTolerance == 0
    fixStoichiometry = 0;
else
    fixStoichiometry = 1;
end

% Cut-off tolerance for reactions that can/ can not carry flux
if nargin < 4 || isempty(fluxTol)
   fluxTol = 1e-6;           
end

removedMetIDs = [];
removedRxnIDs = [];

if fixStoichiometry
    fbaOrgMin = optimizeCbModel(model, 'min');
    fbaOrgMax = optimizeCbModel(model, 'max');

    % Set small absolute values in the stoichiometry matrix to zero,
    % to improve the quality of the warmup points
    [removedMetIDs,removedRxnIDs] = find(model.S ~=0 & abs(model.S) < stoichiometryTolerance);
    nElems = length(unique(find(model.S ~=0 & abs(model.S) < stoichiometryTolerance)));
    model.S(removedMetIDs,removedRxnIDs) = 0;

    fbaChangedMin = optimizeCbModel(model, 'min');
    fbaChangedMax = optimizeCbModel(model, 'max');

    if verbose
        fprintf('---------- STOICHIOMETRIC MATRIX CHANGED -------------\n');
        fprintf('%d stoichiometric coefficients with abs(S_{i,j}) < %4.4f have been set to zero\n', nElems, stoichiometryTolerance);
        fprintf('Minimal FBA original model = %4.4f\n', fbaOrgMin.f);
        fprintf('Maximum FBA original model = %4.4f\n', fbaOrgMax.f);
        fprintf('-----------------------------------\n');
        fprintf('Minimal FBA changed model = %4.4f\n', fbaChangedMin.f);
        fprintf('Maximum FBA changed model = %4.4f\n', fbaChangedMax.f);
        fprintf('------------------------------------------------------\n');
    end
end

[reducedModel, fixedRxnIds, fixedFluxValues] = cReduceModel(model, solverName, fluxTol);
end