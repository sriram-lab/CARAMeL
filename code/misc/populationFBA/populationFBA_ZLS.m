function [growthRates, fluxes, EC, OGS] = populationFBA_ZLS(MODEL, nCell, nSample, display, nested, oGS, geneKO)

% THIS COPY WAS ADAPTED FOR INTEGRATION WITH CARAMeL %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   %
% University of Illinois Open Source License                                        %  
% Copyright 2010-2011 Luthey-Schulten Group,                                        %
% All rights reserved.                                                              %
%                                                                                   %
% Developed by: Luthey-Schulten Group                                               %
% 			     University of Illinois at Urbana-Champaign                         %
% 			     http://www.scs.uiuc.edu/~schulten                                  %
%                                                                                   %
% Permission is hereby granted, free of charge, to any person obtaining a copy of   %
% this software and associated documentation files (the Software), to deal with     % 
% the Software without restriction, including without limitation the rights to      % 
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies     % 
% of the Software, and to permit persons to whom the Software is furnished to       % 
% do so, subject to the following conditions:                                       %
%                                                                                   %
% - Redistributions of source code must retain the above copyright notice,          % 
% this list of conditions and the following disclaimers.                            %
%                                                                                   %
% - Redistributions in binary form must reproduce the above copyright notice,       % 
% this list of conditions and the following disclaimers in the documentation        % 
% and/or other materials provided with the distribution.                            %
%                                                                                   %
% - Neither the names of the Luthey-Schulten Group, University of Illinois at       %
% Urbana-Champaign, nor the names of its contributors may be used to endorse or     %
% promote products derived from this Software without specific prior written        %
% permission.                                                                       %
%                                                                                   %
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR        % 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          % 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL          % 
% THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR         % 
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,             %
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR             %
% OTHER DEALINGS WITH THE SOFTWARE.                                                 %
%                                                                                   %
% Author(s): Piyush Labhsetwar, John Cole, Elijah Roberts,                          %
%            Nathan D Price, and Zaida Luthey-Schulten                              %
%                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 1: Load Data Files                      %
load('populationFBA_ZLS_Data.mat', 'vmax_info', 'xie_data');        %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Account for user inputs
%   model
if ~exist('MODEL', 'var') || isempty(MODEL)
    load('populationFBA_ZLS_Data.mat', 'iJO1366_aerobic'); 
    MODEL = iJO1366_aerobic; % default
end
%   cell number (populatio size)
if ~exist('nCell', 'var') || isempty(nCell)
    nCell = 1000; % default
end
%   sample number (enzyme counts)
if ~exist('nSample', 'var') || isempty(nSample)
    nSample = 352; % default 
end
%   display (Boolean)
if ~exist('display', 'var') || isempty(display)
    display = true; % default 
end
%   nested (Boolean)
if ~exist('nested', 'var') || isempty(nested)
    nested = false; % default 
end
%   implement optGpSampler? (Boolean)
if ~exist('oGS', 'var') || isempty(oGS)
    oGS = false; % default 
end
%   simulate geneKO (list of genes)
if ~exist('geneKO', 'var') || isempty(geneKO)
    geneKO = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 2: Setup Output Array                   %
growthRates = zeros(nCell, 1);                                      %
fluxes = zeros(numel(MODEL.rxns), nCell);                           %
EC = zeros(nSample, nCell);
if oGS
    nSamples = 100; 
    OGS = zeros(numel(MODEL.rxns), nSamples, nCell); 
else
    OGS = []; 
end
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Progress bar
if ~nested
    progressbar('Running populationFBA...')
end

for cellNumber = 1:nCell  % We will loop over nCell simulated cells
    while all(fluxes(:, cellNumber) == 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Step 3: For each cell, sample 352 random enzyme counts                 %   %
        enzymeCounts = zeros(352, 1);                                              %  %
        for enzyme = 1:nSample                                                     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enzymeCounts(enzyme) = gamrnd(xie_data.A(enzyme), xie_data.B(enzyme)); %  %                         %
        end                                                                        %   %                        %
        % apply geneKO 
        ix = ismember(xie_data.genes, geneKO); 
        enzymeCounts(ix) = 0; 
        % populate EC output
        EC(:, cellNumber) = enzymeCounts;                                          % 
    %                                                                              %                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %
                                                                                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %
    %       Step 4: For each cell, use the sampled enzyme counts to                                       %     %
    %               impose constraints on the metabolic model                                             %     %
        model = MODEL;                                                                                    %     %
                                                                                                          %     %
        for i=1:length(model.rxns)                                                                        %     %
            switch char(vmax_info(i).cat)                                                                 %     %
                case {'SINGLE', 'SINGLE AND'}                                                             %     %
                    [isinxie, idx] = ismember(vmax_info(i).geneid,xie_data.genes);                        %     %
                    if(isinxie)                                                                           %     %
                        vmax = vmax_info(i).TN*enzymeCounts(idx)*3600/(258e-15*6.023e20);                 %     %
                        model.ub(i) = vmax;                                                               %     %
                    end                                                                                   %     %
                case {'OR'}                                                                               %     %
                    [~, idx] = ismember(vmax_info(i).geneid,xie_data.genes);                              %     %
                    vmax = vmax_info(i).TN*sum(enzymeCounts(idx))*3600/(258e-15*6.023e20);                %     %
                    model.ub(i) = vmax;                                                                   %     %
                case {'MULT AND'}                                                                         %     %
                    [~, idx] = ismember(vmax_info(i).geneid,xie_data.genes);                              %     %
                    vmax = vmax_info(i).TN*min(enzymeCounts(idx))*3600/(258e-15*6.023e20);                %     %
                    model.ub(i) = vmax;                                                                   %     %
                case {'COMPLEX AND'}                                                                      %     %
                    vmax = 0;                                                                             %     %
                    for j=1:vmax_info(i).numgene                                                          %     %
                        eval(sprintf('[~, idx] = ismember(vmax_info(i).geneid%d(:),xie_data.genes);',j))  %     %
                        count = min(enzymeCounts(idx));                                                   %     %
                        vmax = vmax + vmax_info(i).TN*(count *3600)/(258e-15*6.023e20);                   %     %
                    end                                                                                   %     %
                    model.ub(i) = vmax;                                                                   %     %
            end                                                                                           %     %
        end                                                                                               %     %
    %                                                                                                     %     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %
                                                                                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %
    %       Step 5: For each cell, optimize growth rate and store the result       %                            %
    %     soln = optimizeCbModel(model);                                             %                            %
        try
            soln = optimizeCbModel(model, 'max', 'one');
            % fprintf('Solver status: %d\n', soln.stat)
            growthRates(cellNumber) = soln.f;
            fluxes(:, cellNumber) = soln.v;
        catch
            soln = optimizeCbModel(model, 'max', 'one');
            warning('Solver status: %d. Re-optimizing model', soln.stat)
        end
    end
    % add optGpSampler
    if oGS
        if numel(model.rxns) < 500
            nSteps = 50; 
        else
            nSteps = 2500; 
        end
        oGSsoln = optGpSampler(model, [], nSamples, nSteps, 1, 'gurobi', 1);
        OGS(:, :, cellNumber) = oGSsoln.points; 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nested
        progressbar([], cellNumber / nCell)
    else
        progressbar(cellNumber / nCell)
    end

end % Signifies the end of the for loop over the 1,000 cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Step 6: Histogram results and return data                   %
if display                                                          %
    [heights, centers] = hist(growthRates, 40);                     %
    dx = centers(2) - centers(1);                                   %
    figure;                                                         %
    plot(centers, heights / (nCell * dx));                          %
    xlabel('Growth rate (hr^{-1})');                                %
    ylabel('Relative Probability');                                 %
    return;                                                         %
end                                                                 %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %Signifies the end of the function