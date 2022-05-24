function [fluxstate, grate, solverobj] = constrain_flux_regulation(...
    model1, onreactions, offreactions, ...
    kappa, rho, epsilon, mode, epsilon2, minfluxflag)

    % DESCRIPTION 
    % This function simulates flux through reactions annotated in a
    % genome-scale metabolic model (GEM) by constraining the model using
    % information on "on" and "off" reactions (can be based directly or by 
    % differentially expressed gene (DEG) data. 
    % 
    % STEPS 
    % 1. Inspect inputs and set default values
    % 2. Convert model to gurobi format
    % 3. Maximize flux through "on" reactions
    % 4. Minimize flux through "off" reactions
    % 5. Define function outputs
    % 
    % Author:   Sriram Chandrasekaran
    % Created:  October 12, 2010
    % Updates:  September 20, 2020 (Carolina H. Chung)
    %           November 6, 2019 (Scott Campit)
    %           November 8, 2018 (Fangzhou Shen)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. model1:          genome-scale metabolic model (GEM)
           -> type:         MATLAB structure
    
        2. onreactions:     list of "on" reactions (when 'mode' = 1) or 
                            list of up-regulated genes (when 'mode' = 0)
           -> type:         cell array 
           -> dim:          m x 1 (m = # of "on" reactions/up-reg genes)
    
        3. offreactions:    list of  "off" reactions (when 'mode' = 1) or 
                            list of down-regulated genes (when 'mode' = 0)
           -> type:         cell array
           -> dim:          m x 1 (m = # of "off" reactions/down-reg genes)
    
    OPTIONAL INPUTS:
        1. kappa:       relative weight for "on" reactions/up-reg genes
           -> default:  kappa = 1
           -> note:     use small value (e.g. 1e-3) for large networks
    
        2. rho:         relative weight for "off" reactions/down-reg genes
           -> default:  rho = 1
           -> note:     use small value (e.g. 1e-3) for large networks
    
        3. epsilon:     minimum flux for "on" reactions
           -> default:  epsilon = 1e-3
           
        4. mode:        specifies whether inputs are reaction or gene lists 
                        (mode = 0 for genes     mode = 1 for reactions)
           -> default:  mode = 1
    
        5. epsilon2:    minimum flux for "off" reactions 
           -> default:  epsilon2 = 0
    
        6. minfluxflag: Boolean to minimize flux sum through all reactions
                        (i.e. apply PFBA)
           -> default:  minfluxflag = true
    
    OUTPUTS:
        1. fluxstate:       flux solution through reactions in GEM
        2. grate:           simulated growth rate
        3. solverobj:       flux through GEM objective function
    
    EXAMPLE USAGE: 
        1. Default usage: 
           >> [fluxstate, grate, solverobj] = constrain_flux_regulation(...
                    model1, onreaction_list, offreaction_list)
    
        2. Use gene lists instead with defaults: 
           >> [fluxstate, grate, solverobj] = constrain_flux_regulation(...
                    model1, up_genes, down_genes, [], [], [], 0)
    %}

%% INSPECT INPUTS AND SET DEFAULT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % mode (reactions or genes)
    if (~exist('mode','var')) || (isempty(mode))
        mode = 1;
    end
    if mode == 0 
        [~, ~, onreactions, ~] =  deleteModelGenes(model1, onreactions);
        [~, ~, offreactions, ~] =  deleteModelGenes(model1, offreactions);
    end

    % Function parameters
    if (~exist('epsilon','var')) || (isempty(epsilon))
        epsilon = ones(size(onreactions)) * 1E-3;
    end
    if numel(epsilon) == 1
        epsilon = repmat(epsilon, size(onreactions));
    end
    if (~exist('rho','var')) || (isempty(rho))
        rho = ones(size(onreactions));
    end
    if numel(rho) == 1
        rho  = repmat(rho, size(onreactions));
    end
    if (~exist('kappa','var')) || (isempty(kappa))
        kappa = ones(size(offreactions));
    end
    if numel(kappa) == 1
        kappa  = repmat(kappa, size(offreactions));
    end
    if (~exist('epsilon2','var')) || (isempty(epsilon2))
        epsilon2 = zeros(size(offreactions));
    end
    if numel(epsilon2) == 1
        epsilon2  = repmat(epsilon2, size(offreactions));
    end
    if (~exist('minfluxflag','var')) || (isempty(minfluxflag))
        minfluxflag = true; 
    else
        minfluxflag = logical(minfluxflag); 
    end

    % Don't display output after running FBA
    params.outputflag = 0;

    % Apply PFBA (if prompted)
    if minfluxflag
        kappa = [kappa(:); ...
            ones(size(setdiff(model1.rxns, offreactions))) * 1E-6]; 
        epsilon2 = [epsilon2; ...
            zeros(size(setdiff(model1.rxns, offreactions)))];
        offreactions = [offreactions(:); ...
            setdiff(model1.rxns, offreactions)];
    end

%% CONVERT TO GUROBI FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Re-define model
    model = model1;
    model.A = model1.S;
    model.obj = model1.c;
    model.rhs = model1.b;
    if exist('model1.csense','var') && ~isempty(model1.csense)
        model.sense = model1.csense;
        model.sense(ismember(model.sense, 'E')) = '=';
        model.sense(ismember(model.sense, 'L')) = '<';
        model.sense(ismember(model.sense, 'G')) = '>';
    else
        model.sense = repmat('=', [size(model1.S, 1), 1]);
    end
    model.lb = model1.lb;
    model.ub = model1.ub;
    model.vtype = repmat('C', size(model1.S, 2), 1);
    model.modelsense = 'max';
%     nrows = size(model.A, 1);
%     ncols = size(model.A, 2);
    M = 10000;
    objpos = logical(model1.c);
    nrxns = length(model1.rxns);

%% MAXIMIZE FLUX THROUGH "ON" REACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:length(onreactions)
        rxnpos = find(ismember(model1.rxns, onreactions(j)));

        %         xi - (eps + M)ti >= -M
        % ti = 0 or 1.

        rowpos = size(model.A, 1) + 1;
        colpos = size(model.A, 2) + 1;

        model.A(rowpos, rxnpos) = 1;
        model.A(rowpos, colpos) = -(1*epsilon(j) + M);
        model.rhs(rowpos) = -M;
        model.sense(rowpos) = '>';
        model.vtype(colpos) = 'B';
        model.obj(colpos) = 1*rho(j);
        model.lb(colpos) = 0;
        model.ub(colpos) = 1;

        % xi + (eps + M)ri <= M
        % ri = 0 or 1.

        rowpos = size(model.A,1) + 1;
        colpos = size(model.A,2) + 1;

        model.A(rowpos,rxnpos) = 1;
        model.A(rowpos,colpos) = (1*epsilon(j) + M);
        model.rhs(rowpos) = M;
        model.sense(rowpos) = '<';
        model.vtype(colpos) = 'B';
        model.obj(colpos) = 1*rho(j);
        model.lb(colpos) = 0;
        model.ub(colpos) = 1;

    end

%% MINIMIZE FLUX THROUGH "OFF" REACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for jj = 1:length(offreactions)
        rxnpos = find(ismember(model1.rxns,offreactions(jj)));
        %        xi + si >= -eps2
        %     si >= 0
        %     rho(ri + si)
        % constraint 1
        rowpos = size(model.A,1) + 1;
        colpos = size(model.A,2) + 1;
        model.A(rowpos,rxnpos) = 1;
        model.A(rowpos,colpos) = 1;
        model.rhs(rowpos) = -epsilon2(jj);
        model.sense(rowpos) = '>';
        model.vtype(colpos) = 'C';
        % set si to be positive
        model.lb(colpos) = 0;
        model.ub(colpos) = 1000;
        model.obj(colpos) = -1*kappa(jj); % minimized

        % constraint 2
        %     xi - ri <= eps2
        %     ri >= 0
        % new row and column
        rowpos = size(model.A,1) + 1;
        colpos = size(model.A,2) + 1;
        model.A(rowpos,rxnpos) = 1;
        model.A(rowpos,colpos) = -1;
        model.rhs(rowpos) = epsilon2(jj);
        model.sense(rowpos) = '<';
        model.vtype(colpos) = 'C';
        % set ri to be positive
        model.lb(colpos) = 0;
        model.ub(colpos) = 1000;
        model.obj(colpos) = -1*kappa(jj); % minimized
    end

%% DEFINE FUNCTION OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solg1 = gurobi(model, params);
    try
        fluxstate = solg1.x(1:nrxns);
        grate = solg1.x(objpos);
        solverobj = solg1.objval;

    catch
        warning("Could not determine FBA solution. Returning NaN values.")
        fluxstate = nan(nrxns, 1);
        grate = NaN;
        solverobj = NaN;
    end

end

