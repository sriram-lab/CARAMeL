function model = change_media(model, filename, queried_medium)
% @author: Scott Campit
% Updated: Carolina H. Chung 
%% media.m defines the medium constraints we will impose on the Genome-scale metabolic model. 
% By default, the substrate uptake rates were set to RPMI conditions by default. 
    % Other medium conditions were scaled w.r.t RPMI amounts (using ratios
    % from concentrations as the scaling factor).
    
file = which(filename);
if verLessThan('matlab', '9.6.0.1072779')
    [~, sheetNames] = xlsfinfo(file);
    for sheets = 1:length(sheetNames)
        if ismember(sheetNames(sheets), queried_medium)
            [adjustedLB, rxn_ids] = xlsread(file, sheetNames{sheets});
            rxn_ids(1,:) = [];
            rxn_ids(:,1) = [];
            
            for rxn=1:length(rxn_ids)
                model.lb(ismember(model.rxns, rxn_ids(rxn, 1))) = ...
                    adjustedLB(rxn, 4);
            end
            
        elseif ismember({'nan'}, queried_medium)
            [adjustedLB, rxn_ids] = xlsread(file, 'RPMI');
            rxn_ids(1,:) = [];
            rxn_ids(:,1) = [];
            
            for rxn=1:length(rxn_ids)
                model.lb(ismember(model.rxns, rxn_ids(rxn, 1))) = ...
                    adjustedLB(rxn, 4);
            end
            
        end
    end
    
else
    [~, sheetNames] = xlsfinfo(file);
    for sheets = 1:length(sheetNames)
        if ismember(sheetNames(sheets), queried_medium)
            dataArray = readcell(file,...
                'Sheet', sheetNames{sheets});
            dataArray(1,:) = [];
%             dataArray(:,1) = [];
            
            for rxn=1:length(dataArray)
                model.lb(strcmpi(model.rxns, ...
                    dataArray{rxn,3})) = cell2mat(dataArray(rxn, 10));
            end
            
        elseif ismember({'nan'}, queried_medium)
            dataArray = readcell(file,...
                'Sheet', 'GEM Default');
            dataArray(1,:) = [];
%             dataArray(:,1) = [];
            
            for rxn=1:length(dataArray)
                model.lb(strcmpi(model.rxns, ...
                    dataArray{rxn,3})) = cell2mat(dataArray(rxn, 10));
            end
        end
    end
end

end