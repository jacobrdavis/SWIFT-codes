%% valuation_plots.m
% Script to plot evaluation metrics between microSWIFT and CDIP (or 
% similar) wave spectral estimates
% 
% J. Davis, 2022-03-07
%
% Dependencies: 
%   Mapping Toolbox (deg2km)

dataDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/matfiles/';
evalflist = dir([dataDir,'*.mat']); 


ids = unique(cellfun(@(name) extract(name,"microSWIFT"+digitsPattern(3)), {evalflist.name}));

% n_SWIFTS = length(ids);
% n_metrics = 6;

% T = array2table(zeros(6,n_metrics));
% T.Properties.RowNames = {'RMSE_total','RMSE_swell','RMSE_sea','bias_total','bias_swell','bias_sea'};
% T.Properties.VariableNames = {'energy','a1','b1','a2','b2','check'};

rowNames = {'RMSE_total','RMSE_swell','RMSE_sea','bias_total','bias_swell','bias_sea'};
variableNames = {'energy','a1','b1','a2','b2','check'};


for i = 1:length(evalflist)
    disp(evalflist(i).name)
    load(fullfile(evalflist(i).folder,evalflist(i).name));
    id = SWIFT(1).id{1};
    method = SWIFT(1).metrics.method;
    referenceFrame = SWIFT(1).metrics.referenceFrame;
    filterType = SWIFT(1).metrics.filterType;
    config = [method,'_',referenceFrame,'_',filterType];

    for j=1:length(SWIFT)
        tmp(:,:,j) = table2array(SWIFT(j).metrics.spectral(:,:));
    end

    mean_tmp = mean(tmp,3); 
    std_tmp  = std(tmp,0,3);
    
    results.(config).(id).mean = mean_tmp;
    results.(config).(id).std = std_tmp;
end

%%
configs = transpose(fieldnames(results));

for f = 1:length(configs)
    for m = 1:length(ids)
        tmp2(:,:,m) = results.(configs{f}).(ids{m}).mean;   
    end

    mean_tmp2 = array2table(mean(tmp2,3));
    mean_tmp2.Properties.VariableNames = variableNames;
    mean_tmp2.Properties.RowNames = rowNames;

    std_tmp2 = array2table(std(tmp2,0,3));
    std_tmp2.Properties.VariableNames = variableNames;
    std_tmp2.Properties.RowNames = rowNames;

    results.(configs{f}).combined.mean = mean_tmp2;
    results.(configs{f}).combined.std = std_tmp2;

end

%%


for v = 1:length(variableNames)
    variable = variableNames{v};
    figure(v); hold on
    title(variable)
    for f = 1:length(configs)
        
        b =  bar(f,...
            [results.(configs{f}).combined.mean.(variable)('RMSE_total'); ...
            results.(configs{f}).combined.mean.(variable)('RMSE_sea');...
            results.(configs{f}).combined.mean.(variable)('RMSE_swell')]);
    
        b(1).FaceColor = [.2 .6 .2];
        b(2).FaceColor = [.6 .2 .2];
        b(3).FaceColor = [.2 .2 .6];
    
    end
    xticks(1:length(configs))
    xticklabels(strrep(configs,'_',' '))
    xtickangle(gca,90)
    ylabel('RMSE (m^2/Hz)')

% set(gca, 'YScale', 'log')
end

%%

%         for metric = T.Properties.VariableNames 
%             metric = metric{1};
%             tmp(j) = SWIFT(j).metrics.spectral.energy(metric);
%             T.(metric)(id) = 
%         end
            
%         RMSE_total(j) = SWIFT(j).metrics.spectral.energy('RMSE_total');
%         RMSE_swell(j) = SWIFT(j).metrics.spectral.energy('RMSE_swell');
%         RMSE_sea(j)   = SWIFT(j).metrics.spectral.energy('RMSE_sea');
%         bias_total(j) = SWIFT(j).metrics.spectral.energy('RMSE_total');
%         bias_swell(j) = SWIFT(j).metrics.spectral.energy('RMSE_total');
%         bias_sea(j)   = SWIFT(j).metrics.spectral.energy('RMSE_total');
    

