% querySofarAPI.m - beta function to query Sofar spotter API
% J. Davis
% created on: 2021-10-18
% example usage:
%
%   First create SofarAPIopts structure as input:
%         SofarAPIopts = struct();
%         SofarAPIopts.token                      ='your_token_here';
%         SofarAPIopts.limit                      =100;
%         SofarAPIopts.spotterId                  ='SPOT-1177';
%         SofarAPIopts.startDate                  ='2021-05-04';
%         SofarAPIopts.endDate                    ='2021-05-08';
%         SofarAPIopts.includeWaves               ='true';
%         SofarAPIopts.incudeWind                 ='true';
%         SofarAPIopts.includeSurfaceTempData     ='true';
%         SofarAPIopts.includeTrack               ='true';
%         SofarAPIopts.includeFrequencyData       ='true';
%         SofarAPIopts.includeDirectionalMoments  ='true';
%
%   Call querySofarAPI.m:
%         [spotter1] = querySofarAPI(SofarAPIopts);
%
% updates:
%          
% J. Davis, 2021-12-23, cleaned up and stored in Spotter utility folder in
%                       SWIFT-codes directory


function [out] = querySofarAPI(SofarAPIopts,varargin)

if length(nargin) == 1
    concatcols=true;
end

baseurl='https://api.sofarocean.com/api/wave-data?'; % see: https://docs.sofarocean.com/spotter-and-smart-mooring/spotter-data/wave-data

options = weboptions("ContentType", "json"); % see: https://www.mathworks.com/help/matlab/ref/webread.html#bue6uid-options

data=webread(baseurl,                                                   ...
    'token',SofarAPIopts.token,                                         ...
    'limit',SofarAPIopts.limit,                                         ...
    'spotterId',SofarAPIopts.spotterId,                                 ...
    'startDate',SofarAPIopts.startDate,                                 ...
    'endDate',SofarAPIopts.endDate,                                     ...
    'includeWaves',SofarAPIopts.includeWaves,                           ...
    'incudeWind',SofarAPIopts.incudeWind,                               ...
    'includeSurfaceTempData',SofarAPIopts.includeSurfaceTempData,       ...
    'includeTrack',SofarAPIopts.includeTrack,                           ...
    'includeFrequencyData',SofarAPIopts.includeFrequencyData,           ...
    'includeDirectionalMoments',SofarAPIopts.includeDirectionalMoments, ...
    options);

    data = data.data;
    
    % Concatenate columns into a single structure (optional, but recommended
    % and required for conversion into a SWIFT-compliant structure)

    if concatcols == true
        if strcmp(data.waves(1).timestamp, data.frequencyData(1).timestamp)==true
            ID = repmat(data.spotterId,length(data.waves),1);
            T0 = table(ID);
            T1 = struct2table(rmfield(data.waves,{'timestamp','latitude','longitude'}));       
            T2 = struct2table(data.frequencyData);
            T3 = [T0 T1 T2];
            T = movevars(T3,'timestamp','After',1);
            T.Properties.VariableNames{2} = 'datetime';
            
            C = cellfun(@(x) transpose(datetime(x,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSX','TimeZone','UTC')),...
                T.datetime,'UniformOutput',false);
            [T.datetime] = C;
            
            out =  table2struct(T);
        end
    else % Output not concatenated into structure
       out = data;
    end
end

    
    
    
    
    