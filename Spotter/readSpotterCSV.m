% readSpotterCSV.m - beta function to read spotter CSVs
% J. Davis
% created on: 2021-10-19
% updates:
%          
% J. Davis, 2021-12-23, cleaned up and stored in Spotter utility folder in
%                       SWIFT-codes directory

function [spotter] = readSpotterCSV(filepath,varargin)
if length(nargin) == 1
    s0 = split(filepath,'/');
    s1 = split(s0{end},'_');
    id = s1{contains(s1,'spot','IgnoreCase',true)};
end

format='%s%f%f%f%f%f%f%f%s%s%s%s%s%s%s%s%s%f%f'; % repmat('%f',1,7)

T1 = readtable(filepath,...
             'Format',format,...
             'ReadVariableNames',true);
ID =repmat(id,height(T1),1);
T0 = table(ID);
spotter= table2struct([T0 T1]);
         
% loop through spotter structure fields:
for field = transpose(fieldnames(spotter))
    % if the first element of a field contains a string, evaluate it:
    if isstring(spotter(1).(field{1})) || ischar(spotter(1).(field{1}))
        
        % if a datetime field, evaluate string as datetime
        if strcmp(field{1},'datetime') 
            C = cellfun(@(x) transpose(datetime(x,'InputFormat','yyyy-MM-dd HH:mm:ssXXX','TimeZone','UTC')),...
                {spotter.(field{1})},'UniformOutput',false);
            [spotter.(field{1})] = C{:};
            
        elseif strcmp(field{1},'ID')
            
        % otherwise, evaluate remaining fields as normal:
        else 
            C = cellfun(@(x) transpose(eval(x)),{spotter.(field{1})},'UniformOutput',false);
            [spotter.(field{1})] = C{:};
        end
    end
end

end 