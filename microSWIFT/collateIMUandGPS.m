function [IMU,GPS] = collateIMUandGPS(IMU,GPS)
%% Crop IMU and GPS

% set crop window based on minimum structure array length:
startCrop = max([IMU.time(1),GPS.time(1)]); %datestr(startCrop, 'YYYY-mm-DD hh:MM:ss.fff'); %
endCrop   = min([IMU.time(end),GPS.time(end)]);

% create logical arrays:
imuCrop = IMU.time >= startCrop & IMU.time <= endCrop;
gpsCrop = GPS.time >= startCrop & GPS.time <= endCrop;

% perform cropping for each field:
IMU = cropStructFields(IMU,imuCrop);
GPS = cropStructFields(GPS,gpsCrop);

%% Interpolate GPS onto master clock
masterTime = IMU.time;

GPS.lat = interp1(GPS.time.',GPS.lat.',masterTime);
GPS.lon = interp1(GPS.time.',GPS.lon.',masterTime);
GPS.sog = interp1(GPS.time.',GPS.sog.',masterTime);
GPS.cog = interp1(GPS.time.',GPS.cog.',masterTime);
GPS.z   = interp1(GPS.time.',GPS.z.'  ,masterTime);
GPS.u   = interp1(GPS.time.',GPS.u.'  ,masterTime);
GPS.v   = interp1(GPS.time.',GPS.v.'  ,masterTime);
GPS.x   = interp1(GPS.time.',GPS.x.'  ,masterTime);
GPS.y   = interp1(GPS.time.',GPS.y.'  ,masterTime);

GPS.time = masterTime;

%% Crop NaN values
nonNaN = ~isnan(GPS.lat) & ~isnan(GPS.lon) & ...
         ~isnan(GPS.sog) & ~isnan(GPS.cog) & ~isnan(GPS.z) & ...
         ~isnan(GPS.u)   & ~isnan(GPS.v)   & ...
         ~isnan(GPS.x)   & ~isnan(GPS.y)   ;

IMU = cropStructFields(IMU,nonNaN);
GPS = cropStructFields(GPS,nonNaN);

end

%% Subfunctions
function x = interp1NaN(x)
    nanx = isnan(x);
    t    = 1:numel(x);
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
end

function S = cropStructFields(S,cropIdx)
    fields = fieldnames(S);
    for f = 1:length(fields)
        if length(S.(fields{f})) == length(cropIdx)
            S.(fields{f}) = S.(fields{f})(cropIdx,:);
        end
    end
end



