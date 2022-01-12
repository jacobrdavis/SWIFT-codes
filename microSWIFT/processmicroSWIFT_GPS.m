function [GPS] = processmicroSWIFT_GPS(GPSfile)
%% processmicroSWIFT_GPS.m
% Function to process microSWIFT GPS data, in the form of a <.dat> file, 
% for a single recording burst.
% 
% Usage:
%   [GPS] = processmicroSWIFT_GPS(GPSfile);
%
% Inputs: 
%   GPSfile - location of a raw GPS data file (in the form of a <.dat> file)
%
% Outputs:
%   GPS - structure with fields:
%       time - UTC time (MATLAB datenum)
%       lat  - GPS-reported latitude (decimal degrees)
%       lon  - GPS-reported longitude (decimal degrees)
%       sog  - GPS-reported speed over ground (m/s)
%       cog  - GPS-reported course over ground (deg)
%       x    - buoy x-position (km); converted from lat-lon
%       y    - buoy y-position (km); converted from lon
%       z    - GPS-reported altitude (m)
%       u    - GPS-reported E-W velocity (m/s)
%       v    - GPS-reported N-S velocity (m/s)
%
% Outputs will be '9999' for invalid results.
%
% Dependencies: 
%   SWIFT-codes (readNMEA)
%   Mapping Toolbox (deg2km)
%
% J. Davis 2021-01-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020
%%

% checks
if exist('GPSfile','var') && isfile(GPSfile)

    % read raw GPS data 
    [lat, lon, sog, cog, depth, time, z] = readNMEA(GPSfile);
    
    % assign to structure as column vector, trimming to the shortest field:
    shortest = min([length(sog),length(cog),length(z)]); % shortest field length
    GPS.lat  = transpose(lat(1:shortest));
    GPS.lon  = transpose(lon(1:shortest));
    GPS.sog  = transpose(sog(1:shortest));
    GPS.cog  = transpose(cog(1:shortest));
    GPS.time = transpose(time(1:shortest));
    GPS.z    = transpose(z(1:shortest));

    % compute velocitites from speed and course over ground
    GPS.u = GPS.sog .* sind(GPS.cog);
    GPS.v = GPS.sog .* cosd(GPS.cog);
    
    % convert GPS signals to km and detrend:
    GPS.x = detrend(deg2km(GPS.lon,cosd(median(GPS.lat))*6371)*1000);
    GPS.y = detrend(deg2km(GPS.lat)*1000);
    
    % approximate sample rate:
    GPS.samplingrate = length(GPS.time)./((max(GPS.time)-min(GPS.time))*24*3600); % Hz

else % if conditions are not satisfied, report an empty structure (also useful for initialization)
    
    GPS = struct('lat',[],'lon',[],'sog',[],'cog',[],'time',[],'z',[],'u',[],'v',[],'x',[],'y',[],'samplingrate',[]);

end



