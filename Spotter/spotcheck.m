% spotcheck.m - script to check and demonstrate spotter utility codes
% J. Davis
% created on: 2021-10-18
% updates:
%          
% J. Davis, 2021-12-23, cleaned up and stored in Spotter utility folder in
%                       SWIFT-codes directory

%% Obtaining spotter data - method 1: query SofarAPI

% create SofarAPIopts structure as input to uerySofarAPI.m
SofarAPIopts = struct();
SofarAPIopts.token                      ='206de34bbc07fce0507cf9791aacc2';
SofarAPIopts.limit                      =100;
SofarAPIopts.spotterId                  ='SPOT-1177';
SofarAPIopts.startDate                  ='2021-05-04';
SofarAPIopts.endDate                    ='2021-05-08';
SofarAPIopts.includeWaves               ='true';
SofarAPIopts.incudeWind                 ='true';
SofarAPIopts.includeSurfaceTempData     ='true';
SofarAPIopts.includeTrack               ='true';
SofarAPIopts.includeFrequencyData       ='true';
SofarAPIopts.includeDirectionalMoments  ='true';

% call querySofarAPI.m
[spotter1] = querySofarAPI(SofarAPIopts);

%% Obtaining spotter data - method 2: load Spotter data previously stored as a CSV
file = '../data/SPOT-1177_2021-05-04_to_2021-10-17.csv';
[spotter2] = readSpotterCSV(file);

%% Convert to "SWIFT-compliant" structure for use with SWIFT-codes
[spotter] = spotter2SWIFTcompliant(spotter2);

%% Test directional moments match
directiontest(spotter2,1)
 
%% Test use with SWIFTdirectionalspectra
pre_hurricane = spotter([spotter.time]>datenum(2021,8,20,18,0,0) & [spotter.time]<datenum(2021,8,20,19,0,0));
in_hurricane = spotter([spotter.time]>datenum(2021,8,21,14,0,0) & [spotter.time]<datenum(2021,8,21,15,0,0));
post_hurricane = spotter([spotter.time]>datenum(2021,8,22,12,0,0) & [spotter.time]<datenum(2021,8,22,13,0,0));

SWIFTdirectionalspectra(pre_hurricane)
SWIFTdirectionalspectra(in_hurricane)
SWIFTdirectionalspectra(post_hurricane)