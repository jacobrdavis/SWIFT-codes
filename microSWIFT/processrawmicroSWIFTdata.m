%% processrawmicroSWIFTdata.m
% Script to process microSWIFT raw data from GPS and IMU; 
% 
% J. Davis, 2021-01-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020

%% Settings and inputs:
dataDir   = '/Users/jacob/Dropbox/Projects/microSWIFT/OffshoreExample2/RawData/microSWIFT002/';
outputDir = '/Users/jacob/Dropbox/Projects/microSWIFT/OffshoreExample2/';

%% GPS processing

% Obtain GPS file list from specified data directory:
GPSflist = dir([dataDir,'*GPS*.dat']); 

% Initialize GPS structure array:
clear GPS; GPS(length(GPSflist)) = processmicroSWIFT_GPS(); 

% Process bursts:
for gi = 1:length(GPSflist)
    disp(['GPS file ' num2str(gi) ' of ' num2str(length(GPSflist))])
    GPS(gi) = processmicroSWIFT_GPS(fullfile(GPSflist(gi).folder,GPSflist(gi).name));
end

%% IMU processing

% Obtain IMU file list from specified data directory:
IMUflist = dir([dataDir,'*IMU*.dat']);

% Initialize IMU structure array:
clear IMU; IMU(length(IMUflist)) = processmicroSWIFT_IMU(); 

% Process bursts:
for ii = 1:length(IMUflist)
    disp(['IMU file ' num2str(ii) ' of ' num2str(length(IMUflist))])
    IMU(ii) = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),'body');
end

%% Collation and interpolation onto master clock
for i = 1:length(IMU)
    [IMU(i),GPS(i)] = collateIMUandGPS(IMU(i),GPS(i));
end

%% Wave processing
for i = 1:length(IMU)
    x  = GPS(i).x;
    y  = GPS(i).y;
    z  = IMU(i).pos(:,3);
    fs = round(IMU(i).samplingrate);

    SWIFT(i).time = median(IMU(i).time);
    SWIFT(i).lat  = median(GPS(i).lat);
    SWIFT(i).lon  = median(GPS(i).lon);

    [SWIFT(i).sigwaveheight, SWIFT(i).peakwaveperiod, SWIFT(i).peakwavedirT, ...
     SWIFT(i).wavespectra.energy, SWIFT(i).wavespectra.freq, ...
     SWIFT(i).wavespectra.a1, SWIFT(i).wavespectra.b1, ...
     SWIFT(i).wavespectra.a2, SWIFT(i).wavespectra.b2, ...
     SWIFT(i).wavespectra.check ] ...
        = XYZwaves(x,y,z,fs);
  
end

%% Save
save([outputDir,'microSWIFT002.mat'],'SWIFT')

%% scraps

% for gi = 1:length(GPSflist)
%         matfile = dir([GPSflist(gi).name(1:end-4) '.mat']);
%        if readraw | isempty(matfile)
%         disp(['GPS file ' num2str(gi) ' of ' num2str(length(GPSflist))])
%        GPS(gi) = processmicroSWIFT_GPS(fullfile(GPSflist(gi).folder,GPSflist(gi).name));
%             save([GPSflist(gi).name(1:end-4)],'GPS')
%             save([GPSflist(gi).name(1:end-4) '.mat'],'GPS','GPSsamplingrate');
%             save([GPSflist(end).name(1:13) '_' GPSflist(end).name(19:27) '_results'],'GPSresults');
%        else
%             load([GPSflist(gi).name(1:end-4) '.mat']);
%             disp('loading existing mat file')
%        end
% end

% for ii = 1:length(IMUflist)
%     disp(['IMU file ' num2str(ii) ' of ' num2str(length(IMUflist))])
%    matfile = dir([IMUflist(ii).name(1:end-4) '.mat']);
%    if readraw | isempty(matfile)
%       
%     IMU(ii) = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),'body');
%      
    %         save([IMUflist(ii).name(1:end-4) '.mat'],'IMU*')
%     else
%         load([IMUflist(ii).name(1:end-4) '.mat']) % load(matfile)?
%         disp('loading existing mat file')
%     end
% end