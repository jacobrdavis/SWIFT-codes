%% processrawmicroSWIFTdata.m
% Script to process microSWIFT raw data from GPS and IMU; 
% 
% J. Davis, 2022-01-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020
%
% Dependencies: 
%   SWIFT-codes (processmicroSWIFT_GPS,processmicroSWIFT_IMU,collateIMUandGPS)
%   Mapping Toolbox (deg2km)

%% Settings and inputs:
% runs to do: body frame XYZwaves, GPSwaves, GPSandIMUwaves
                    % RC@0.4 and Dynamic cutoff? ...maybe not
% XYZwaves - RC done!
% XYZwaves - bwdyn done!
% GPSwaves - RC done!
% GPSwaves - bwdyn done!
% GPSandIMUWaves - RC
% GPSandIMUWaves - bwdyn

% dataDir   = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/StaticMobileBay/';
% dataDir  = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/DuckFRF/microSWIFT002/';
dataDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/Benchtop_2021-03-14/';
outputDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/';

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
clear gi;

%% IMU processing
% Methods
referenceFrame = 'body';
filterType = 'RC';%'butterworth_dynamic'; %'butterworth_highpass'; %'RC_dynamic'; %'butterworth_highpass'
detrendSignals = true; 

% Obtain IMU file list from specified data directory:
IMUflist = dir([dataDir,'*IMU*.dat']);

% Initialize IMU structure array:
clear IMU; IMU(length(IMUflist)) = processmicroSWIFT_IMU([],referenceFrame); 

% Process bursts:
for ii = 1:length(IMUflist)
    disp(['IMU file ' num2str(ii) ' of ' num2str(length(IMUflist))])
    if strcmp(filterType,'butterworth_dynamic')
        fc = dynamicFilterCutoff(sqrt(GPS(ii).u.^2+GPS(ii).v.^2),GPS(ii).samplingrate,0.70,true);
        [IMU(ii),filterFunc] = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),referenceFrame,filterType,detrendSignals,fc);
        
        GPS(ii).x = filterFunc(GPS(ii).x);
        GPS(ii).y = filterFunc(GPS(ii).y);
        GPS(ii).u = filterFunc(GPS(ii).u);
        GPS(ii).v = filterFunc(GPS(ii).v);
        disp('dyn')
    else
        [IMU(ii),filterFunc] = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),referenceFrame,filterType,detrendSignals);
        disp('not dyn')

        GPS(ii).x = filterFunc(GPS(ii).x);
        GPS(ii).y = filterFunc(GPS(ii).y);
        GPS(ii).u = filterFunc(GPS(ii).u);
        GPS(ii).v = filterFunc(GPS(ii).v);

    end
end 
clear ii
%%
figure
plot(IMU.time,IMU.acc(:,2:3))

figure
plot(IMU.time,IMU.pos(:,2:3))

%%
figure
plot(IMU.time,IMU.gyro)

%% Collation and interpolation onto master clock
for i = 1:length(IMU) 
    [IMU(i),GPS(i)] = collateIMUandGPS(IMU(i),GPS(i));
end

%% Wave processing
for i = 1:length(IMU)
    id = extract(IMUflist(i).name,"microSWIFT"+digitsPattern(3));
% XYZwaves:
%     x  = GPS(i).x;
%     y  = GPS(i).y;
%     z  = IMU(i).pos(:,3);
%     fs = round(IMU(i).samplingrate);
%     method = @() XYZwaves_microSWIFT(x,y,z,fs);

% GPSwaves:
%     u  = GPS(i).u;
%     v  = GPS(i).v;
%     z  = GPS(i).z;
%     fs = round(IMU(i).samplingrate);
%     method = @() GPSwaves_microswift(u,v,z,fs);

% GPSandIMUwaves:
    u  = GPS(i).u;
    v  = GPS(i).v;
    z  = IMU(i).pos(:,3);
    fs = round(IMU(i).samplingrate);
    method = @() GPSandIMUwaves_microswift(u,v,z,fs);
    SWIFT(i) = processWaves(method,initSWIFT(IMU(i),GPS(i),id));
end
%%
% figure; hold on
% for i = 1:length(SWIFT)
% plot(SWIFT(i).wavespectra.freq,SWIFT(i).wavespectra.energy,'color',0.6*[1 1 1])
% set(gca,'YScale','log','XScale','log')
% %     legend('Location','northeast')
%     xlabel('frequency (Hz)')
%     ylabel('Energy density (m^2/Hz)')
% end
% plot([0.01 0.5],10^(-8).*[0.01 0.5].^(-4),'k','LineWidth',2)
% set(gca,'FontSize',14)
% 
% figure
% [Ez,fz] = pwelch(IMU(3).pos(:,3),[],[],[],fs);
% [Ex,fx] = pwelch(GPS(3).x,[],[],[],fs);
% [Ey,fy] = pwelch(GPS(3).y,[],[],[],fs);
% plot(f,E,'color',0.6*[1 1 1]); hold on
% plot(f,Ex+Ey,'color',0*[1 1 1])
% set(gca,'YScale','log','XScale','log')
%     xlabel('frequency (Hz)')
%     ylabel('Energy density (m^2/Hz)')

%% Evaluation
load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/DuckFRF/CDIP/CDIP192_Oct2021.mat')
load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/DuckFRF/CDIP/CDIP433_Oct2021.mat')

figureVisibility = 'off';
printFigures = true;

for i = 1:length(SWIFT)
    if day(SWIFT(i).time) == 26
        CDIP = CDIP433;
    elseif day(SWIFT(i).time) == 27
        CDIP = CDIP192;
    end

    [~,j]=min(abs([CDIP.time] - SWIFT(i).time));
    SWIFT(i).metrics= compute_evaluation_metrics(SWIFT(i),CDIP(j),initMetrics(func2str(method),referenceFrame,filterType,IMU(i).cutoff),outputDir,figureVisibility,printFigures);
    
%     in another file, plot RMSE as a function of dataset for each method (load in files)
end


%% Save
% figure; hold on
% for i = 1:length(SWIFT)
%     scatter(SWIFT(i).metrics.benchmarkTime,SWIFT(i).metrics.spectral.energy('RMSE_total'))
% end
% summary plot of RMSE for each method for every buoy (or combine buoys)!

filenameStr = join([SWIFT(1).metrics.method,SWIFT(1).metrics.referenceFrame,SWIFT(1).metrics.filterType,SWIFT(1).id],'_'); % [char(method),'_',char(SWIFT.id),'_',char(CDIP.id),'_',filenamedate]
save([outputDir,'matfiles/',char(filenameStr),'.mat'],'SWIFT')
;
%%
function [SWIFT] = initSWIFT(IMU,GPS,id)
   SWIFT.id = id;
   SWIFT.time = datetime(median(IMU.time),'ConvertFrom','datenum'); % datestr(median(IMU.time),'yyyy-mm-dd HH:MM:SS UTC');
   SWIFT.lat  = median(GPS.lat);
   SWIFT.lon  = median(GPS.lon);
end

function metrics = initMetrics(methodStr,referenceFrame,filterType,cutOffFrequency)
    if contains(methodStr,'@') % remove function handle characters
        methodStr = extractBetween(methodStr,')','(');
        if contains(methodStr,'_')
            methodStrTemp = split(methodStr,'_');
            methodStr = methodStrTemp{1};
        end
    end
    metrics.method = methodStr;
    metrics.referenceFrame = referenceFrame;
    metrics.filterType = filterType;
    metrics.cutOffFrequency = cutOffFrequency;
end

function [SWIFT] = processWaves(method,SWIFT)

    [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check ]= method();

    SWIFT.sigwaveheight = Hs;
    SWIFT.peakwaveperiod = Tp;
    SWIFT.peakwavedirT = Dp;
    SWIFT.wavespectra.energy = E;
    SWIFT.wavespectra.freq = f;
    SWIFT.wavespectra.a1 = a1;
    SWIFT.wavespectra.b1 = b1;
    SWIFT.wavespectra.a2 = a2;
    SWIFT.wavespectra.b2 = b2;
    SWIFT.wavespectra.check = check;

end


%% Sensor fusion
% burst = 2;
% filt = insfilterNonholonomic('IMUSampleRate',12);
% filt.DecimationFactor = 2;
% filt.ZeroVelocityConstraintNoise =1;
% imuSamplesPerGPS = round(IMU(burst).samplingrate)/round(GPS(burst).samplingrate);
% 
% time = datetime([IMU(burst).time],'ConvertFrom','datenum');
% 
% gpsdt = 1/12; 
% gpsw = [diff(GPS(burst).z)/gpsdt; 0];
% 
% accelData = IMU(burst).acc;
% gyroData = IMU(burst).gyro;
% gpsLLA = [GPS(burst).lon GPS(burst).lat GPS(burst).z]; % geodetic latitude, longitude, and altitude (LLA); try x and y? 
% filt.ReferenceLocation = gpsLLA(1,:);
% gpsVel = [GPS(burst).u GPS(burst).v gpsw]; % gradient(GPS(burst).z,gpsdt)
% positionCov = 0;
% velocityCov = 0;
% numIMUSamples = size(accelData,1);
% estOrient = quaternion.ones(numIMUSamples,1);
% estPos = zeros(numIMUSamples,3);
%     
% gpsIdx = 1;
% for idx = 1:numIMUSamples
%     predict(filt,accelData(idx,:),gyroData(idx,:));       %Predict filter state
%     if (mod(idx,imuSamplesPerGPS) == 0)                   %Correct filter state
%       % fusegps(filt,gpsLLA(gpsIdx,:),Rpos,gpsVel(gpsIdx,:),Rvel);
%         fusegps(filt,gpsLLA(gpsIdx,:),positionCov,gpsVel(gpsIdx,:),velocityCov);
%         gpsIdx = gpsIdx + 1;
%     end
%     
%     [estPos(idx,:),estOrient(idx,:)] = pose(filt);        %Log estimated pose
% end
% 
% figure
% subplot(3,1,1)
% plot(time,estPos(:,1),'DisplayName','x'); hold on
% plot(time,GPS(burst).x)
% subplot(3,1,2)
% plot(time,estPos(:,2),'DisplayName','y')
% subplot(3,1,3)
% plot(time,estPos(:,3),'DisplayName','z')

% accelData   = [0 0 9.8];
% gyroData    = [0 0 0];
% position    = [0 0 0];
% positionCov = 0;
% velocity    = rand(1,3);
% velocityCov = 0;
% predictDataSampleRate = 100;
% fuseDataSampleRate = 2;
% predictSamplesPerFuse = predictDataSampleRate/fuseDataSampleRate;
% duration = 5;
% for i = 1:duration*fuseDataSampleRate
%     
%     for j = 1:predictSamplesPerFuse
%         
%         predict(FUSE,accelData,gyroData);
%         
%     end
%     
%     fusegps(FUSE,position,positionCov,velocity,velocityCov);
%     
% end

% [Ek,fk] = pwelch(estPos(:,3),[],[],[],fs);







