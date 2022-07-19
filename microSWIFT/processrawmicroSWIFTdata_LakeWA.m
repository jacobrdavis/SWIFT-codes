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
dataDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/microSWIFT/Deployment_1/';
% dataDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/microSWIFT/Deployment_2/';
outputDir = '/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/';

%% GPS processing

% Obtain GPS file list from specified data directory:
% GPSflist = dir([dataDir,microSwiftFolders{m},'/','*GPS*.dat']);
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
referenceFrame = 'earth'; %TRY EARTH!
filterType = 'RC'; %'butterworth_dynamic'; %'butterworth_highpass'; %'RC_dynamic'; %'butterworth_highpass'
detrendSignals = true; 

% Obtain IMU file list from specified data directory:
IMUflist = dir([dataDir,'*IMU*.dat']);

% Initialize IMU structure array:
clear IMU; IMU(length(IMUflist)) = processmicroSWIFT_IMU([],referenceFrame); 

% Process bursts:
for ii = 1:length(IMUflist)
    disp(['IMU file ' num2str(ii) ' of ' num2str(length(IMUflist))])

    if strcmp(referenceFrame,'body')
        [IMU(ii),filterFunc] = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),referenceFrame,filterType,detrendSignals);
    else
        [IMU(ii)] = processmicroSWIFT_IMU(fullfile(IMUflist(ii).folder,IMUflist(ii).name),referenceFrame);
    end
    
end  
clear ii
%%

% plot 1: vertical accelerations for foam
figure; hold on
plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).acc(:,3),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum'),IMU(5).acc(:,3),'Marker','*')
legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('acc z')
set(gca,'FontSize',14)
% print(gcf,[char(outputDir),'figures/','foam_vertical_accelerations','.png'],'-dpng');

% plot 2: vertical accelerations for foam, phase-shifted
figure; hold on
plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).acc(:,3),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum')+seconds(0.13),IMU(5).acc(:,3),'Marker','*')
legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('acc z'); title('phase shifted')
set(gca,'FontSize',14)

% plot 3: gyro about x for foam
figure; hold on
plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).gyro(:,1),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum'),IMU(5).gyro(:,1),'Marker','*')
legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('gyro x')
set(gca,'FontSize',14)

% plot 4: vertical accelerations for bottle-only, with foam for reference
figure; hold on
plot(datetime(IMU(1).time,'ConvertFrom','datenum'),IMU(1).acc(:,3),'Marker','*')
plot(datetime(IMU(4).time,'ConvertFrom','datenum'),IMU(4).acc(:,3),'Marker','*')
plot(datetime(IMU(3).time,'ConvertFrom','datenum'),IMU(3).acc(:,3),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum'),IMU(5).acc(:,3),'Marker','*')
legend({'microSWIFT040 (12Hz, SWIFT V4)','microSWIFT064 (12Hz, bottle-only)','microSWIFT043 (48Hz, bottle-only','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('acc z')
set(gca,'FontSize',14)

% plot 5: gyro about x for bottle-only, with foam for reference
figure; hold on
plot(datetime(IMU(4).time,'ConvertFrom','datenum'),IMU(4).gyro(:,1),'Marker','*')
plot(datetime(IMU(3).time,'ConvertFrom','datenum'),IMU(3).gyro(:,1),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum'),IMU(5).gyro(:,1),'Marker','*')
legend({'microSWIFT064 (12Hz, bottle-only)','microSWIFT043 (48Hz, bottle-only','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('gyro x')
set(gca,'FontSize',14)

% plot 6: z position foam
figure; hold on
plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).pos(:,3),'Marker','*')
plot(datetime(IMU(5).time,'ConvertFrom','datenum')+seconds(0.6),IMU(5).pos(:,3),'Marker','*')
legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
xlabel('time (s)'); ylabel('pos z'); title('phase shifted')
set(gca,'FontSize',14)

% plot 7 GPS time series
% GPS
i=3
figure
plot(datetime(GPS(i).time,'ConvertFrom','datenum'),GPS(i).v,'.')
GPSt = datetime(GPS(i).time,'ConvertFrom','datenum');
title('microSWIFT 043 (1L) in Lake WA')
ylabel('GPS u (m/s)')
ylim([-1 1])
xlim([GPSt(1) GPSt(end)])
xline(GPSt(1:100))
figure
hist(milliseconds(diff(GPSt)))
xlabel('dt (ms)')
ylabel('count')


%% trim IMU signal
for i = [1,2,4] %1:length(IMU) 
    signalLength = length(IMU(i).time);
    IMU(i).pos = IMU(i).pos(1:signalLength*0.7,:);
    IMU(i).time = IMU(i).time(1:signalLength*0.7);
    [IMU(i),GPS(i)] = collateIMUandGPS(IMU(i),GPS(i));
end

%% Collation and interpolation onto master clock
for i = 1:length(IMU) 
    [IMU(i),GPS(i)] = collateIMUandGPS(IMU(i),GPS(i));
end

%% Wave processing
for i = 1:length(IMU)
    id = extract(IMUflist(i).name,"microSWIFT"+digitsPattern(3));

% XYZwaves:
    x  = GPS(i).x;
    y  = GPS(i).y;
    z  = IMU(i).pos(:,3);
    fs = round(IMU(i).samplingrate);
    method = @() XYZwaves_microSWIFT(x,y,z,fs);

    mSWIFT.XYZwaves(i) = processWaves(method,initSWIFT(IMU(i),GPS(i),id));

% GPSwaves:
    u  = GPS(i).u;
    v  = GPS(i).v;
    z  = GPS(i).z;
    fs = round(IMU(i).samplingrate);
    method = @() GPSwaves_microswift(u,v,z,fs);

    mSWIFT.GPSwaves(i) = processWaves(method,initSWIFT(IMU(i),GPS(i),id));

% GPSandIMUwaves:
    u  = GPS(i).u;
    v  = GPS(i).v;
    z  = IMU(i).pos(:,3);
    fs = round(IMU(i).samplingrate);
    method = @() GPSandIMUwaves_microswift(u,v,z,fs);

    mSWIFT.GPSandIMUwaves(i) = processWaves(method,initSWIFT(IMU(i),GPS(i),id));


end

%% Evaluaion
% load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT22_15Mar2022.mat'); SWIFT22 = SWIFT;
% load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT23_15Mar2022.mat'); SWIFT23 = SWIFT;
% load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT24_15Mar2022.mat'); SWIFT24 = SWIFT;

load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT22_15Mar2022_reprocessedSBG.mat'); SWIFT22 = SWIFT;
load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT23_15Mar2022_reprocessedSBG.mat'); SWIFT23 = SWIFT;
load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT24_15Mar2022_reprocessedSBG.mat'); SWIFT24 = SWIFT;
load('/Users/jacob/Dropbox/Projects/microSWIFT/onboard_development/data/2022-03-15 Lake WA/SWIFT/SWIFT25_15Mar2022_reprocessedSBG.mat'); SWIFT25 = SWIFT; clear SWIFT


figure; hold on
si = 20; %20-1700, 25
datetime([SWIFT22(si).time],'ConvertFrom','datenum')
plot(SWIFT22(si).wavespectra.freq,SWIFT22(si).wavespectra.energy,'k','LineWidth',2,'DisplayName','SWIFT22')
plot(SWIFT23(si).wavespectra.freq,SWIFT23(si).wavespectra.energy,'color',[0.5 0.5 0.5],'LineWidth',2,'DisplayName','SWIFT23')
plot(SWIFT24(si).wavespectra.freq,SWIFT24(si).wavespectra.energy,'color',[0.8 0.8 0.8],'LineWidth',2,'DisplayName','SWIFT24')
set(gca,'YScale','log'); set(gca,'XScale','log')
xlabel('frequency (Hz)')
ylabel('Energy density (m^2/Hz)')

mSWIFTlabels = {'microSWIFT040 (12Hz, SWIFTv4)','microSWIFT042 (12Hz, foam)','microSWIFT043 (48Hz, bottle-only)','microSWIFT064 (12Hz, bottle-only)','microSWIFT066 (48Hz, foam)'};
colors = [0.1 0.7 0.85; 0.83 0.14 0.14; 1.00 0.54 0.00; 0.47 0.25 0.80; 0.25 0.80 0.54];
styles = {'-','--',':'};
methodNames = fields(mSWIFT);
for m = 1:length(methodNames)
    for i = 1:length(mSWIFT.(methodNames{m}))
        mSWIFT.(methodNames{m})(i).time
        plot(mSWIFT.(methodNames{m})(i).wavespectra.freq,mSWIFT.(methodNames{m})(i).wavespectra.energy,'LineStyle',styles{m},'LineWidth',1.5,'color',colors(i,:),'DisplayName',mSWIFTlabels{i})
        set(gca,'YScale','log'); set(gca,'XScale','log')
        xlabel('frequency (Hz)')
        ylabel('Energy density (m^2/Hz)')
    end
end

set(gca,'FontSize',14)

legend('Location','northeast')

% f = gcf; exportgraphics(f,[outputDir,datestr(mSWIFT(1).time,'ddmmmyyyy_HHMM'),'_LakeWA_scalar_spectra','.png'],'Resolution',600)

% legend({'SWIFT22','microSWIFT040 (12Hz, SWIFTv4','microSWIFT042 (12Hz, foam',...
%        'microSWIFT043 (48Hz, bottle-only)','microSWIFT064 (12Hz, bottle-only)','microSWIFT066 (48Hz, foam'},'Location','northeast')



%% Save
% figure; hold on
% for i = 1:length(SWIFT)
%     scatter(SWIFT(i).metrics.benchmarkTime,SWIFT(i).metrics.spectral.energy('RMSE_total'))
% end
% summary plot of RMSE for each method for every buoy (or combine buoys)!

% filenameStr = join([SWIFT(1).metrics.method,SWIFT(1).metrics.referenceFrame,SWIFT(1).metrics.filterType,SWIFT(1).id],'_'); % [char(method),'_',char(SWIFT.id),'_',char(CDIP.id),'_',filenamedate]
% save([outputDir,'matfiles/',char(filenameStr),'.mat'],'SWIFT')



%%
function [SWIFT] = initSWIFT(IMU,GPS,id)
   SWIFT.id = id;
   SWIFT.time = datetime(median(IMU.time),'ConvertFrom','datenum'); % datestr(median(IMU.time),'yyyy-mm-dd HH:MM:SS UTC');
   SWIFT.lat  = median(GPS.lat);
   SWIFT.lon  = median(GPS.lon);
   SWIFT.samplingrate = IMU.samplingrate;
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



% 
% figureVisibility = 'off';
% printFigures = true;
% 
% figure; hold on
% 
% plot(SWIFT22(24).wavespectra.freq,SWIFT22(24).wavespectra.energy,'k','LineWidth',2,'DisplayName','SWIFT22')
% 
% for i = 1:length(mSWIFT)
% 
%     plot(mSWIFT(i).wavespectra.freq,mSWIFT(i).wavespectra.energy,'LineWidth',1.5,'DisplayName',mSWIFT(i).id{1})
%     set(gca,'YScale','log','XScale','log')
%     xlabel('frequency (Hz)')
%     ylabel('Energy density (m^2/Hz)')
% end
% 
% legend('Location','northeast')


% %%
% fs = 200;
% 
% t0 = datetime(IMU(4).time,'ConvertFrom','datenum');
% tint = datenum(t0(1):seconds(1/fs):t0(end));
% t0 = datenum(t0);
% yspl = spline(IMU(4).time,IMU(4).acc(:,3),tint);
% 
% fc = 2;
% [b,a] = butter(3,fc/(fs/2),'low');
% ysplf = filtfilt(b,a,yspl);
% 
% dt = fs^(-1);
% 
% % integrate to vel
% vel = cumtrapz(ysplf)*dt ;
% vel_filt = filtfilt(b, a, vel);
% 
% % integrate to pos
% pos = cumtrapz(vel_filt)*dt ;
% pos_filt = filtfilt(b,a,pos);
% 
% [E1,f1] = pwelch(pos_filt,[],[],[],fs);
% 
% figure;
% loglog(f1, E1)
% hold on
% loglog(SWIFT22(si).wavespectra.freq,SWIFT22(24).wavespectra.energy,'k','LineWidth',2,'DisplayName','SWIFT22')
% 
% 
% figure; hold on
% plot(t0,IMU(4).acc(:,3))
% plot(tint,ysplf)
% 
% [E1,f1] = pwelch(ysplf,[],[],[],fs);
% [E2,f2] = pwelch(IMU(4).acc(:,3),[],[],[],IMU(4).samplingrate);
% 
% figure;
% loglog(f1, E1)
% hold on
% loglog(f2, E2)
% loglog(SWIFT22(si).wavespectra.freq,SWIFT22(24).wavespectra.energy,'k','LineWidth',2,'DisplayName','SWIFT22')
% 
% 
% figure; hold on
% plot(tint,yspl)
% plot(tint,ysplf)
% 
% var(yspl)
% var(ysplf)
% 
% 
% var(IMU(4).acc(:,3))
% var(IMU(3).acc(:,3))
% var(IMU(5).acc(:,3))
% 
% 
% 
% 
% figure; hold on
% plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).pos(:,3),'Marker','*')
% plot(datetime(IMU(5).time,'ConvertFrom','datenum')+seconds(0.6),IMU(5).pos(:,3),'Marker','*')
% legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
% xlabel('time (s)'); ylabel('pos z'); title('phase shifted')
% 
% figure; hold on
% plot(datetime(IMU(2).time,'ConvertFrom','datenum'),IMU(2).gyro(:,1),'Marker','*')
% plot(datetime(IMU(5).time,'ConvertFrom','datenum')+seconds(0.6),IMU(5).gyro(:,1),'Marker','*')
% legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
% xlabel('time (s)'); ylabel('gyro x'); title('phase shifted')
% 
% 
% [E1,f1] = pwelch(IMU(2).pos(:,3),[],[],[],IMU(2).samplingrate);
% [E2,f2] = pwelch(IMU(5).pos(:,3),[],[],[],IMU(5).samplingrate);
% figure; hold on
% plot(f1,E1,'Marker','*')
% plot(f2,E2,'Marker','*')
% legend({'microSWIFT042 (12Hz, foam','microSWIFT066 (48Hz, foam'})
% set(gca,'YScale','log','XScale','log')
% xlabel('frequency (Hz)')
% ylabel('Energy density (m^2/Hz)')
% 
% 
% figure; hold on
% for i = 1:length(IMU)
%     plot(IMU(i).time,IMU(i).pos(:,3))
% end
% xlabel('time (s)'); ylabel('z (m)')
% legend({'microSWIFT040 (12Hz, SWIFTv4','microSWIFT042 (12Hz, foam',...
%        'microSWIFT043 (48Hz, bottle-only)','microSWIFT064 (12Hz, bottle-only)','microSWIFT066 (48Hz, foam'})






