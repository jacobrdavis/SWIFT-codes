function [IMU] = processmicroSWIFT_IMU(IMUfile,referenceFrame)
%% processmicroSWIFT_IMU.m
% Function to process microSWIFT IMU data, in the form of a <.dat> file, 
% for a single recording burst.
%
% Usage:
%   [IMU] = processmicroSWIFT_IMU(IMUfile, mo, Wd);
%
% Inputs:
%   IMUfile - location of a raw IMU data file (in the form of a <.dat> file)
%   mo      - 1x3 matrix of magnetometer offsets [mxo, myo, mzo]
%   Wd      - weight coef for filtering the gyro
%
% Outputs: 
%   IMU - structure with fields:
%       time - UTC time (MATLAB datenum)
%       acc  - vector of accelerations [ax, ay, ax] (m/s^2)
%       mag  - vector of magnetometer readings [mx, my, mz] (uTesla)
%       gyro - vector of gyroscope readings [gx gy gz] (deg/s)
%       x    - buoy x-position 
%       y    - buoy y-position 
%       z    - buoy z-position, altitude
%       u    - buoy x-velocity (m/s)
%       v    - buoy y-velocity (m/s)
%       
%       (x,y, and z coordinates depend on the chosen reference frame, body or Earth) 
%
% Outputs will be '9999' for invalid results.
%
% Notes:
%   - The input weight coef Wd must be between 0 and 1, with 0 as default 
%     (this controls importantce dynamic angles in a complimentary filter)
%   - The default magnetomoter offsets are mxo = 60, myo = 60, mzo = 120
%   - The sampling rate is usually 4 Hz
%   - The body reference frame for the inputs is
%       x: along bottle (towards cap), roll around this axis
%       y: across bottle (right hand sys), pitch around this axis
%       z: up (skyward, same as GPS), yaw around this axis
%
% Dependencies: 
%   SWIFT-codes (readmicroSWIFT_IMU)
%
% J. Davis 2021-01-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020
%%

% checks
if exist('IMUfile','var') && isfile(IMUfile) % && IMUfile.bytes > 0
    
    % read raw IMU files:
    IMU = readmicroSWIFT_IMU(IMUfile, false);
    
    % Trim last entry in each field:
    IMU.acc(end,:)   = []; 
    IMU.mag(end,:)   = []; 
    IMU.gyro(end,:)  = [];
    IMU.clock(end)   = [];  
    IMU.time(end)    = []; 

    % Compute sample rate:
    IMU.samplingrate = length(IMU.acc)./((max(IMU.time)-min(IMU.time))*24*3600); % usually 12 Hz
    
    % Sort IMU onto master clock
    [IMU] = sortIMU(IMU);

    % Integrate to position
    switch referenceFrame
        case 'earth'
        % TODO: make these inputs:
        mxo = 60; myo = 60; mzo = 120; % magnetometer offsets
        Wd = .5;  % weighting in complimentary filter, 0 to 1

        % call IMUtoXYZ
        [IMU.x, IMU.y, IMU.z, IMU.roll, IMU.pitch, IMU.yaw, IMU.heading] = ...
            IMUtoXYZ(IMU.acc(:,1), IMU.acc(:,2), IMU.acc(:,3), IMU.gyro(:,1), IMU.gyro(:,2), ...
            IMU.gyro(:,3), IMU.mag(:,1), IMU.mag(:,2), IMU.mag(:,3), mxo, myo, mzo, Wd, IMUsamplingrate);

        case 'body'
        RC = 4;  % high pass RC filter constant, T > (2 * pi * RC)
        [IMU] = IMUtoXYZbody(IMU,RC);

    end
else % if conditions are not satisfied, report an empty structure (also useful for initialization)
    
    IMU = struct('clock',[],'acc',[],'mag',[],'gyro',[],'time',[],'samplingrate',[],'pos',[]); % [],'z',[],'u',[],'v',[],'x',[],'y',[]

end

end

function [IMU] = IMUtoXYZbody(IMU,RC)
% filter and integrate linear accelerations to get linear velocities
ax = IMU.acc(:,1);
ay = IMU.acc(:,2);
az = IMU.acc(:,3);
t  = IMU.time;
fs = IMU.samplingrate;
dt = fs^(-1);

ax = detrend(ax);
ay = detrend(ay);
az = detrend(az);

axf = RCfilter(ax, RC, fs);
ayf = RCfilter(ay, RC, fs);
azf = RCfilter(az, RC, fs);

% figure(); 
% subplot(3,1,1)
%     plot(t,ax,'color','b','DisplayName','unfiltered'); hold on
%     plot(t,axf,'color','r','DisplayName','filtered')
%     ylabel('ax (m/s)')
%     legend()
% subplot(3,1,2)
%     plot(t,ay,'color','b') ; hold on
%     plot(t,ayf,'color','r')
%     ylabel('ay (m/s)')
% subplot(3,1,3)
%     plot(t,az,'color','b') ; hold on
%     plot(t,azf,'color','r')
%     ylabel('az (m/s)')
% 
% figure;
% [f,Ma,~] = pltFFT(az,fs);
% semilogx(f,Ma,'color','b'); hold on 
% [f,Ma,~] = pltFFT(azf,fs);
% semilogx(f,Ma,'color','r'); hold on 
% xline(0.05,'LineWidth',1.5); xline(0.5,'LineWidth',1.5); 
% xline(1/(2*pi*RC),'LineWidth',1.5,'color','r')
% xlabel('f (Hz)')
% ylabel('FFT(az) (m/s^2)')

vx = cumtrapz(axf)*dt; % m/s
vy = cumtrapz(ayf)*dt; % m/s
vz = cumtrapz(azf)*dt; % m/s

vx = detrend(vx,'omitnan');
vy = detrend(vy,'omitnan');
vz = detrend(vz,'omitnan');

vxf = RCfilter(vx, RC, fs);
vyf = RCfilter(vy, RC, fs);
vzf = RCfilter(vz, RC, fs);

% figure(); 
% subplot(3,1,1)
%     plot(t,vx,'color','b','DisplayName','unfiltered'); hold on
%     plot(t,vxf,'color','r','DisplayName','filtered')
%     ylabel('vx (m/s)')
%     legend()
% subplot(3,1,2)
%     plot(t,vy,'color','b') ; hold on
%     plot(t,vyf,'color','r')
%     ylabel('vy (m/s)')
% subplot(3,1,3)
%     plot(t,vz,'color','b') ; hold on
%     plot(t,vzf,'color','r')
%     ylabel('vz (m/s)')
% 
% figure;
% [f,Ma,~] = pltFFT(vz,fs);
% semilogx(f,Ma,'color','b'); hold on 
% [f,Ma,~] = pltFFT(vzf,fs);
% semilogx(f,Ma,'color','r'); hold on 
% xline(0.05); xline(0.5); 
% xline(1/(2*pi*RC),'color','r')
% xlabel('f (Hz)')
% ylabel('FFT(vz) (m/s)')

x = cumtrapz(vxf)*dt;
y = cumtrapz(vyf)*dt;
z = cumtrapz(vzf)*dt;

x = detrend(x,'omitnan');
y = detrend(y,'omitnan');
z = detrend(z,'omitnan');

xf = RCfilter(x, RC, fs);
yf = RCfilter(y, RC, fs);
zf = RCfilter(z, RC, fs);
% 
% figure(); 
% subplot(3,1,1)
%     plot(t,x,'color','b','DisplayName','unfiltered'); hold on
%     plot(t,xf,'color','r','DisplayName','filtered')
%     ylabel('x (m)')
%     legend()
% subplot(3,1,2)
%     plot(t,y,'color','b') ; hold on
%     plot(t,yf,'color','r')
%     ylabel('y (m)')
% subplot(3,1,3)
%     plot(t,z,'color','b') ; hold on
%     plot(t,zf,'color','r')
%     ylabel('z (m)')
% 
% figure;
% [f,Ma,~] = pltFFT(z,fs);
% semilogx(f,Ma,'color','b'); hold on 
% [f,Ma,~] = pltFFT(zf,fs);
% semilogx(f,Ma,'color','r'); hold on 
% xline(0.05); xline(0.5); 
% xline(1/(2*pi*RC),'color','r')
% xlabel('f (Hz)')
% ylabel('FFT(z) (m)')

% remove first portion, which can has initial oscillations from filtering
xf(1:round(RC./dt*10),:) = 0;
yf(1:round(RC./dt*10),:) = 0;
zf(1:round(RC./dt*10),:) = 0;



% figure()
% subplot(3,1,1)
%     plot(t,xf)
%     ylabel('x (m)')
% subplot(3,1,2)
%     plot(t,yf)
%     ylabel('y (m)')
% subplot(3,1,3)
%     plot(t,zf)
%     ylabel('z (m)')
% drawnow

IMU.pos = [xf yf zf];

end

%% EMBEDDED RC FILTER function (high pass filter) %%
function a = RCfilter(b, RC, fs)

alpha = RC / (RC + 1./fs);
a = b;

for ui = 2:length(b)
    a(ui) = alpha * a(ui-1) + alpha * ( b(ui) - b(ui-1) );
end

end


function [f,Ma,Ph] = pltFFT(y,fs)
Y = fft(y);
L = length(y);

% magnitude
Ma2 = abs(Y/L);
Ma = Ma2(1:round(L/2+1));
Ma(2:end-1) = 2*Ma(2:end-1);

% phase
Ph2 = angle(Y);
Ph = Ph2(1:round(L/2+1));

% frequency
f = fs*(0:round(L/2))/L;
end