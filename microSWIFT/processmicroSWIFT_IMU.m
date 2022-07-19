function [IMU,filterFunc] = processmicroSWIFT_IMU(IMUfile,referenceFrame,filterType,detrendSignals,fc)
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
% J. Davis 2022-01-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020
%%

% checks
if exist('IMUfile','var') && ~isempty(IMUfile) && isfile(IMUfile) % && IMUfile.bytes > 0
    
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
    
%     % Record mean time of the current burst
%     IMU.meanBurstTime = datetime(mean(IMU.time),'ConvertFrom','datenum');

    % Integrate to position
    switch referenceFrame
        case 'earth'
        % TODO: make these inputs:
        mxo = 60; myo = 60; mzo = 120; % magnetometer offsets
        Wd = .5;  % weighting in complimentary filter, 0 to 1

        % call IMUtoXYZ
        [x,y,z,vx,vy,vz,roll,pitch,yaw,heading] = ...
        IMUtoXYZ(-IMU.acc(:,3), IMU.acc(:,2), IMU.acc(:,1), -IMU.gyro(:,3), IMU.gyro(:,2), ...
            IMU.gyro(:,1), -IMU.mag(:,3), IMU.mag(:,2), IMU.mag(:,1), mxo, myo, mzo, Wd, IMU.samplingrate);
        
        IMU.pos = [x y z];
        IMU.vel = [vx vy vz];
        IMU.angles = [roll pitch yaw];
        IMU.heading = heading;

        case 'body'
        
        if ~exist('fc',"var")
            fc = [];
        end

        [IMU,filterFunc] = IMUtoXYZbody(IMU,filterType,detrendSignals,fc);

    end
else % if conditions are not satisfied, report an empty structure (also useful for initialization)
    switch referenceFrame
        case 'earth'
            IMU = struct('clock',[],'acc',[],'mag',[],'gyro',[],'time',[],'samplingrate',[],'pos',[],'vel',[],'angles',[],'heading',[]); % [],'z',[],'u',[],'v',[],'x',[],'y',[]
        case 'body'
            IMU = struct('clock',[],'acc',[],'mag',[],'gyro',[],'time',[],'samplingrate',[],'cutoff',[],'pos',[],'vel',[]); % [],'z',[],'u',[],'v',[],'x',[],'y',[]
    end
end

end
