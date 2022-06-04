% reprocess SWIFT v3 wave results using a surface reconstruction
% and acounting for listing or capsizing during icing conditions
% loop thru raw data for a given SWIFT deployment, then
% replace values in the SWIFT data structure of results
% (assuming concatSWIFTv3_processed.m has already been run.
%
% J. Thomson, Oct 2015
%   edited version to add spectral check to reprocessed SWIFT data - Maddie 04/2016
%   cleaned and revised with IMU read function, Thomson, Jun 2016
%   subroutine to recalc all spectral moments based on displacements, Jun 2016
%   use RC filter in displacements, Oct 2016
%   revert to original directional moments, Oct 2017
%   optional filters and GPS reprocessing, Jun 2022

%% set up
clear all; close all
parentdir = pwd;  % change this to be the parent directory of all raw raw data (CF card offload from SWIFT)

%% choose a filter
%filter = str2cell('RC');
filter = str2cell('elliptic');
dB = 10;

%% load existing SWIFT structure created during concatSWIFTv3_processed, replace only the new wave results
cd(parentdir);
wd = pwd;
wdi = find(wd == '/',1,'last');
wd = wd((wdi+1):length(wd));

load([wd '.mat'])

prune = false(1,length(SWIFT)); % initialize logical array for later pruning of bad data

%cd('IMU/Raw/') % v3.2
cd('COM-6/Raw/') % v3.3


%% loop thru raw data

dirlist = dir('20*');

for di = 1:length(dirlist),
    
    cd([dirlist(di).name])
    filelist = dir('*.dat');
    
    for fi=1:length(filelist),
        
        % read or load raw data
        if isempty(dir([filelist(fi).name(1:end-4) '.mat'])),
            [ AHRS GPS ] = readSWIFTv3_IMU( filelist(fi).name );
        else
            load([filelist(fi).name(1:end-4) '.mat']),
        end
        
        % make sure there is data to work with
        if ~isempty(GPS) & ~isempty(AHRS) & length(AHRS.Accel) > 12000 & ~isempty(SWIFT),
            
            % find matching time index
            % use median to get burst time, because first entries are bad (no satellites acquired yet)
            time = nanmedian(datenum(GPS.UTC.Yr,GPS.UTC.Mo, GPS.UTC.Da, GPS.UTC.Hr, GPS.UTC.Mn, GPS.UTC.Sec));
            % match time to SWIFT structure and replace values
            [tdiff tindex] = min(abs([SWIFT.time]-time));
            if tdiff>1/48,
                disp('time gap too large at '),
                datestr(time)
                continue
            else
            end
            
            f = SWIFT(tindex).wavespectra.freq;  % original frequency bands
            
            % reconstruct sea surface by double integrating (and filtering) the accelerations
            [y,x,z, hs ] = rawdisplacements(AHRS, filter); % call is [y,x,z] to get output in east, north, up
            save([filelist(fi).name(1:end-4) '.mat'],'z','-APPEND')
            
            % make new scalar energy spectra from sea surface heights
            dt = median(diff(AHRS.Timestamp_sec));  % time step should be 0.04 s
            if isnan(dt),
                dt = 600 ./ length(AHRS.Accel);
            else
            end
            fs = 1/dt; % should be 25 Hz
            
            % make wave spectra based on non-nan displacements
            good = ~isnan( x + y + z);
            [ newHs, newTp, newDp, newE, newf, newa1, newb1, newa2, newb2, check ] = XYZwaves(x(good),y(good),z(good),fs);

            % also reprocess GPS raw data, noting that GPS waves already has an RC filter
            u = GPS.NED_Vel.Velocity_NED(end-2047:end,2); ... 
            v =  GPS.NED_Vel.Velocity_NED(end-2047:end,1);
            z_gps = GPS.Geodetic_Pos.H_above_MSL(end-2047:end); 
            fs_gps = 1000./median(diff(GPS.UTC.mSec));
            if strcmp( cellstr( filter ), 'elliptic')  % eliptic filter option
                [B,A] = ellip(3, .5, dB, 0.05/(fs_gps/2), 'high'); % original is ellip(3, .5, 20, 0.05/(fs/2), 'high');
                u = filtfilt(B, A, double(u));
                v = filtfilt(B, A, double(v));
                z_gps = filtfilt(B, A, double(z_gps));
            end
            [ Hs, Tp, Dp, E_gps, f_gps, a1, b1, a2, b2 ] = GPSwaves(u, v, z_gps, fs_gps ); 
            
            if newHs~=9999 & ~isnan(newHs), % only replace valid results
                
                % interp to the original freq bands
                E = interp1(newf,newE,f);
                Egps = interp1(f_gps,E_gps,f);
                a1 = interp1(newf,newa1,f);
                b1 = interp1(newf,newb1,f);
                a2 = interp1(newf,newa2,f);
                b2 = interp1(newf,newb2,f);
                
                % replace scalar values, but not directional moments
                SWIFT(tindex).sigwaveheight = newHs;
                SWIFT(tindex).peakwaveperiod = newTp;
                SWIFT(tindex).peakwaveperiod = newTp;
                %SWIFT(tindex).peakwavedirT = newDp;
                SWIFT(tindex).wavespectra.energy = E';
                SWIFT(tindex).wavespectra.energy_gps = E_gps';
                %SWIFT(tindex).wavespectra.a1 = a1;
                %SWIFT(tindex).wavespectra.b1 = b1;
                %SWIFT(tindex).wavespectra.a2 = a2;
                %SWIFT(tindex).wavespectra.b2 = b2;
                SWIFT(tindex).wavespectra.check = check';
                
                % include raw displacements (25 Hz)
                SWIFT(tindex).x = x;
                SWIFT(tindex).y = y;
                SWIFT(tindex).z = z;
                
                % include raw GPS velocities (4 Hz)
                gpslength = length(GPS.NED_Vel.Velocity_NED(:,2));
                last2048 = fliplr(gpslength - [0:2047]);
                if last2048 > 0 & isreal(last2048),
                    SWIFT(tindex).u = GPS.NED_Vel.Velocity_NED(last2048,2);
                    SWIFT(tindex).v = GPS.NED_Vel.Velocity_NED(last2048,1);
                    SWIFT(tindex).rawlat = GPS.Geodetic_Pos.Lat_Lon(last2048,1);
                    SWIFT(tindex).rawlon = GPS.Geodetic_Pos.Lat_Lon(last2048,2);
                else
                    SWIFT(tindex).u = NaN(2048,1);
                    SWIFT(tindex).v = NaN(2048,1);
                     SWIFT(tindex).rawlat = NaN(2048,1);
                    SWIFT(tindex).rawlon = NaN(2048,1);
                end
                
            else
                
                prune(tindex) = true; % set for pruning, b/c invalid wave result
                
            end
                        
        else
            
            % not enough raw data

        end
    end
    
    cd('../')
    
end

cd(parentdir)

% Quality control
SWIFT(prune) = [];

for si=1:length(SWIFT),
    if SWIFT(si).peakwavedirT > 9000 , 
        SWIFT(si).peakwavedirT = NaN;
    else
    end
end

%% save a big file with raw displacements, then a small file with stats only

save([ wd '_reprocessedIMU_' filter{1} 'filter_displacements.mat'],'SWIFT')

%% save a small file with stats only

SWIFT = rmfield(SWIFT,'x');
SWIFT = rmfield(SWIFT,'y');
SWIFT = rmfield(SWIFT,'z');
SWIFT = rmfield(SWIFT,'u');
SWIFT = rmfield(SWIFT,'v');
SWIFT = rmfield(SWIFT,'rawlat');
SWIFT = rmfield(SWIFT,'rawlon');

save([ wd '_reprocessedIMU_' filter{1} 'fitler.mat'],'SWIFT')

%% (re)plot

plotSWIFT(SWIFT)

%[Etheta theta f dir spread spread2 spread2alt ] = SWIFTdirectionalspectra(SWIFT, 1);

%% save the GPS results as their own structure
for si=1:length(SWIFT),
    SWIFT(si).wavespectra.energy = SWIFT(si).wavespectra.energy_gps;
end
save([ wd '_reprocessedGPS_' filter{1} 'fitler.mat'],'SWIFT')


