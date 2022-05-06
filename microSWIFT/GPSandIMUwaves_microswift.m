function [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = GPSandIMUwaves(u,v,z,fs) 

% matlab function to read and process GPS and IMU data
%   to estimate wave height, period, direction, and spectral moments
%   assuming deep-water limit of surface gravity wave dispersion relation
%
% Inputs are east velocity [m/s], north velocity [m/s], vertical 
% accelerations [g], pitch [rad], roll [rad], sampling rate [Hz]
%
% Some inputs can be empty variables, in which case the algorithm will use
% whatever non-empty inputs are available, with preference for GPS
% velocities
%
% Required input is sampling rate, which must be at least 1 Hz and the same
% for all variables.  Additionaly, non-empty input time series data must 
% have at least 512 points and all be the same size.
%
% Outputs are significat wave height [m], dominant period [s], dominant direction 
% [deg T, using meteorological from which waves are propagating], spectral 
% energy density [m^2/Hz], frequency [Hz], and 
% the normalized spectral moments a1, b1, a2, b2, 
%
% Outputs will be '9999' for invalid results.
%
% Outputs can be supressed, in order, thus full usage is as follows:
%
%   [ Hs, Tp, Dp, E, f, a1, b1, a2, b2 ] = GPSandIMUwaves(u,v,az,pitch,roll,fs); 
%
% and minimal usage is:
%
%   [ Hs, Tp ] = GPSandIMUwaves(u,v,[],[],[],fs); 
%
%
% J. Thomson,  12/2013, v1, modified from GPSwaves.m
%               1/2014, v2, include pitch and roll to correct vertical accelerations
%               4/2014, v3, N/A 
%               6/2014, v4, high-pass filtering with RC time constant = 20s
%                           correct IMU orientation (+ z acc has proper phase)
%                           use velocity spectra to select peak period
%                           correct direction estimate (rotation error)
%                           use GPS for all scalar spectra (less noise)
%                           changes spectra parameters to report higher f
%                               while keeping same # (=42) of frequencies
%               12/2014, v5, fix directions by using cospectra of az and u v
%                           (rather than quadspectrum, which is only appropriate for displacements                           directional estimate when not upright
%               12/2015, v6, remove the low-freq cutoff for spectra, but keep it for bulk stats
%
%               10/2017, v7, change the RC filter parameter to 3.5, after
%                           realizing that the cuttoff period is 2 pi * RC, not RC
%             
%              
%#codegen

%% tunable parameters

% low frequency noise ratio tolerance
% NOT APPLIED IN v6
    %LFNR = 4 ; 
    
% standard deviations for despiking        
Nstd = 10; 

% time constant [s] for high-pass filter, 
%     RC = 3.5; % cutoff period is 2*pi*RC

% energy ratios (unused as of version 3)
%maxEratio = 5; % max allowed ratio of Ezz to Exx + Eyy, default is 5
%minEratio = .1; % min allowed ratio of Ezz to Exx + Eyy, default is 0.1
    

%% fixed parameters
wsecs = 256;   % window length in seconds, should make 2^N samples
merge = 3;      % freq bands to merge, must be odd?
maxf = 1;       % frequency cutoff for telemetry Hz



    

%% rotate accelerations to realworld coordinates using Euler angles of pitch, roll, yaw
% where IMU pitch is centered at zero and roll is centered at +/- pi
% and note that IMU ouput on SWIFT is in radians 
% IMU is mounted upside down, so gravity is +1 g
% yaw not used (irrelevant for vertical)
% flip = -1;  % SWIFT specific, this accounts for phase relation of heave to vert acceleration

% if angles available, and clearly in radians (by variance), correct acceleration
% if ~isempty(pitch) & ~isempty(roll) & var(abs(pitch)) < 1 & var(abs(roll)) < 1,
%     z =  flip .* z ./ ( cos(pitch).*cos(roll) ) ; % approximation for small tilts (sin is small)
% else
%     z = flip .* z;
% end
    
%% Quality control inputs (despike)

badu = abs(detrend(u)) >= Nstd * std(u); % logical array of indices for bad points
badv = abs(detrend(v)) >= Nstd * std(v); % logical array of indices for bad points
u(badu) = mean( u(~badu) ); %sum(badu)
v(badv) = mean( v(~badv) ); %sum(badv)

%% begin processing, if data sufficient
pts = length(u);       % record length in data points

if pts >= 2*wsecs & fs>=1 & sum(badu)<100 & sum(badv)<100,  % minimum length and quality for processing

    
%% high-pass RC filter, detrend first

% u = detrend(u);
% v = detrend(v);
% z = detrend(z);
% 
% % initialize
% ufiltered = u;
% vfiltered = v;
% zfiltered = z;
% 
% alpha = RC / (RC + 1./fs); 
% 
% for ui = 2:length(u),
%     ufiltered(ui) = alpha * ufiltered(ui-1) + alpha * ( u(ui) - u(ui-1) );
%     vfiltered(ui) = alpha * vfiltered(ui-1) + alpha * ( v(ui) - v(ui-1) );
%     zfiltered(ui) = alpha * zfiltered(ui-1) + alpha * ( z(ui) - z(ui-1) );
% end
% 
% u = ufiltered;
% v = vfiltered;
% z = zfiltered;

%% break into windows (use 75 percent overlap)
w = round(fs * wsecs);  % window length in data points
if rem(w,2)~=0, w = w-1; else end  % make w an even number
windows = floor( 4*(pts/w - 1)+1 );   % number of windows, the 4 comes from a 75% overlap
dof = 2*windows*merge; % degrees of freedom
% loop to create a matrix of time series, where COLUMN = WINDOW 
uwindow = zeros(w,windows);
vwindow = zeros(w,windows);
zwindow = zeros(w,windows);
for q=1:windows 
	uwindow(:,q) = u(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
	vwindow(:,q) = v(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
  	zwindow(:,q) = z(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
end

%% detrend individual windows (full series already detrended)
for q=1:windows
uwindow(:,q) = detrend(uwindow(:,q));
vwindow(:,q) = detrend(vwindow(:,q));
zwindow(:,q) = detrend(zwindow(:,q));
end

%% taper and rescale (to preserve variance)
% form taper matrix (columns of taper coef)
taper = sin ( (1:w) * pi/w )' * ones(1,windows); 
% taper each window
uwindowtaper = uwindow .* taper;
vwindowtaper = vwindow .* taper;
zwindowtaper = zwindow .* taper;
% now find the correction factor (comparing old/new variance)
factu = sqrt( var(uwindow) ./ var(uwindowtaper) );
factv = sqrt( var(vwindow) ./ var(vwindowtaper) );
factz = sqrt( var(zwindow) ./ var(zwindowtaper) );
% and correct for the change in variance
% (mult each window by it's variance ratio factor)
uwindowready = (ones(w,1)*factu).* uwindowtaper;
vwindowready = (ones(w,1)*factv).* vwindowtaper;
zwindowready = (ones(w,1)*factz).* zwindowtaper;


%% FFT
% calculate Fourier coefs
Uwindow = fft(uwindowready);
Vwindow = fft(vwindowready);
Zwindow = fft(zwindowready);
% second half of fft is redundant, so throw it out
Uwindow( (w/2+1):w, : ) = [];
Vwindow( (w/2+1):w, : ) = [];
Zwindow( (w/2+1):w, : ) = [];
% throw out the mean (first coef) and add a zero (to make it the right length)  
Uwindow(1,:)=[]; Vwindow(1,:)=[]; Zwindow(1,:)=[]; 
Uwindow(w/2,:)=0; Vwindow(w/2,:)=0; Zwindow(w/2,:)=0; 
% POWER SPECTRA (auto-spectra)
UUwindow = real ( Uwindow .* conj(Uwindow) );
VVwindow = real ( Vwindow .* conj(Vwindow) );
ZZwindow = real ( Zwindow .* conj(Zwindow) );
% CROSS-SPECTRA 
UVwindow = ( Uwindow .* conj(Vwindow) );
UZwindow = ( Uwindow .* conj(Zwindow) );
VZwindow = ( Vwindow .* conj(Zwindow) );
disp('MAKE SURE CORRECT PHASING')

%% merge neighboring freq bands (number of bands to merge is a fixed parameter)
% initialize
UUwindowmerged = zeros(floor(w/(2*merge)),windows);
VVwindowmerged = zeros(floor(w/(2*merge)),windows);
ZZwindowmerged = zeros(floor(w/(2*merge)),windows);
UVwindowmerged = 1i*ones(floor(w/(2*merge)),windows);
UZwindowmerged = 1i*ones(floor(w/(2*merge)),windows);
VZwindowmerged = 1i*ones(floor(w/(2*merge)),windows);

for mi = merge:merge:(w/2) 
	UUwindowmerged(mi/merge,:) = mean( UUwindow((mi-merge+1):mi , : ) );
	VVwindowmerged(mi/merge,:) = mean( VVwindow((mi-merge+1):mi , : ) );
   	ZZwindowmerged(mi/merge,:) = mean( ZZwindow((mi-merge+1):mi , : ) );
	UVwindowmerged(mi/merge,:) = mean( UVwindow((mi-merge+1):mi , : ) );
  	UZwindowmerged(mi/merge,:) = mean( UZwindow((mi-merge+1):mi , : ) );
	VZwindowmerged(mi/merge,:) = mean( VZwindow((mi-merge+1):mi , : ) );
end
% freq range and bandwidth
n = (w/2) / merge;                         % number of f bands
Nyquist = .5 * fs;                % highest spectral frequency 
bandwidth = Nyquist/n ;                    % freq (Hz) bandwitdh
% find middle of each freq band, ONLY WORKS WHEN MERGING ODD NUMBER OF BANDS!
f = 1/(wsecs) + bandwidth/2 + bandwidth.*(0:(n-1)) ; 


%% ensemble average windows together
% take the average of all windows at each freq-band
% and divide by N*samplerate to get power spectral density
% the two is b/c Matlab's fft output is the symmetric FFT, and we did not use the redundant half (so need to multiply the psd by 2)
UU = mean( UUwindowmerged.' ) / (w/2 * fs  );
VV = mean( VVwindowmerged.' ) / (w/2 * fs  );
ZZ = mean( ZZwindowmerged.' ) / (w/2 * fs  );
UV = mean( UVwindowmerged.' ) / (w/2 * fs  ); 
UZ = mean( UZwindowmerged.' ) / (w/2 * fs  ); 
VZ = mean( VZwindowmerged.' ) / (w/2 * fs  ); 


%% convert to displacement spectra (from velocity and acceleration)
% assumes perfectly circular deepwater orbits
% could be extended to finite depth by calling wavenumber.m 
Exx = ( UU )  ./ ( (2*pi*f).^2 );  %[m^2/Hz]
Eyy = ( VV )  ./ ( (2*pi*f).^2 );  %[m^2/Hz]
Ezz = ( ZZ ); %./ ( (2*pi*f).^4 ) .* (9.8^2);  %[m^2/Hz]

% Qxz = imag(UZ) ./ ( (2*pi*f).^3 ) .* (9.8); %[m^2/Hz], quadspectrum of vertical acc and horizontal velocities
% Cxz = real(UZ) ./ ( (2*pi*f).^3 ) .* (9.8); %[m^2/Hz], cospectrum of vertical acc and horizontal velocities
% 
% Qyz = imag(VZ) ./ ( (2*pi*f).^3 ) .* (9.8); %[m^2/Hz], quadspectrum of vertical acc and horizontal velocities
% Cyz = real(VZ) ./ ( (2*pi*f).^3 ) .* (9.8); %[m^2/Hz], cospectrum of vertical acc and horizontal velocities
% 
% Cxy = real(UV) ./ ( (2*pi*f).^2 );  %[m^2/Hz]

Qxz = imag(UZ) ./ ( (2*pi*f).^1 ); %[m^2/Hz], quadspectrum of vertical displacement and horizontal velocities
Cxz = real(UZ) ./ ( (2*pi*f).^1 ); %[m^2/Hz], cospectrum of vertical displacement and horizontal velocities

Qyz = imag(VZ) ./ ( (2*pi*f).^1 ); %[m^2/Hz], quadspectrum of vertical displacement and horizontal velocities
Cyz = real(VZ) ./ ( (2*pi*f).^1 ); %[m^2/Hz], cospectrum of vertical displacement and horizontal velocities

Cxy = real(UV) ./ ( (2*pi*f).^2 );  %[m^2/Hz]


%% wave spectral moments 
% wave directions from Kuik et al, JPO, 1988 and Herbers et al, JTech, 2012

a1 = Cxz ./ sqrt( (Exx+Eyy).* Ezz );  %[], would use Qxz for actual displacements
b1 = Cyz ./ sqrt( (Exx+Eyy).* Ezz );  %[], would use Qyz for actual displacements
a2 = (Exx - Eyy) ./ (Exx + Eyy);
b2 = 2 .* Cxy ./ ( Exx + Eyy );

%% wave directions
% note that 0 deg is for waves headed towards positive x (EAST, right hand system)
dir1 = atan2(b1,a1) ;  % [rad], 4 quadrant
dir2 = atan2(b2,a2)/2 ; % [rad], only 2 quadrant
spread1 = sqrt( 2 * ( 1 - sqrt(a1.^2 + b2.^2) ) );
spread2 = sqrt( abs( 0.5 - 0.5 .* ( a2.*cos(2.*dir2) + b2.*cos(2.*dir2) )  ));



%% use orbit shape as check on quality
check = Ezz ./ (Eyy + Exx);

%% apply LFNR tolerance 
%Exx(LFNR*(UU) < Exx ) = 0;  % quality control based on LFNR of swell
%Eyy(LFNR*(VV) < Eyy ) = 0;  % quality control based on LFNR of swell
%Ezz(LFNR*(AZAZ.* (9.8^2)) < Ezz ) = 0;  % quality control based on LFNR of swell


%% Scalar energy spectra (a0)

% E = Exx + Eyy;
E = Ezz;

%E = zeros(1,length(f));
%if azdummy ==1,
%    E = Exx + Eyy;
%else
%fchange = 0.1; 
%E(f>fchange) = Exx(f>fchange) +Eyy(f>fchange) ; % use GPS for scalar energy of wind waves
%E(f<=fchange) = Ezz(f<=fchange); % use heave acceleratiosn for scalar energy of swell
%end

% testing bits
%E = nanmean([Ezz' (Exx+Eyy)'],2)';
%E = Eyy+Exx; % pure GPS version (for testing)
%E( check > maxEratio | check < minEratio ) = 0; 
%figure, loglog(f,check)
%clf, loglog(f,UU+VV,'g',f,Exx+Eyy,'b',f,AZAZ.*9.8^2,'m',f,Ezz,'r'),legend('UU+VV','XX+YY','AZAZ','ZZ') % for testing
%loglog(f,abs(Cxz),f,abs(Cyz))
%drawnow

%% wave stats
% fwaves = f>0.04 & f<1; % frequency cutoff for wave stats, 0.4 is specific to SWIFT hull
% 
% E( ~fwaves ) = 0;
fwaves = f>0;

% significant wave height
Hs  = 4*sqrt( sum( E(fwaves) ) * bandwidth);

%  energy period
fe = sum( f(fwaves).*E(fwaves) )./sum( E(fwaves) );
[~ , feindex] = min(abs(f-fe));
Ta = 1./fe;

% peak period
[~ , fpindex] = max(UU+VV); % can use velocity (picks out more distint peak)
%[~ , fpindex] = max(E);
Tp = 1./f(fpindex);


%% spectral directions
dir = - 180 ./ 3.14 * dir1;  % switch from rad to deg, and CCW to CW (negate)
dir = dir + 90;  % rotate from eastward = 0 to northward  = 0
dir( dir < 0 ) = dir( dir < 0 ) + 360;  % take NW quadrant from negative to 270-360 range
westdirs = dir > 180;
eastdirs = dir < 180;
dir( westdirs ) = dir ( westdirs ) - 180; % take reciprocal such wave direction is FROM, not TOWARDS
dir( eastdirs ) = dir ( eastdirs ) + 180; % take reciprocal such wave direction is FROM, not TOWARDS

% directional spread
spread = 180 ./ 3.14 .* spread1;


%% dominant direction

% or peak direction (very noisy)
%Dp = dir(fpindex); % dominant (peak) direction, use peak f

% or average
Dp = dir(fpindex); % dominant (peak) direction, use peak f

%% screen for bad direction estimate, or no heave data    

inds = fpindex + [-1:1]; % pick neighboring bands
if all(inds>0), 
    
  dirnoise = std( dir(inds) );

  if dirnoise > 45 
      Dp = 9999;
  else
      Dp = Dp;
  end
  
else
    Dp =9999;
end


%% prune high frequency results
E( f > maxf ) = [];
dir( f > maxf ) = [];
spread( f > maxf ) = [];
a1( f > maxf ) = [];
b1( f > maxf ) = [];
a2( f > maxf ) = [];
b2( f > maxf ) = [];
check( f > maxf ) = [];
f( f > maxf ) = [];


else % if not enough points or sufficent sampling rate or data, give 9999
  
     Hs = 9999;
     Tp = 9999; 
     Dp = 9999; 
     E = 9999; 
     f = 9999;
     a1 = 9999;
     b1 = 9999;
     a2 = 9999;
     b2 = 9999;
     check = 9999;

end


% quality control
if Tp>20,   
     Hs = 9999;
     Tp = 9999; 
     Dp = 9999; 
else 
end

