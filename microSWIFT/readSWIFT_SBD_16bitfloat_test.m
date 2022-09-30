% J. Davis
% Read 16bit floats with MATLABs fread

SWIFT = readSWIFT_SBD_16bitfloat('microSWIFT019_TX_12Sep2022_165147UTC.dat','microSWIFT021');

%%
figure
plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.energy,'k','LineWidth',1.5)
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('frequency (Hz)')
    ylabel('E (m^2/Hz)')
%     ylim([10^(-3) 10^(1)])
    xlim([10^(-2) 10^(0)])
    title([SWIFT.ID,' (',num2str(SWIFT.lat),',',num2str(SWIFT.lon),'): ', newline(), ... %' at ',datestr(SWIFT.time) 
        'Hsig = ',num2str(SWIFT.sigwaveheight)  ,'; ', ...
        'Tp = '  ,num2str(SWIFT.peakwaveperiod) ,'; ', ...
        'Dp = '  ,num2str(SWIFT.peakwavedirT)   , ...
         ])

% exportgraphics(gcf,'TX_16bitfloat_test_Westport_12Jul2021_scalar_matlab.png','Resolution',400)

figure
subplot(4,1,1)
    plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.a1,'k','LineWidth',1.5)
    ylabel('a1') 
    ylim([-1 1])
    set(gca,'XTick',[])
    yline(0,'b-.')
subplot(4,1,2)
    plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.b1,'k','LineWidth',1.5)
    ylabel('b1')
    ylim([-1 1])
    set(gca,'XTick',[])
    yline(0,'b-.')
subplot(4,1,3)
    plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.a2,'k','LineWidth',1.5)
    ylabel('a2')
    ylim([-1 1])
    set(gca,'XTick',[])
    yline(0,'b-.')
subplot(4,1,4)
    plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.b2,'k','LineWidth',1.5)
    ylabel('b2')
    ylim([-1 1])
    yline(0,'b-.')
    xlabel('frequency (Hz)')

% exportgraphics(gcf,'TX_16bitfloat_test_Westport_12Jul2021_dirmom_matlab.png','Resolution',400)
%%
function [SWIFT] = readSWIFT_SBD_16bitfloat(fname,ID)

SWIFT.time = [];
SWIFT.lat = [];
SWIFT.lon = [];
SWIFT.ID = ID;

fid = fopen(fname); % open file for reading

payloadtype = fread(fid,1,'uint8=>char');
type = fread(fid,1,'uint8');
port = fread(fid,1,'uint8');
size = fread(fid,1,'uint16');

if type == 52 && size > 0 % microSWIFT, size should be 327 bytes
        disp('reading microSWIFT 52')
        SWIFT.sigwaveheight      = half.typecast(fread(fid, 1,'*uint16')).double; % sig wave height
        SWIFT.peakwaveperiod     = half.typecast(fread(fid, 1,'*uint16')).double; % dominant period
        SWIFT.peakwavedirT       = half.typecast(fread(fid, 1,'*uint16')).double; % dominant wave direction
        SWIFT.wavespectra.energy = half.typecast(fread(fid,42,'*uint16')).double; % spectral energy density of sea surface elevation
        fmin                     = half.typecast(fread(fid, 1,'*uint16')).double; 
        fmax                     = half.typecast(fread(fid, 1,'*uint16')).double; 
        fstep                    = (fmax - fmin) / (length(SWIFT.wavespectra.energy)- 1);
        SWIFT.wavespectra.freq   = fmin:fstep:fmax; % frequency
        SWIFT.wavespectra.a1     = double(fread(fid,42,'*int8'))/100; % spectral moment
        SWIFT.wavespectra.b1     = double(fread(fid,42,'*int8'))/100; % spectral moment
        SWIFT.wavespectra.a2     = double(fread(fid,42,'*int8'))/100; % spectral moment
        SWIFT.wavespectra.b2     = double(fread(fid,42,'*int8'))/100; % spectral moment
        SWIFT.wavespectra.check  = double(fread(fid,42,'*uint8'))/10; % spectral check factor (should be unity)
        SWIFT.lat                = fread(fid, 1,'float'); % Latitude
        SWIFT.lon                = fread(fid, 1,'float'); % Longitude
        SWIFT.watertemp          = half.typecast(fread(fid, 1,'*uint16')).double; % water temp
        SWIFT.salinity           = half.typecast(fread(fid, 1,'*uint16')).double; % salinity
        BatteryVoltage           = half.typecast(fread(fid, 1,'*uint16')).double; % battery level
        epochTime                = fread(fid, 1,'float'); % epoch time
        asDatetime               = datetime(epochTime, 'ConvertFrom', 'posixtime', 'TimeZone','UTC');
        SWIFT.time               = datenum(asDatetime); % time at end of burst
end
fclose(fid);        
end