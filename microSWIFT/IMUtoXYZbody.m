function [IMU,filterFunc] = IMUtoXYZbody(IMU,filterType,detrendSignals,fc)
plotFlag = true;

t  = IMU.time;
fs = IMU.samplingrate;
dt = fs^(-1);

% if ~exist('fc',"var") || isempty(fc)
%     fc = [];
% end
[fc,filterFunc] = filterTypeEvaluation(filterType,fs,fc);

IMU.cutoff = fc; fc = fc(1);

figure
plot( datetime(IMU.time,'ConvertFrom','datenum'),IMU.acc)

for i = 1:3

    a = IMU.acc(:,i) - 10.025;
    if detrendSignals == 1; a = detrend(a); end %a  = a - mean(a); % uncommented on 5/6/2022

 % % % % zero out nonzero initial and end conditions
%     ZeroCrossings = find(diff(sign(a)));
%     
%     figure; hold on
%     plot(IMU.time,a)
%    	plot(IMU.time(ZeroCrossings(1)),a(ZeroCrossings(1)),'x')
%     plot(IMU.time(ZeroCrossings(end)),a(ZeroCrossings(end)),'x')
%     
% 
%     a(1:ZeroCrossings(1)) = 0;
%     a(ZeroCrossings(end):end) = 0;
%     
%     figure; hold on
%     plot(IMU.time,a)
% % % %    

    a_f = filterFunc(a);
%     fc2 = 6;
%     [b2,a2] = butter(3,fc2/(fs/2),'low');
%     a_f = filtfilt(b2,a2,a)
    if plotFlag == 1 && i==1; plotsignal(a,a_f,t,fs,fc,'acc (m/s^2)'); end

    v = cumtrapz(a_f)*dt; % m/s
    if detrendSignals == 1; v = detrend(v,'omitnan'); end
    v_f = filterFunc(v);
    IMU.vel(:,i) = v_f;
    if plotFlag == 1 && i==1; plotsignal(v,v_f,t,fs,fc,'vel (m/s)'); end

    p = cumtrapz(v_f)*dt; % m
    if detrendSignals == 1; p = detrend(p,'omitnan'); end
    p_f = filterFunc(p);
    if ~strcmp(filterType,'no_filter')
        p_f(1:round(4*fs/fc),:) = 0; % round(RC./dt*20) remove first portion, which can have initial oscillations from filtering
    end
    IMU.pos(:,i) = p_f;
    if plotFlag == 1 && i==1;  plotsignal(p,p_f,t,fs,fc,'pos (m)'); end

end


end

%% plot signals
function [] = plotsignal(y,y_f,t,fs,fc,label)

[E,f]     = pwelch(y,[],[],[],fs);

if ~isempty(y_f)
    [E_f,f_f] = pwelch(y_f,[],[],[],fs);
end

t = datetime(t,'ConvertFrom','datenum');

figure
subplot(2,1,1)
plot(t,y,'color','b','DisplayName','original'); hold on
if ~isempty(y_f)
    plot(t,y_f,'color','r','DisplayName','filtered')
end
ylabel(label)
xlabel('time (s)')
legend()

subplot(2,1,2)
plot(f,E,'DisplayName','original'); hold on
if ~isempty(y_f)
    plot(f_f,E_f,'DisplayName','filtered'); hold on
end
xline(fc,'DisplayName','fc')
set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
xlim([.01,max(f)])
ylabel('energy density')
xlabel('frequency (Hz)')
legend()

end

%% RC filter
function a = RCfilter(b, fc, fs)
RC = (2*pi*fc)^(-1);
alpha = RC / (RC + 1./fs);
a = b;
    for ui = 2:length(b) % speed this up 
        a(ui) = alpha * a(ui-1) + alpha * ( b(ui) - b(ui-1) );
    end

end

%% FFT function
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

%% Evaluate filterType
function [fc,filterFunc] = filterTypeEvaluation(filterType,fs,fc)
switch filterType % TODO: PUT IN A FUNCTION
    case 'no_filter'
        if ~exist('fc',"var") || isempty(fc)
            fc = fs;
        end
        filterFunc = @(y) y;
    case 'butterworth_highpass'
        if ~exist('fc',"var") || isempty(fc)
            fc = 0.04;
        end
        [b,a] = butter(4,fc/(fs/2),'high');
        filterFunc = @(y) filtfilt(b,a,y);
%     case 'butterworth_bandpass'
%         if ~exist('fc',"var") || isempty(fc)
%             fc = [0.04 0.5];
%         end
%         [b,a] = butter(4,fc/(fs/2),'bandpass');
%         filterFunc = @(y) filtfilt(b,a,y);
    case 'butterworth_bandpass'
        if ~exist('fc',"var") || isempty(fc)
            fc = [0.04 1];
        end
        [b,a] = butter(2,fc/(fs/2),'bandpass');
        filterFunc = @(y) filtfilt(b,a,y);
    case 'RC'
        if ~exist('fc',"var") || isempty(fc)
            fc = 0.04;
        end
        filterFunc = @(y) RCfilter(y,fc,fs);
    case 'RC_dynamic'
        if ~exist('fc',"var") || isempty(fc)
            [fc] = dynamicFilterCutoff(IMU.acc(:,3),fs,0.70,true);
        end
        filterFunc = @(y) RCfilter(y,fc,fs);
    case 'butterworth_dynamic'
        if ~exist('fc',"var") || isempty(fc)
            [fc] = dynamicFilterCutoff(IMU.acc(:,3),fs,0.70,true);
        end
        [b,a] = butter(3,fc(1)/(fs/2),'high');
        filterFunc = @(y) filtfilt(b,a,y);
end
end

%% SCRAPS

% design butterworth filter:
% fc = 0.04;
% [b,a] = butter(4,fc/(fs/2),'high');
% fc2 = 0.5;
% [b2,a2] = butter(3,fc2/(fs/2),'low');
% bandpass filter...


% ax = IMU.acc(:,1);
% ay = IMU.acc(:,2);
% az = IMU.acc(:,3);

% ax = detrend(ax,3);
% ay = detrend(ay,3);
% az = detrend(az,3);
% 
% ax_f = filterFunc(ax);
% ay_f = filterFunc(ay);
% az_f = filterFunc(az);
% 
% plotsignal(az,az_f,t,fs,'acc (m/s^2)')
% 
% vx = cumtrapz(ax_f)*dt; % m/s
% vy = cumtrapz(ay_f)*dt; % m/s
% vz = cumtrapz(az_f)*dt; % m/s
% 
% vx = detrend(vx,3,'omitnan');
% vy = detrend(vy,3,'omitnan');
% vz = detrend(vz,3,'omitnan');
% 
% vx_f = filterFunc(vx);
% vy_f = filterFunc(vy);
% vz_f = filterFunc(vz);
% 
% plotsignal(vz,vz_f,t,fs,'vel (m/s)')
% 
% x = cumtrapz(vx_f)*dt;
% y = cumtrapz(vy_f)*dt;
% z = cumtrapz(vz_f)*dt;
% 
% x = detrend(x,3,'omitnan');
% y = detrend(y,3,'omitnan');
% z = detrend(z,3,'omitnan');
% 
% x_f = filterFunc(b,a,x);
% y_f = filterFunc(b,a,y);
% z_f = filterFunc(b,a,z);
% 
% plotsignal(z,z_f,t,fs,'pos (m)')
% 
% % remove first portion, which can have initial oscillations from filtering
% x_f(1:round(RC./dt*20),:) = 0;
% y_f(1:round(RC./dt*20),:) = 0;
% z_f(1:round(RC./dt*20),:) = 0;
% 
% plotsignal(z,z_f,t,fs,'pos (m)')
% 
% IMU.vel = [vx_f vy_f vz_f];
% IMU.pos = [x_f y_f z_f];




% [h_filt,f_filt] = freqz(b,a,ff2,fs);
% h_ma_filt = abs(h_filt);
% h_ph_filt = unwrap(angle(h_filt))*180/pi;

% subplot(2,1,1)
% subplot(2,1,2)
% yyaxis left
% semilogx(f_filt,mag2db(h_ma_filt),'DisplayName','Ma','LineWidth',2); hold on
% ylabel('Magnitude (dB)')
% xlim([min(f),max(f)])
% yyaxis right
% semilogx(f_filt,h_ph_filt,'DisplayName','Ph','LineWidth',2);
% xlabel('f(Hz)')
% ylabel('Phase (deg)')
% legend()
% xlim([.01,max(f)])

% figure;
% [f,Ma,~] = pltFFT(az,fs);
% semilogx(f,Ma,'color','b'); hold on 
% [f,Ma,~] = pltFFT(azf,fs);
% semilogx(f,Ma,'color','r'); hold on 
% xline(0.05,'LineWidth',1.5); xline(0.5,'LineWidth',1.5); 
% xline(1/(2*pi*RC),'LineWidth',1.5,'color','r')
% xlabel('f (Hz)')
% ylabel('FFT(az) (m/s^2)')
