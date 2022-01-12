function [IMU] = sortIMU(IMU)
%% sortmicroSWIFT_IMU.m
% Function containing methods to sort microSWIFT IMU data onto a master
% clock for a single recording burst.
%
% Usage:
%   [IMU] = sortIMU(IMU);
%
% Inputs:
%   IMU - structure with fields:
%       time - UTC time (MATLAB datenum)
%       acc  - vector of accelerations [ax, ay, ax] (m/s^2)
%       mag  - vector of magnetometer readings [mx, my, mz] (uTesla)
%       gyro - vector of gyroscope readings [gx gy gz] (deg/s)
%
% Outputs: 
%   IMU - structure with fields sorted and interpolated onto a master clock
%
% J. Davis 2021-01-11

%% Sort IMU onto master clock

% Convert to datetime format with ms accuracy; prevents rounding problems
IMU.time = datetime(IMU.time,'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm:ss.SSS');

% Extract rounded sample rate and start-stop times
fs = round(IMU.samplingrate);
t0 = IMU.time(1);
tf = IMU.time(end);

masterTime = transpose(t0:seconds(fs^(-1)):tf+seconds(1-fs^(-1))); % datestr(masterTime, 'YYYY-mm-DD hh:MM:ss.fff') % datestr(masterTime(end), 'YYYY-mm-DD hh:MM:ss.fff') % masterTime = datenum(transpose(t0:seconds(fs^(-1)):tf+seconds(1-fs^(-1)))); % datestr(masterTime, 'YYYY-mm-DD hh:MM:ss.fff') % datestr(masterTime(end), 'YYYY-mm-DD hh:MM:ss.fff')
wholeSecs = t0:seconds(1):tf;

for second = wholeSecs
    inSec = IMU.time==second;
    len = sum(inSec);
    fs0 = seconds(len^(-1));
    secSteps = datenum(transpose(seconds(0):fs0:seconds(1)-fs0)); % top of second; better   %  secSteps = datenum(transpose(seconds(0)+fs0/2:fs0:seconds(1)-fs0/2)); % centered; not good
    IMU.time(inSec) = IMU.time(inSec) + secSteps;
end

%% Interpolate IMU onto master clock
IMU.acc  = interp1(IMU.time,IMU.acc,masterTime); % TODO: generalize to more fields?
IMU.mag  = interp1(IMU.time,IMU.mag,masterTime);
IMU.gyro = interp1(IMU.time,IMU.gyro,masterTime);

%% Overwrite time field with master time array
IMU.time = datenum(masterTime); % convert back to datenum?

%% Crop NaN values
nonNaN = ~isnan(IMU.acc(:,1))  & ~isnan(IMU.acc(:,2))  & ~isnan(IMU.acc(:,3)) & ...
         ~isnan(IMU.mag(:,1))  & ~isnan(IMU.mag(:,2))  & ~isnan(IMU.mag(:,3)) & ...
         ~isnan(IMU.gyro(:,1)) & ~isnan(IMU.gyro(:,2)) & ~isnan(IMU.gyro(:,3));

IMU = cropStructFields(IMU,nonNaN);

end

%% Subfunctions
function S = cropStructFields(S,cropIdx)
    fields = fieldnames(S);
    for f = 1:length(fields)
        if length(S.(fields{f})) == length(cropIdx)
            S.(fields{f}) = S.(fields{f})(cropIdx,:);
        end
    end
end

%%% scraps

%% wo for loop?

% https://www.mathworks.com/matlabcentral/answers/37196-count-number-of-specific-values-in-matrix
% https://www.mathworks.com/help/matlab/ref/repelem.html
% https://www.mathworks.com/help/matlab/ref/histcounts.html#d123e600422
% https://www.mathworks.com/matlabcentral/answers/490729-linspace-with-varying-increment
% https://www.mathworks.com/matlabcentral/answers/35811-convert-matrix-in-single-column

% uv = unique(IMU.time);
% cnts  = histc(IMU.time,uv);
% t0 = repelem(zeros(size(cnts)),cnts)
% for c = 1:length(cnts)
%     n = cnts(c)
%     stepsize = 1/(n);
%     A = (0:n-1).*stepsize;
%     t0
% end
% 
% 
% figure
% histogram(n)
% 
% t0 = repelem(zeros(size(n)),n)
%%


% 
% a = [ 1; 2; 3];
%  n = [10; 9; 10];
%  n= 10
%  stepsizes = (max(a,0)-min(a,0))/(n-1);
%  stepsizes = 1./(n-1);
%  stepsizes = 1./(n);
%  A = 0 + (0:n).*stepsizes;


% n = -2:0.5:0;
% x = 4.^n;
% 
%  a = [ 1; 2; 3];
%  n = 3;
%  stepsizes = (max(a,0)-min(a,0))/(n-1);
%  A = min(a,0) + (0:(n-1)).*stepsizes;
%  a = A(:);

%%

% accSorted = NaN(length(masterTime),3);
% magSorted = NaN(length(masterTime),3);
% gyroSorted = NaN(length(masterTime),3);
% 
% masterIdx = 1;



% for second = wholeSecs 
%    
%     inSec = IMU.time==second;
%     len = sum(inSec); 
% 
% %     tc = IMU.time(inSec);
% %     len   = length(tc);
%     
%     fs0 = seconds(len^(-1));
%     secSteps = datenum(transpose(seconds(0)+fs0/2:fs0:seconds(1)-fs0/2)); % centered; not good
%     secSteps = datenum(transpose(seconds(0):fs0:seconds(1)-fs0)); % top of second; better
%     
%     % Q: what happens if no readings are present in a second?
%     % A: it interpolates...
%     IMU.time(inSec) = IMU.time(inSec) + secSteps;
% % 
% %     if len == fs
% %         sortedIdx = masterIdx:masterIdx+len-1;
% %         unsortedIdx = find(inSec).';
% % 
% %         accSorted(sortedIdx,:)  = IMU.acc(unsortedIdx,:);
% %         magSorted(sortedIdx,:)  = IMU.mag(unsortedIdx,:);
% %         gyroSorted(sortedIdx,:) = IMU.gyro(unsortedIdx,:);
% % 
% %     elseif len < fs || len > fs 
% %         
% %         %TO DO: just center and interpolate later?
% %         % method 1: throwout:
% %         tmp = find(inSec).'; %unsorted?
% %         sortedIdx = masterIdx:masterIdx+fs-1;
% % %         unsortedIdx = tmp(1:end-(len-fs));
% % 
% %         t1 = datenum(transpose(seconds(0):seconds(len^(-1)):seconds(1-len^(-1))));
% %         t2 = datenum(transpose(seconds(0):seconds(fs^(-1)):seconds(1-fs^(-1))));
% %         
% %         accSorted(sortedIdx,:)  = interp1(t1,IMU.acc(tmp,:),t2);
% %         magSorted(sortedIdx,:)  = interp1(t1,IMU.mag(tmp,:),t2);
% %         gyroSorted(sortedIdx,:) = interp1(t1,IMU.gyro(tmp,:),t2);
% 
% %         if isnan(accSorted(masterIdx-(len-fs),:)) % fill previous
% %             sortedIdx = masterIdx-(len-fs):masterIdx+fs-1;
% %             unsortedIdx = tmp;
% %         else % downsample
% %             sortedIdx = masterIdx:masterIdx+fs-1;
% %             unsortedIdx = tmp(1:end-(len-fs));
% %             
% %             t1 = datenum(transpose(seconds(0):seconds(len^(-1)):seconds(1-len^(-1))))
% %             t2 = datenum(transpose(seconds(0):seconds(fs^(-1)):seconds(1-fs^(-1))))
% %             Interp = interp1(t1,IMU.acc(tmp,:),t2)
% %             
% %             figure
% %             plot(t2,Interp); hold on
% %             plot(t1,IMU.acc(tmp,:))
% %             plot(t2,IMU.acc(unsortedIdx,:))
% %             datestr(masterTime(masterIdx:masterIdx+fs-1), 'YYYY-mm-DD hh:MM:ss.fff')
% %         end
% %         
% %         datestr(IMU.time(tmp), 'YYYY-mm-DD hh:MM:ss.fff')
% %     end
%     
%     
% %     masterIdx = masterIdx+fs;
% end
% 
% 
% % for i = 1:3
% %     accInterp(:,i) = interp1(IMU.time,IMU.acc,masterTime);
% %     magInterp(:,i) = interp1NaN(magSorted(:,i));
% %     gyroInterp(:,i) = interp1NaN(gyroSorted(:,i));
% % end
% 
% % if len <= fs
% %         sortedIdx = masterIdx:masterIdx+len-1;
% %         unsortedIdx = find(inSec).';
% %     elseif len > fs 
% %         %TODO: push overfilled seconds to previous or next second,
% %         %currently it is being thrown out!
% %         
% %         %TO DO: just center and interpolate later?
% %         % method 1: throwout:
% %         tmp = find(inSec).';
% %         sortedIdx = masterIdx:masterIdx+fs-1;
% %         unsortedIdx = tmp(1:end-(len-fs));
% 
% %         if isnan(accSorted(masterIdx-(len-fs),:)) % fill previous
% %             sortedIdx = masterIdx-(len-fs):masterIdx+fs-1;
% %             unsortedIdx = tmp;
% %         else % downsample
% %             sortedIdx = masterIdx:masterIdx+fs-1;
% %             unsortedIdx = tmp(1:end-(len-fs));
% %             
% %             t1 = datenum(transpose(seconds(0):seconds(len^(-1)):seconds(1-len^(-1))))
% %             t2 = datenum(transpose(seconds(0):seconds(fs^(-1)):seconds(1-fs^(-1))))
% %             Interp = interp1(t1,IMU.acc(tmp,:),t2)
% %             
% %             figure
% %             plot(t2,Interp); hold on
% %             plot(t1,IMU.acc(tmp,:))
% %             plot(t2,IMU.acc(unsortedIdx,:))
% %             datestr(masterTime(masterIdx:masterIdx+fs-1), 'YYYY-mm-DD hh:MM:ss.fff')
% %         end
% %         
% %         datestr(IMU.time(tmp), 'YYYY-mm-DD hh:MM:ss.fff')
% 
% figure
% plot(IMU.time,IMU.acc,'-','LineWidth',1.5);  hold on
% plot(masterTime,accInterp,'--')
% legend({'x','y','z'})
% 
% figure
% plot(IMU.time,IMU.mag);  hold on
% plot(masterTime,magInterp)
% legend({'x','x','y','y','z','z'})
% 
% figure
% plot(IMU.time,IMU.gyro);  hold on
% plot(masterTime,gyroInterp)
% legend({'x','y','z','x','y','z'})
% 
% %
% figure
% plot(masterTime,accSorted)
% legend({'x','y','z'})
% 
% figure
% plot(masterTime,magSorted)
% legend({'x','y','z'})
% 
% figure
% plot(masterTime,gyroSorted)
% legend({'x','y','z'})
% 
% 

% 
% end




%%%%%%%%

% for second = wholeSecs 
%     inSec = IMU.time==second;
%     tc = IMU.time(inSec);
%     len   = length(tc);
%     
%     if len <= fs
%         sortedIdx = masterIdx:masterIdx+len-1;
%         unsortedIdx = inSec;
% 
% %         accSorted(sortedIdx,:) = IMU.acc(unsortedIdx,:);
% 
%     elseif len > fs 
%         %TODO: push overfilled seconds to previous or next second,
%         %currently it is being thrown out!
% 
%         tmp = find(inSec);
%         sortedIdx = masterIdx:masterIdx+fs-1;
%         unsortedIdx = tmp(1:end-1);
% %         acc = IMU.acc(inSec,:);
% %         unsortedIdx0 = 1:fs;
% %         accSorted(sortedIdx,:) = acc(unsortedIdx0,:);
% 
%         
%     end
%     accSorted(sortedIdx,:) = IMU.acc(unsortedIdx,:);
%     magSorted(sortedIdx,:) = IMU.mag(unsortedIdx,:);
%     gyroSorted(sortedIdx,:) = IMU.gyro(unsortedIdx,:);
%     masterIdx = masterIdx+fs;
% end


%         baseNaT(1:len) = tc;
%         ti = baseNaT + base;     
%         baseNaN(1:len) = IMU.acc(inSec,1);
%         IMU.time(IMU.time==second) = inSec + base;    
%         accSorted(masterIdx:len,:) = IMU.acc(inSec,:);


%         end    
%         inSec(len+1:fs) = NaT(fs-len,1);
%         +base

% baseNaT = NaT(fs,1);
% baseNaN = NaN(fs,1);
% base= seconds(0:fs^(-1):(1-fs^(-1))).';


% times = unique(IMU.time)

% 
% 
% 
% minFs=4;
% fs = minFs;
% t0 = GPS.time(1);
% tf = GPS.time(end);
% 
% masterTime = transpose(t0:seconds(fs^(-1)):tf+seconds(1-fs^(-1))); 
% datestr(GPS.time, 'YYYY-mm-DD hh:MM:ss.fff')
% wholeSecs = t0:seconds(1):tf;
% 
% latSorted = NaN(length(masterTime),1);
% lonSorted = NaN(length(masterTime),1);
% sogSorted = NaN(length(masterTime),1);
% u = NaN(length(masterTime),1);
% v = NaN(length(masterTime),1);
% z = NaN(length(masterTime),1);
% 
% masterIdx = 1;
% 
% for second = wholeSecs 
%     inSec = IMU.time==second;
%     tc = IMU.time(inSec);
%     len   = length(tc);
%     
%     if len <= fs
%         sortedIdx = masterIdx:masterIdx+len-1;
%         unsortedIdx = inSec;
% 
%     elseif len > fs 
%         %TODO: push overfilled seconds to previous or next second,
%         %currently it is being thrown out!
%         tmp = find(inSec);
%         sortedIdx = masterIdx:masterIdx+fs-1;
%         unsortedIdx = tmp(1:end-1);
%         
%     end
%     
%     lat(sortedIdx,:) = GPS.lat(unsortedIdx,:);
%     lon(sortedIdx,:) = GPS.lon(unsortedIdx,:);
%     sog(sortedIdx,:) = GPS.sog(unsortedIdx,:);
%     u(sortedIdx,:) = GPS.u(unsortedIdx,:);
%     v(sortedIdx,:) = GPS.v(unsortedIdx,:);
%     z(sortedIdx,:) = GPS.z(unsortedIdx,:);
% 
%     masterIdx = masterIdx+fs;
% end
% 

% %% Interpolate IMU NaNs
% 
% accInterp  = zeros(size(accSorted));
% magInterp  = zeros(size(magSorted));
% gyroInterp = zeros(size(gyroSorted));
% 
% for i = 1:3
%     accInterp(:,i) = interp1NaN(accSorted(:,i));
%     magInterp(:,i) = interp1NaN(magSorted(:,i));
%     gyroInterp(:,i) = interp1NaN(gyroSorted(:,i));
% end

% figure
% plot(masterTime,latSorted)
% 
% figure
% plot(masterTime,lonSorted)
% 
% figure
% plot(masterTime,sogSorted)
% 
% figure
% plot(masterTime,uSorted); hold on
% plot(masterTime,vSorted); 
% legend({'u','v'})
