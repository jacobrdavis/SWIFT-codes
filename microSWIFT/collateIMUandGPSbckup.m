function [collatedOut] = collateIMUandGPS(IMU,GPS,maxFs)

% IMU.time = IMU.time(1:30);
% IMU.acc = IMU.acc(1:30,:);
startCrop = max(IMU.time(1),GPS.time(1));
endCrop = min(IMU.time(end),GPS.time(end));

datestr(startCrop, 'YYYY-mm-DD hh:MM:ss.fff'); %

imuCrop = IMU.time > startCrop & IMU.time < endCrop;
gpsCrop = GPS.time > startCrop & GPS.time < endCrop;

IMU.time = IMU.time(imuCrop);
IMU.acc = IMU.acc(imuCrop,:);
IMU.mag = IMU.mag(imuCrop,:);
IMU.gyro = IMU.gyro(imuCrop,:);

GPS.time = GPS.time(gpsCrop);
GPS.lat = GPS.lat(gpsCrop);
GPS.lon = GPS.lon(gpsCrop);
GPS.sog = GPS.sog(gpsCrop);
GPS.u(end+1) = GPS.u(end); 
GPS.v(end+1) = GPS.v(end); 
GPS.z(end+1) = GPS.z(end); 
GPS.u = GPS.u(gpsCrop);
GPS.v = GPS.v(gpsCrop);
GPS.z = GPS.z(gpsCrop);

datestr(IMU.time(1), 'YYYY-mm-DD hh:MM:ss.fff') %
datestr(GPS.time(1), 'YYYY-mm-DD hh:MM:ss.fff')

datestr(IMU.time(end), 'YYYY-mm-DD hh:MM:ss.fff')
datestr(GPS.time(end), 'YYYY-mm-DD hh:MM:ss.fff')

fs = round(maxFs);
t0 = IMU.time(1);
tf = IMU.time(end);

masterTime = datenum(transpose(t0:seconds(fs^(-1)):tf+seconds(1-fs^(-1)))); % datestr(masterTime, 'YYYY-mm-DD hh:MM:ss.fff')
wholeSecs = datenum(t0:seconds(1):tf);

accSorted = NaN(length(masterTime),3);
magSorted = NaN(length(masterTime),3);
gyroSorted = NaN(length(masterTime),3);

masterIdx = 1;

for second = wholeSecs 
   
    inSec = IMU.time==second;
    tc = IMU.time(inSec);
    len   = length(tc);
   
    datestr(IMU.time(inSec), 'YYYY-mm-DD hh:MM:ss.fff')
    datestr(second, 'YYYY-mm-DD hh:MM:ss.fff')
    

    if len <= fs
        sortedIdx = masterIdx:masterIdx+len-1;
        unsortedIdx = find(inSec).';
        if len < fs
            ;
        end
    elseif len > fs 
        %TODO: push overfilled seconds to previous or next second,
        %currently it is being thrown out!
        tmp = find(inSec).';
        if isnan(accSorted(masterIdx-(len-fs),:))
            sortedIdx = masterIdx-(len-fs):masterIdx+fs-1;
            unsortedIdx = tmp;
        else
            sortedIdx = masterIdx:masterIdx+fs-1;
            unsortedIdx = tmp(1:end-(len-fs));
            
            t1 = datenum(transpose(seconds(0):seconds(len^(-1)):seconds(1-len^(-1))))
            t2 = datenum(transpose(seconds(0):seconds(fs^(-1)):seconds(1-fs^(-1))))
            Interp = interp1(t1,IMU.acc(tmp,1),t2)
            
            figure
            plot(t2,Interp); hold on
            plot(t1,IMU.acc(tmp,1))
            plot(t2,IMU.acc(unsortedIdx,1))
            datestr(masterTime(masterIdx:masterIdx+fs-1), 'YYYY-mm-DD hh:MM:ss.fff')
        end
        
        datestr(IMU.time(tmp), 'YYYY-mm-DD hh:MM:ss.fff')
    end
    
    accSorted(sortedIdx,:)  = IMU.acc(unsortedIdx,:);
    magSorted(sortedIdx,:)  = IMU.mag(unsortedIdx,:);
    gyroSorted(sortedIdx,:) = IMU.gyro(unsortedIdx,:);
    masterIdx = masterIdx+fs;
end


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


%%

accInterp  = zeros(size(accSorted));
magInterp  = zeros(size(magSorted));
gyroInterp = zeros(size(gyroSorted));

for i = 1:3
    accInterp(:,i) = interp1NaN(accSorted(:,i));
    magInterp(:,i) = interp1NaN(magSorted(:,i));
    gyroInterp(:,i) = interp1NaN(gyroSorted(:,i));
end

% 
% figure
% plot(masterTime,accInterp); hold on
% plot(masterTime,accSorted(:,:))

latInterp = interp1(GPS.time.',GPS.lat.',masterTime);
lonInterp = interp1(GPS.time.',GPS.lon.',masterTime);
sogInterp = interp1(GPS.time.',GPS.sog.',masterTime);
uInterp = interp1(GPS.time.',GPS.u.',masterTime);
vInterp = interp1(GPS.time.',GPS.v.',masterTime);
zInterp = interp1(GPS.time.',GPS.z.',masterTime);

nonNaN = ~isnan(accInterp(:,1)) & ~isnan(magInterp(:,1)) & ~isnan(gyroInterp(:,1)) & ...
~isnan(latInterp) & ~isnan(lonInterp) & ~isnan(sogInterp) & ...
~isnan(uInterp)   & ~isnan(vInterp)   & ~isnan(zInterp);



% figure
% plot(GPS.time,GPS.lat,'o'); hold on
% plot(masterTime,latInterp,'o')



collatedOut.IMU.time = masterTime(nonNaN);
collatedOut.IMU.acc = accInterp(nonNaN,:);
collatedOut.IMU.mag = magInterp(nonNaN,:);
collatedOut.IMU.gyro = gyroInterp(nonNaN,:);

collatedOut.GPS.time = masterTime(nonNaN);
collatedOut.GPS.lat = latInterp(nonNaN);
collatedOut.GPS.lon = lonInterp(nonNaN);
collatedOut.GPS.sog = sogInterp(nonNaN);
collatedOut.GPS.u   = uInterp(nonNaN);
collatedOut.GPS.v   = vInterp(nonNaN);
collatedOut.GPS.z   = zInterp(nonNaN);


end


function x = interp1NaN(x)
    nanx = isnan(x);
    t    = 1:numel(x);
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
end



%% scraps

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
