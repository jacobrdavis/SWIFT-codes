%% dynamic filtering cutoff
function [fc,fp] = dynamicFilterCutoff(z,fs,ratio,concatenateFrequencies)
    fc_bp = [0.05 0.25];
    [b,a] = butter(4,fc_bp/(fs/2),'bandpass');
    z = filtfilt(b,a,z);
    [E,f] = pwelch(z,[],[],[],fs);
    [~,locs] = findpeaks(E,f,'MinPeakHeight',0.9*max(E));
    fp = min(locs);
    fc = ratio*fp;

    figure;
    plot(f,E,'color','b'); hold on 
    findpeaks(E,f,'MinPeakHeight',0.9*max(E))
    yline(0.9*max(E))
    xline(fc,'--')
    xline(fp,'-')
    set(gca, 'XScale', 'log')

    if concatenateFrequencies == 1
        fc = [fc fp];
    end
end