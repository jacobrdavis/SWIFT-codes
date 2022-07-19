%% compute_evaluation_metrics.m
% Script to compute evaluation metrics between microSWIFT and CDIP (or 
% similar) wave spectral estimates. 
% 
% J. Davis, 2022-02-05
% adapted from <explorerawmicroSWIFTdata.m script by J. Thomson, 10/2020
%
% Dependencies: 
%   Mapping Toolbox (deg2km)

function [metrics] = compute_evaluation_metrics(SWIFT,CDIP,metrics,outputDir,figureVisibility,printFigures)

metrics.benchmark = char(CDIP.id);
metrics.benchmarkTime = CDIP.time;
metrics.separation = separation_dist(SWIFT,CDIP);
titleStr = createTitleStr(SWIFT,CDIP,metrics.separation);
filenameStr = createFilenameStr(SWIFT,CDIP,metrics.method,metrics.referenceFrame,metrics.filterType);
metrics.bulk.sigwaveheight_err  = CDIP.sigwaveheight  - SWIFT.sigwaveheight;
metrics.bulk.peakwaveperiod_err = CDIP.peakwaveperiod - SWIFT.peakwaveperiod;
metrics.bulk.peakwavedirT_err   = CDIP.peakwavedirT   - SWIFT.peakwavedirT;
metrics.spectral                = spectral_metrics(SWIFT,CDIP,titleStr,filenameStr,outputDir,figureVisibility,printFigures);

spectral_comparison_plots(SWIFT,CDIP,metrics,titleStr,filenameStr,outputDir,figureVisibility,printFigures)

end

%% Collate and interpolate SWIFT spectral statistics onto CDIP
function [SWIFT,CDIP]  = collate_wavespectra(SWIFT,CDIP)
fields = fieldnames(SWIFT.wavespectra);
fields(contains(fields,'freq')) = [];
for fi = 1:length(fields)
    interpolatedField = interp1(SWIFT.wavespectra.freq,SWIFT.wavespectra.(fields{fi}),CDIP.wavespectra.freq);
    SWIFT.wavespectra.(fields{fi}) = interpolatedField(~isnan(interpolatedField));
    CDIP.wavespectra.(fields{fi}) = CDIP.wavespectra.(fields{fi})(~isnan(interpolatedField));
end
freq = CDIP.wavespectra.freq(~isnan(interpolatedField));
SWIFT.wavespectra.freq = freq;
CDIP.wavespectra.freq  = freq;
end

%% Create plot title string
function titleStr = createTitleStr(SWIFT,CDIP,separation)
titleStr = [char(SWIFT.id),': ',datestr(SWIFT.time),                             ...
    ' (',num2str(SWIFT.lat,'%.4f'),',',num2str(SWIFT.lon,'%.4f'),')',newline,    ...
    '      ',char(CDIP.id),': ',datestr(CDIP.time),' (',num2str(CDIP.lat,'%.4f'),...
    ',',num2str(CDIP.lon,'%.4f'),')',newline,                                    ...
    'separation: ',num2str(separation,'%.4f'),' km'];
end

%% Create figure filename string
function filenameStr = createFilenameStr(SWIFT,CDIP,methodStr,referenceFrame,filterType)
filenamedate = datestr(CDIP.time); filenamedate = strrep(strrep(filenamedate(1:end-3),':',''),' ','_');
filenameStr = join([methodStr,referenceFrame,filterType,SWIFT.id,CDIP.id,filenamedate],'_'); % [char(method),'_',char(SWIFT.id),'_',char(CDIP.id),'_',filenamedate]
end
%% Compute spectral metrics
function [spectralMetricsTable] = spectral_metrics(SWIFT,CDIP,titleStr,filenameStr,outputDir,figureVisibility,printFigures)
[SWIFT,CDIP] = collate_wavespectra(SWIFT,CDIP);
f = SWIFT.wavespectra.freq;
swell_logical = f>=0.04 & f<=0.1;
sea_logical = f>=0.1;

fields = fieldnames(SWIFT.wavespectra);
fields(contains(fields,'freq')) = [];
metricArray = nan(6,length(fields));

for fi = 1:length(fields)
    if ~strcmp(fields{fi},'freq')
        yhat = SWIFT.wavespectra.(fields{fi});
        y    = CDIP.wavespectra.(fields{fi});

        RMSE_total = root_mean_square_err(yhat,y);
        bias_total = measurement_bias(yhat,y);
        RMSE_swell = root_mean_square_err(yhat(swell_logical),y(swell_logical));
        bias_swell = measurement_bias(yhat(swell_logical),y(swell_logical));
        RMSE_sea = root_mean_square_err(yhat(sea_logical),y(sea_logical));
        bias_sea = measurement_bias(yhat(sea_logical),y(sea_logical));
        
        metricArray(:,fi) = [RMSE_total;RMSE_swell;RMSE_sea;bias_total;bias_swell;bias_sea];
    end
end
spectralMetricsTable = array2table(round(metricArray,2));
spectralMetricsTable.Properties.VariableNames = fields;
spectralMetricsTable.Properties.RowNames = {'RMSE_total','RMSE_swell','RMSE_sea','bias_total','bias_swell','bias_sea'};

spectral_evaluation_plots(SWIFT,CDIP,f,titleStr,filenameStr,outputDir,figureVisibility,printFigures)

end

    
%% Compute root mean square error
function [RMSE] = root_mean_square_err(yhat,y)
    err = y - yhat;
    RMSE = sqrt(mean(err.^2));
end

%% Compute bias
function [b] = measurement_bias(yhat,y)
    b = mean(yhat) - mean(y);
end

function [r] = separation_dist(SWIFT,CDIP)
    dlon = SWIFT.lon - CDIP.lon;
    dlat = SWIFT.lat - CDIP.lat;
    dy = deg2km(dlat);
    dx = deg2km( dlon, 6371 * cosd( SWIFT.lat) );
    r = sqrt( dx.^2 + dy.^2 );
end

%% Spectral evaluation plots
function [] = spectral_evaluation_plots(SWIFT,CDIP,f,titleStr,filenameStr,outputDir,figureVisibility,printFigures)

figure('visible',figureVisibility);
x0=6; y0=4; width=6.5;height=10;
set(gcf, 'Units', 'Inches', 'Position', [x0, y0, width, height], ...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
t = tiledlayout('flow','TileSpacing','tight','Padding','loose');
cmap = 'parula'; % getPyPlot_cMap('hot')
title(t,titleStr,'FontName','FixedWidth')
xlabel(t,CDIP.id)
ylabel(t,SWIFT.id)
nexttile(t); hold on
    scatter(log10(CDIP.wavespectra.energy),log10(SWIFT.wavespectra.energy),50,f,'filled'); hold on
    colormap(cmap) 
    caxis(round([min(f),max(f)],1))
    xl = [-4 2]; yl = [-4 2];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('log(E)')
nexttile(t); hold on
    scatter(CDIP.wavespectra.check,SWIFT.wavespectra.check,50,f,'filled'); hold on
    colormap(cmap) 
    caxis(round([min(f),max(f)],1))
    xl = [0 10]; yl = [0 10];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('check')
nexttile(t); hold on
    scatter(CDIP.wavespectra.a1,SWIFT.wavespectra.a1,50,f,'filled'); hold on
    colormap(cmap) 
    caxis(round([min(f),max(f)],1))
    xl = [-1 1]; yl = [-1 1];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('a1')
nexttile(t); hold on
    scatter(CDIP.wavespectra.b1,SWIFT.wavespectra.b1,50,f,'filled'); hold on
    colormap(cmap) 
    caxis(round([min(f),max(f)],1))
    xl = [-1 1]; yl = [-1 1];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('b1')
nexttile(t); hold on
    scatter(CDIP.wavespectra.a2,SWIFT.wavespectra.a2,50,f,'filled'); hold on
    colormap(cmap) 
    caxis(round([min(f),max(f)],1))
    xl = [-1 1]; yl = [-1 1];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('a2')
nexttile(t); hold on
    scatter(CDIP.wavespectra.b2,SWIFT.wavespectra.b2,50,f,'filled'); hold on
    colormap(cmap) % getPyPlot_cMap('hot')
    caxis(round([min(f),max(f)],1))
    xl = [-1 1]; yl = [-1 1];
    plot(xl,yl,'k') % 1:1 line
    ylim(xl); xlim(xl)
    set(gca,'FontSize',12)
    pbaspect([1 1 1])
    title('b2')
cb = colorbar;
cb.Label.String = 'frequency (Hz)';
cb.Label.FontSize = 12;
cb.Layout.Tile = 'east';
if printFigures == 1; print(gcf,[char(outputDir),'figures/',char(filenameStr),'_scatter','.png'],'-dpng'); end

end

%% Spectral comparison plots
function [] = spectral_comparison_plots(SWIFT,CDIP,metrics,titleStr,filenameStr,outputDir,figureVisibility,printFigures)

figure('visible',figureVisibility); hold on
    x0=6; y0=4; width=6.5;height=5;
    set(gcf, 'Units', 'Inches', 'Position', [x0, y0, width, height], ...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
    plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.energy,...
        'DisplayName','microSWIFT')
    plot(CDIP.wavespectra.freq,CDIP.wavespectra.energy,...
        'DisplayName','CDIP','color','k','LineWidth',1.5)
    if length(metrics.cutOffFrequency) > 1
        xline(metrics.cutOffFrequency(1),'--','DisplayName','cutoff')
        xline(metrics.cutOffFrequency(2),'-','DisplayName','peak')
    else
        xline(metrics.cutOffFrequency(1),'--','DisplayName','cutoff')
    end
    set(gca,'YScale','log','XScale','log')
    legend('Location','northeast')
    xlabel('frequency (Hz)')
    ylabel('E (m^2s)')
    title(titleStr,'FontName','FixedWidth','FontWeight','normal')

if printFigures==1; print(gcf,[char(outputDir),'figures/',char(filenameStr),'_scalar_spectra','.png'],'-dpng'); end
% f = gcf; exportgraphics(f,[outputDir,char(SWIFT.id),'_',char(CDIP.id),'_',filenamedate,'_scalar_spectra','.png'],'Resolution',300)

figure('visible',figureVisibility)
    x0=6; y0=4; width=6.5;height=7.5;
    set(gcf, 'Units', 'Inches', 'Position', [x0, y0, width, height], 'PaperUnits', 'Inches', 'PaperSize', [width, height])
    xl = [min(SWIFT.wavespectra.freq) max(SWIFT.wavespectra.freq)];
    yl = [-1 1];
    t = tiledlayout(5,1,'TileSpacing','tight','Padding','loose');
    title(t,titleStr,'FontName','FixedWidth')
    nexttile(1); hold on
        plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.a1,...
            'DisplayName','microSWIFT')
        plot(CDIP.wavespectra.freq,CDIP.wavespectra.a1,...
            'DisplayName','CDIP','color','k','LineWidth',1.5)
        xlim(xl); ylim(yl); set(gca,'XTick',[])
        ylabel('a1')
     nexttile(2); hold on
        plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.b1,...
            'DisplayName','XYZwaves')
        plot(CDIP.wavespectra.freq,CDIP.wavespectra.b1,...
            'DisplayName','CDIP','color','k','LineWidth',1.5)
        xlim(xl); ylim(yl); set(gca,'XTick',[])
        ylabel('b1')
     nexttile(3); hold on
        plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.a2,...
            'DisplayName','XYZwaves')
        plot(CDIP.wavespectra.freq,CDIP.wavespectra.a2,...
            'DisplayName','CDIP','color','k','LineWidth',1.5)
        xlim(xl); ylim(yl); set(gca,'XTick',[])
        ylabel('a2')
     nexttile(4); hold on
        plot(SWIFT.wavespectra.freq,SWIFT.wavespectra.b2,...
            'DisplayName','XYZwaves')
        plot(CDIP.wavespectra.freq,CDIP.wavespectra.b2,...
            'DisplayName','CDIP','color','k','LineWidth',1.5)
        xlim(xl); ylim(yl);
        ylabel('b2')
        xlabel('frequency (Hz)')

    lgd = legend;
    lgd.Layout.Tile = 5;

if printFigures==1; print(gcf,[char(outputDir),'figures/',char(filenameStr),'_directional_coeff','.png'],'-dpng'); end
% f = gcf; exportgraphics(f,[outputDir,char(SWIFT.id),'_',char(CDIP.id),'_',filenamedate,'_directional_coeff','.png'],'Resolution',300)

end