% directiontest.m - function to evaluate direct comparison of gps waves directions with spotter 
% J. Davis
% created on: 2021-10-18
% updates:
%          
% J. Davis, 2021-12-23, cleaned up and stored in Spotter utility folder in
%                       SWIFT-codes directory


function [] = directiontest(spotter,i)
spot_spread = spotter(i).directionalSpread;
spot_dir = spotter(i).direction;
spot_f = spotter(i).frequency;
a1 = spotter(i).a1;
b1 = spotter(i).b1;
a2 = spotter(i).a2;
b2 = spotter(i).b2;

% GPSwaves routine:

% wave directions
dir1 = atan2(b1,a1) ;  % [rad], 4 quadrant
dir2 = atan2(b2,a2)/2 ; % [rad], only 2 quadrant
spread1 = sqrt( 2 * ( 1 - sqrt(a1.^2 + b2.^2) ) );
spread2 = sqrt( abs( 0.5 - 0.5 .* ( a2.*cos(2.*dir2) + b2.*cos(2.*dir2) )  ));

% spectral directions
dir = - 180 ./ 3.14 * dir1;  % switch from rad to deg, and CCz to Cz (negate)
dir = dir + 90;  % rotate from eastward = 0 to northward  = 0
dir( dir < 0 ) = dir( dir < 0 ) + 360;  % take Nz quadrant from negative to 270-360 range
westdirs = dir > 180;
eastdirs = dir < 180;
dir( westdirs ) = dir ( westdirs ) - 180; % take reciprocal such wave direction is FROM, not TOWARDS
dir( eastdirs ) = dir ( eastdirs ) + 180; % take reciprocal such wave direction is FROM, not TOWARDS

% directional spread
spread = 180 ./ 3.14 .* spread1;
%% plot of results

figure
ax1 = subplot(1,2,1); 
    scatter(spot_f,dir,15,...
        'markerFaceColor','k','markerEdgeColor','k','marker','o','markerFaceAlpha',0.8,'markerEdgeAlpha',0.8,...
        'DisplayName','GPSWaves'); hold on
    scatter(spot_f,spot_dir,15,...
        'markerFaceColor','r','markerEdgeColor','r','marker','o','markerFaceAlpha',0.8,'markerEdgeAlpha',0.8,...
        'DisplayName','Spotter on-board')
    xlabel('frequency')
    ylabel('direction')
    legend()
    pbaspect(ax1,[1 1 1])
    
ax2 = subplot(1,2,2); 
    scatter(dir,spot_dir,...
            'markerFaceColor','b','marker','o','markerFaceAlpha',0.5,'markerEdgeAlpha',0.5); hold on
    plot([0 100],[0 100],'k')
    % ylim([30 90]); xlim([30 90])
    xlabel('GPSWaves')
    ylabel('Spotter on-board')
    pbaspect(ax2,[1 1 1])

%

figure
ax1 = subplot(1,2,1);
    scatter(spot_f,spread,15,...
        'markerFaceColor','k','markerEdgeColor','k','marker','o','markerFaceAlpha',0.8,'markerEdgeAlpha',0.8,...
        'DisplayName','GPSWaves'); hold on
    scatter(spot_f,spot_spread,15,...
        'markerFaceColor','r','markerEdgeColor','r','marker','o','markerFaceAlpha',0.8,'markerEdgeAlpha',0.8,...
        'DisplayName','Spotter on-board')
    ylim([30 90])
    xlabel('frequency')
    ylabel('spread')
    legend()
    pbaspect(ax1,[1 1 1])
ax2 = subplot(1,2,2); 
    scatter(spread,spot_spread,...
            'markerFaceColor','b','marker','o','markerFaceAlpha',0.5,'markerEdgeAlpha',0.5); hold on
    plot([0 100],[0 100],'k')
    ylim([30 90]); xlim([30 90])
    xlabel('GPSWaves')
    ylabel('Spotter on-board')
    pbaspect(ax2,[1 1 1])