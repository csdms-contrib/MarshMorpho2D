function [E] = totalsedimenterosionMUDsine(U, MANN, VEG, fTide, taucr, ...
                                           tcrgradeint, taucrVEG, me, ...
                                           h, lev, TrangeVEG)
    % Total sediment erosion of mud?
    % Inputs:
        % U: water velocity at each cell [m/s]
        % MANN: Manning's coefficient at each cell
        % VEG: water depth above lower limit for veg growht with respest to
               % mean sea level
        % fTide: hydroperiod
        % taucr: critical shear stress unvegetated [Pa]
        % tcrgradeint: linear increase in critical shear stress below mean
                     % water level [Pa/m]
        % taucrVEG: critical shear stress vegetated [Pa]
        % me: mud erodability [kg/m2/day]
        % h: reference water depths
        % lev: water depth with respect of mean sea level
        % TrangeVEG: tidal Trange [m] (usually same as tidal range)
    % Output:
        % E: tide-averaged erosion flux                     
                            
    fUpeak = pi/2; % term for intra-tidal velocity equation
    taucro = U*0+taucr; % create array of critical shear stress values
    taucro(VEG==1) = taucrVEG; % where there is vegetation update critical shear stress

    % increase tcr with depth (asusming an existing vertical distribution. 
    % USE with caution, only ok for simulation of small marsh domain
    xi = -lev-TrangeVEG/2;
    xi(xi<0)=0;
    taucro = taucro+xi*tcrgradeint;

    % tidal current erosion
    ncyc = 10; % number of tidal cycles
    E = 0;
    for i = 0:ncyc
        Utide = U*fUpeak*sin(i/ncyc*pi/2); % intra-tidal velocity
        tauC = 1030*9.81*MANN.^2.*h.^(-1/3).*Utide.^2; % instantaneous bed shear stress
        E = E+1/(ncyc+1)*me.*(sqrt(1+(tauC./taucro).^2)-1); % tide-averaged erosion flux for a tidal cycle
    end

end

