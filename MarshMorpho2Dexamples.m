clear; close all; clc
% NOTES

% A
% 0: not in the domain
% 1: a normal cell
% 2: the open sea boundary

% Various initialization for plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cdata = [-1 205 165 0 0 0 0 0 0 255 255 255 1 0 0 0]; dlmwrite('mycmap.cpt', cdata, ' ');
cdata = [-1 205 165 0 0 0 0 0 0 255 255 255 1 0 0    0]; dlmwrite('mycmapCLR.cpt', cdata, ' ');
cdata = [-1 205 165 0 0 0 0 0]; dlmwrite('mycmap1.cpt', cdata, ' ');

%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
P = struct;
P.kro = 0.01; % minimum water depth [m]
P.DiffS = 0.1; % coefficient for tidal dispersion [-]

P.RSLR = 1/1000/365; % Relative Sea Level Rise from mm/yr to m/day (the time unit is the day!)

% tide
P.Ttide = 12.5/24; % tidal period [day]
P.Trange = 3.1; % 1.5; % tidal Trange [m]
P.TrangeVEG = 3.1; % 2.7; % 1.5; % tidal Trange [m] MOST OF THE TIME THIS IS THE SAME OF TIDAL RANGE

% SSC at the sea boundary
P.co2 = 5/1000; % Sea boundary SSC for mud [g/l]

P.ws2 = 0.2/1000; % mud settling velocity [m/s]
P.rbulk2 = 800; % dry bulk density [kg/m3]

% mud parameters
P.me = 10^-5*24*3600;  % mud erodability kg/m2/day!!!
P.tcrgradeint = 0.2; % linear increase in tcr below MLW [pa/m]
P.taucr = 0.2; % critical shear stress unvegetated [Pa]
P.crMUD = 3.65/365; % [m2/day] unvegetated creep coefficient
P.crMARSH = 0.2/365; % [m2/day] vegetated creep coefficient

% Vegetation parameters
P.dBlo = -(0.237*P.TrangeVEG-0.092); % lower limit for veg growth [m]  see McKee, K.L., Patrick, W.H., Jr., 1988. 
P.dBup = P.TrangeVEG/2; % upper limit for veg growth [m] ee McKee, K.L., Patrick, W.H., Jr., 1988. 
P.Cb = 0.02; % manning coefficient unvegegated
P.Cv = 0.1; % mannign coefficient in vegegated areas
P.wsB = 1/1000; % settling velocity in vegetated area [m/s]
P.taucrVEG = 0.5; % Critical stress for vegetated areas

P.Korg = 8/1000/365; % organic sediment production [m/day]

% to impose the boundary condition. 
% if = 1 the boundary cell becomes equal to the closest cell in the domain 
% if = 0 the boundary depth does not change
P.imposeseaboundarydepthmorphoALL = 1;

tmax = 5; % 500; % how many time steps in total
dtO = 2*365; % length of one time step
tINT = 1; % how many time step you should plot
time = [1:tmax]*dtO/365;

% to create a video (no or yes)
makevideo = 1;

%%%%%%%%%%%%%%% Geometry Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here are two examples. chose one or make your own grid
[N,M,dx,A,z,x,y,msl] = initializegeometrySWENEY();
%[N,M,dx,A,z,x,y,msl] = initializegeometryWEST(P);

IO.z = z;
IO.msl = msl;
fIO.FLX = 0; % keep track of sediment flux through open boundary
fIO.KBTOT = 0; % keep track of how much organic was produced

zoriginal = z; % store the original bed elevation
sumY2IN = sum(z(A==1)); % Store value for mass balance check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MAIN LOOP
if makevideo==1
    v = VideoWriter('PUTHEREVIDEONAME','Motion JPEG AVI');
    open(v);
end

s = 0;
for t = 1:tmax  

    % THIS IS THE CORE OF THE MODEL!!!!!!!!!
    [IO,fIO,PLT] = mainevolutionstep(A,P,dx,dtO,IO,fIO);

    %%% PLOT
    if mod(t,tINT)==0
        s=s+1;  

        % read the variables
        names = fieldnames(IO);
        for i = 1:length(names)
            eval([names{i} '=IO.' names{i} ';' ]);
        end

        % read the fluxes
        names = fieldnames(fIO);
        for i = 1:length(names)
            eval([names{i} '=fIO.' names{i} ';' ]);
        end

        % read the plot
        names = fieldnames(PLT);
        for i = 1:length(names)
            eval([names{i} '=PLT.' names{i} ';' ]);
        end

        ax1 = subplot(1,3,1); 
        IM = imagesc(y,x,-(z+msl));
        axis equal;
        set(IM,'alphadata',~(A==0));
        set(gca,'YDir','normal'); % colormap('jet');
        % cmp=demcmap([-3 P.Trange/2],256); %-3 1
        % colormap(ax1,cmp)
        caxis([-3 P.Trange/2]);
        title(strcat(num2str(time(t)),' years'))

        ax1 = subplot(1,3,2);
        IM=imagesc(y,x,-(zoriginal));
        axis equal;
        set(IM,'alphadata',~(A==0));
        set(gca,'YDir','normal'); % colormap('jet');
        % cmp=demcmap([-3 P.Trange/2],256);
        % colormap(ax1,cmp)
        caxis([-3 P.Trange/2]);
        title('original bathymetry')

        ax2 = subplot(1,3,3);
        IM=imagesc(y,x,1000*SSC);
        set(IM,'alphadata',~(A==0));
        axis equal;
        set(gca,'YDir','normal'); %colormap('jet');%colorbar('hori') 
        caxis([0 20]);
        colormap('jet')
        colormap(ax2,flipud(hot))
        title('SSC [mg/l]')

        sumY2 = sum(z(A==1));
        sumFLUX2 = -FLX*dx;
        % NOTE: Thsi is the equivalent volumetric flux, not the mass flux
        checksum = [-(sumY2IN-sumY2)+sumFLUX2/dx^2]+KBTOT;
        % if close to zero it means there are no issues (mass is conserved!!!)

        if makevideo==1
            V = getframe(figure(1));writeVideo(v,V);
        end
        pause(0.1)
    end

end

if makevideo==1
    close(v);
end % UNCOMMENT THIS TO CREATE A VIDEO