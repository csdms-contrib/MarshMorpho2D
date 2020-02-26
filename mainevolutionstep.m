function [IO, fIO, PLT] = mainevolutionstep(A, PARAMS, dx, dt, IO, fIO)

% input/output cell values
z = IO.z;
msl = IO.msl;
% input/output fluxes
FLX = fIO.FLX;
KBTOT = fIO.KBTOT;

% extract all parameters from the structure
names = fieldnames(PARAMS);
for i = 1:length(names)
    eval([names{i} '=PARAMS.' names{i} ';' ]);
end

% pull dimensions of the model grid
[N, M] = size(A);

% update Sea level
msl = msl + RSLR*dt;

% calculate Water depth
[h, ho, fTide, dtide, dHW, wl] = getwaterdepth(Trange, msl, z, kro);

% define Vegetation
lev = -z-msl; % water depth with respect of MSL
VEG = lev>dBlo; % depths above lower limit for veg growth
B = 4*(lev-dBup).*(dBlo-lev)/(dBlo-dBup)^2;
B(lev>dBup)=0;
B(lev<dBlo)=0;

% define Manning's n at grid cells 
MANN = 0*A+Cb;
MANN(VEG==1) = Cv;

% Tide & Surge flow
DHeff = min(Trange, dtide); % the water level excursion for the total water prism min
[U, Uy, Ux] = tidalFLOW(A, MANN, h, dHW, dx, DHeff, Ttide, kro);

%% MORPHODYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT-DRIVEN TRANSPORT (Tide) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total sediment resupension MUD
[E2] = totalsedimenterosionMUDsine(U, MANN, VEG, fTide, taucr, ...
                                   tcrgradeint, taucrVEG, me, ...
                                   h, lev, TrangeVEG);          
E2(A==2) = 0; % needed for b.c.
    
% ADiffusion Sediment transport
WS = A*0+ws2;
WS(VEG==1) = wsB;
[EmD2, SSM, FLX] = sedtran(h, A, DiffS, h, ho, E2, WS, dx, dt, rbulk2, ...
                           co2, Ux, Uy, FLX, fTide, Ttide, kro);     
SSC = SSM./h;
    
% end of compute tide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bed evolution erosion/depositon from tidal transport
z = z + dt*EmD2;

% Organic accretion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z - B*Korg*dt; 
KBTOT = KBTOT + sum(B(A==1))*Korg*dt;

%%%%%%%%%%%%%%%%%%%%% BED EVOLUTION DIVERGENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = bedcreep(z, A, crMUD, crMARSH, dx, dt, VEG); % MUD CREEP

%%%% IMPOSE "NEUMANN" boundary condition for morphodynamics %%%%%%%%%%%%%%%
% Translate the first boundary cell bed elevetion
if imposeseaboundarydepthmorphoALL==1
    [z] = seaboundaryNeumanbedelevationALLBOUNDARY(A, z);
end

%%%%%%%%%%%%%%%%% END OF MORPHODYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% UPDATE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IO.z = z;
IO.msl = msl;
fIO.FLX = FLX;
fIO.KBTOT = KBTOT;

% HERE YOU CAN ADD ANY VARIABLE TO BE SHOWN AS OUTPUT. 
% Just write PLT.name = name
PLT.U = U;
PLT.SSC = SSC;
PLT.EmD2 = EmD2;
PLT.Ux = Ux;
PLT.Uy = Uy;
PLT.VEG = VEG;
PLT.wsB = wsB;
PLT.h = h;
PLT.B = B;
PLT.fTide = fTide;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              