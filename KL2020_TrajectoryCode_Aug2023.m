%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Hail trajectory code
%
%  Updated: August 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%LOAD processed CM1 output file
load 3D_composites_WKsndumax31.mat;  %3D gridded CM1 output; need specific variables with names comp_rhoa, comp_u, etc. --> ask about!
compgridsize = size(comp_qc);

%======================================================================
% Initial Conditions -- Initialize Embryos in Model based on CM1 input grid
% INDEX
init_indx = (20:70);     %Indices of initial trajectory x position 
init_indy1 =(20:70);   %indices of initial trajectory y location poisitions. [20:70]
init_indz1 = (14:43);   %indices of initial trajectory z location [14:43]
nz = numel(init_indz1); %First height level 14 is ~5 km
    %----Microphysical Information------
init_diam =5.0;  %Diameter of initial particle (mm);
rho_init = 917;  %initial embryo density
init_m = pi ./ 6 .* rho_init .* (init_diam.*0.001).^ 3;  %initial mass of embryo (kg)
CCN = 250;  %number concentration of drops per CC from Morrison (250)
%=======================================================================

%======OTHER CONSTANTS AND DEFINITIONS ==============================
%Composite grid in physical units (~50 km x 50 xkm x 20 km)
xgrid = (-25:0.5:25); ygrid = xgrid;
zgrid=(0.125:0.25:19.875);
total_t = 10000;            %Total integration time in s 
delt = 1;                   %time step
nt = floor(total_t/delt);   %number of time steps
nxi = numel(init_indx);     %number of x points for seeds
nyi = numel(init_indy1);    %number of y points for seeds
np = numel(init_indx);      %number of seeds for each x index (= number of y indices)
part_max_si = zeros(np,1);  %Array of maximum sizes
tstgr = zeros(np,nyi,nz);   %allocate array for the max size grid.
max_size_grid = zeros(np,nyi,nz);   %allocate array for the max size grid.
final_posx = max_size_grid;         %allocate array for final x positions as function of starting position x,y,z
final_posy = max_size_grid;         %allocate array for final y positions as function of starting position x,y,z
final_posz = max_size_grid;         %allocate array for final z positions as function of starting position x,y,z
total_traj_time = max_size_grid;    %allocate array for total time of trajectory
%====================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define some thermodynamic constants        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = 9.81;           %gravitational acceleration in m.s2
cp=1004.0;          %specific heat capacity of air at constant pressure (J/kg/K)
cpi = 2108.0;       %specific heat capacity of ice at constant pressure (J/kg/K) -- note function of T, could improve
cpw = 4187.0;       %specific heat capacity of liquid at constant pressure (J/kg/K) -- note function of T
Rv = 461.5;         %gas constant for vapor (J/kg/K)
Rdry = 287.1;         %gas constant for dry air (J/kg/K)
p0 = 100000.0;      %reference surface pressure (Pa);
es0 = 611;          %equilibrium vapor pressure at T=T0=273.15 K (Pa)
T0=273.15;          %reference temperature (K)
lv = 2.5e6;         %enthalpy of vaporization (J/kg)    **NOTE THESE ARE FUNCTIONS OF T, COULD
lf= 3.33e5;         %enthalpy of fusion/melting (J/kg)  **ADD EXPRESSIONS TO MAKE MORE ACCURATE
ls = lv+lf;         %enthalpy of sublimation (J/kg)     *****
sfcdens = 1.1222;   %surface density (from composite, inflow) in kg/m3
rhol = 1000.0;      %density of liquid water kg/m3
rhosolid = 917.0;   %density of solid ice kg/m3
rhomin = 500.;      %minimum density allowed for rime, in kg/m3  ***need to udpate with Prodi data
Cdrag = 0.5;        %Drag coefficient for hailstones. Recent work suggests 0.8-1.0

    %Compute some important thermodynamic quantities aross the whole domain
R = (1-comp_qv).*Rdry + comp_qv.*Rv;    %gas constant for moist air
p = (comp_theta .* comp_rhoa .* R ./ p0.^(R./cp)).^(1./(1-R./cp)); %pressure (Pa)
T = (1./(comp_rhoa .* R)) .* (comp_theta .* comp_rhoa .* R ./ (p0.^(R./cp))).^(1./(1-R./cp)); %temperature (K)
esT = es0 .* exp(lv./Rv .* (1./T0-1./T));   %equilibrium vapor pressure as function of T (Pa)
rhos=esT./Rv./T;                            %equilibrium vapor density (kg/m3)
rhos0 = es0 ./Rv./T0;                       %reference equilibrium vapor density for T=T0 (kg/m3)
rhov = comp_qv .* comp_rhoa;                %Vapor density (kg/m3)
kT = (2.381 + 0.0071 .* (T-T0)) .* 1e-2;                %thermal conductivity (J/m/s/K)
Dv = (2.11e-5) .* (T./T0).^1.94 .* (p0./p);             %difffusivity of vapor in air (m2/s)
eta_a = (0.379565 + 0.0049.*T).*1e-5;                   %1.7e-5;  %dynamic viscosity (kg/m/s)
Kair = 9.1018e-11 .* T.^2 + 8.8197e-8 .* T - 1.0654e-5; %thermal diffusivity of air (m2/s)
nu = eta_a./comp_rhoa;                                  %kinematic viscosity

%Precompute some arrays for speed
guessT = 253.15;
Tarray = (guessT-30.0:0.05:guessT+20.0);  %make an array of previous T+/- 20 deg for minimization
rhossfc1= es0 .* exp(ls./Rv .* (1./T0-1./Tarray))./Rv./Tarray; %Find the rho_v_sfc that matches T from Tarray

Ecc = zeros(compgridsize(1),compgridsize(2),compgridsize(3)) + 1.0;   %Initialize efficiency everywhere equal to one 
Ecr = zeros(compgridsize(1),compgridsize(2),compgridsize(3)) + 1.0;   %Initialize efficiency everywhere equal to one 

    %Get Cloud and Rain Liquid Water Contents
LWC_cloud = comp_qc .* comp_rhoa;            %kg/m3 of cloud mass content
LWC_rain = comp_qr .* comp_rhoa;             %kg/m3 of rain mass content
IWC_total = (comp_qi + comp_qs).*comp_rhoa;  %kg/m3 of IWC (cloud and snow ice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hailstone Fall speed parameterization -- right now from MORR   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah = sqrt(4.0/3.0 .* rho_init ./ Cdrag .* g); %Fall speed coeff = 114.5 in Morr code
bh = 0.5;  %Fall speed exponent
vIMP = 0.65;   %impact velocity reduction following Rasmussen and Heymsfield 1985.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cloud Droplets and Collection Parameterization           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCC = CCN * 1e6; %convert number concentration to 1/m3
D_c = 1e6 .* (6.0 .* comp_rhoa .* comp_qc ./ NCC ./ pi ./ rhol).^(1.0./3.0);  %mean droplet size in micrometers
Ecthresh = 0.1;  %Threshold low value below which collection efficiency linearly drops off to zero
Dcthresh = 5.0;  %Threshold droplet diameter (microns) below which collection efficiency linearly drops off to zero
temp = find(D_c < Dcthresh);
Ecc(temp) = (Ecthresh./Dcthresh).*D_c(temp);   %Collection Efficiency for cloud droplets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rain Drops and Collection Parameterization               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE COMPUTE MASS WEIGHTED FALL SPEED OF RAINDROPS EVERYWHERE
velrain_Q = zeros(compgridsize(1),compgridsize(2),compgridsize(3));  %allocate array for mass-weighted rain velocity
deldr = 0.1;   %increment of raindrop sizes in mm
drmm = (0.1:deldr:8.0);   %raindrop diameters in mm
velrain = -0.102 + 4.932.*drmm - 0.9551.*drmm.^2 + 0.07934.*drmm.^3 - 0.002362.*drmm.^4; %Drop fall velocity (Brandes et al. 2002)
lamr = (0.001) .* (pi .* rhol .* comp_nr ./ comp_qr).^(1.0./3.0);  %find slope parameter for raindrops (1/mm) - leading coef is conversion to 1/mm
n0r = comp_nr .* lamr .* comp_rhoa;      %find intercept parameter for raindrop size distribution. (1/m3 mm)

        for xx=1:compgridsize(1)
            for yy=1:compgridsize(2)
                for zz=1:compgridsize(3)
                   rainND = n0r(xx,yy,zz) .* exp(-lamr(xx,yy,zz) .* drmm);  %rain DSD
                   velrain_Q(xx,yy,zz) = sum(rainND .* velrain .* drmm.^3 .* deldr) ./ sum(rainND .* drmm.^3 .* deldr);
                end
            end
        end
        
        %Clean up velocities
        cleanup = find(isnan(velrain_Q) == 1);
        velrain_Q(cleanup) = 0;
        
        %Set the collection efficiency for rain at 0.8 for now.
        Ecr(:) = 0.8;

%====================================================================%
%========= MAIN LOOPING STRUCTURE BEGINS HERE========================%
%====================================================================%
%%%% Loop over each height level; index = zz %%%%%%
for zz=1:nz
    initz = init_indz1(zz)
    init_indz = init_indx.*0 + initz;  %shaping arrays the right way

    %%%% LOOP OVER Y POSITIONS, index = yi %%%%%%%%
        for yi = 1:nyi  %begin loop over y positions
            init_indy = init_indx.*0 + init_indy1(yi);   

            for pp=1:np  %%%%%LOOP OVER PARTICLES, index = pp %%%%%
                
                %Allocate Arrays/Clear out each time (new particles)
                part_posx = zeros(nt+1,1);      %Particle x at each step
                part_posy = part_posx; part_posz = part_posx;
                part_si = zeros(nt+1,1);  %Particle size at each step
                part_vel = zeros(nt+1,1); %Particle fall speed
                GROWTH = zeros(nt+1,1); %Particle growth rate in m/s
                GROWTH_MICRO = zeros(nt+1,1); %Particle growth rate in m/s
                part_Ts = zeros(nt+1,1);  %Particle surface temperature K
                rhoh = zeros(nt+1,1); %Particle density (instantaneous)
                avg_rhoh = zeros(nt+1,1); %Particle average density (time integrated)
                FF = zeros(nt+1,1);  %frozen fraction trcace
                diffgrow = FF;   %diffusion growth
                massgrow = FF;   %mass accretion growth
                accretgrow =FF;  %accretion growth mass (redundant)
                BUGGYSMALLZ = FF;  %temp thing for debugging
                available_mass = 0;  %declare it up here for debugging
                needed_to_soak = 0;   %declare it up here for debugging
                soaked = 0;       %initially no soaked liquid
                remaining_unfrozen_mass = 0; %declare it up here for debugging
                cur_ind = zeros(3,1);  %Current particle position index
                record_Theta = zeros(nt+1,1);
                record_liqr = record_Theta;
                record_liqc = record_Theta;
                record_ice = record_liqr;
                record_qv = record_liqr;
                record_rhoa = record_liqr;
                record_u = record_liqr;
                record_v = record_liqr;
                record_w = record_liqr;
                record_T = record_liqr;
                record_sfcliq = record_T;
                record_soaked = record_T;
                
                %------------INITIAL CONDITIONS ----------------
                part_posx(1) = xgrid(init_indx(pp));    %initial particle X position (km)
                part_posy(1) = ygrid(init_indy(pp));    %initial particle Y position (km)
                part_posz(1) = zgrid(init_indz(pp));    %initial particle Z position (km)
                
                %Initialize particle fall speed. Based on aD^b equation
                %from Morrison 2-moment scheme. D in expression is in
                %meters, hence conversion factor of 1000. Accounts for
                %density effect.
                part_vel(1) = (ah .* (init_diam./1000).^bh).*sqrt(sfcdens./comp_rhoa(init_indx(pp),init_indy(pp),init_indz(pp))); 
                part_si(1) = init_diam;                 %particle initial size in mm
                tot_m = init_m;                         %initial particle mass is embryo mass (kg)
                cur_ind(1) = init_indx(pp);             %particle current location index (X)
                cur_ind(2) = init_indy(pp);             %particle current location index (Y)
                cur_ind(3) = init_indz(pp);             %particle current location index (Z)
                part_Ts(1) = T(init_indx(pp),init_indy(pp),init_indz(pp)); %initial particle temperature is ambient at initial location
                rhoh(1) = rho_init;                        %particle initial density kg/m3
                GROWTH_MICRO(1) = 0;                          %no growth on first time step
                tt=1;                                   %time begins at index 1
                
                %Calculate growth and position of growing hailstone as long
                %as the hailstone is above the ground:
                while part_posz(tt) > 0.0 && tt < nt   %KELLY bugfix- added max allowable time of integration so dont get infinite trajectory 
                      if T(cur_ind(1),cur_ind(2),cur_ind(3)) >273.15   %No growth below melting layer
                          GROWTH_MICRO(tt) = 0;
                          mass_accret = 0;
                          ice_accret = 0;
                          part_Ts(tt) = 273.15;% 
                          FrozenFrac = 0.0;
                          FF(tt) = 0.0;   %no freezing
                          massgrow(tt) = 0.0;   %%ADDED 4 april 2021
                          which_l = lv;   
                          esTs = es0 .* exp(which_l./Rv .* (1./T0-1./part_Ts(tt)));   %equilibrium vapor pressure as function of Ts (Pa)
                          rhossfc=esTs./Rv./part_Ts(tt); 
                          rhoh(tt) = 900.0;
                          N_Re = part_vel(tt) .* (0.001 .* part_si(tt)) ./ nu(cur_ind(1),cur_ind(2),cur_ind(3));   %Reynolds number
                          N_Pr = nu(cur_ind(1),cur_ind(2),cur_ind(3)) ./ Kair(cur_ind(1),cur_ind(2),cur_ind(3));   %Prandtl Number
                          N_Sc = nu(cur_ind(1),cur_ind(2),cur_ind(3)) ./ Dv(cur_ind(1),cur_ind(2),cur_ind(3));     %Schmidt Number
                            %Ventilation coefficient
                        ventH = 0.78 + 0.308 .* (N_Pr.^(1.0/3.0)) .* sqrt(N_Re); %ventilation coefficient for thermal energy
                        ventV = 0.78 + 0.308 .* (N_Sc.^(1.0/3.0)) .* sqrt(N_Re); %ventilation coefficient for vapor 
                      else             

%%%%%%%%%%%%% MICROPHYSICS %%%%%%%%%%%%%%%%%%%

                            %Some unitless numbers
                        N_Re = part_vel(tt) .* (0.001 .* part_si(tt)) ./ nu(cur_ind(1),cur_ind(2),cur_ind(3));   %Reynolds number
                        N_Pr = nu(cur_ind(1),cur_ind(2),cur_ind(3)) ./ Kair(cur_ind(1),cur_ind(2),cur_ind(3));   %Prandtl Number
                        N_Sc = nu(cur_ind(1),cur_ind(2),cur_ind(3)) ./ Dv(cur_ind(1),cur_ind(2),cur_ind(3));     %Schmidt Number
                            %Ventilation coefficients from P/R
                        ventH = 0.78 + 0.308 .* (N_Pr.^(1.0/3.0)) .* sqrt(N_Re); %ventilation coefficient for thermal energy
                        ventV = 0.78 + 0.308 .* (N_Sc.^(1.0/3.0)) .* sqrt(N_Re); %ventilation coefficient for vapor 
                        FrozenFrac = 1;   %added in to reset to 1 in case it got stuck
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % MASS GROWTH FROM ACCRETION  dm/dt|acc%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            Kerc = pi ./ 4 .* Ecc(cur_ind(1),cur_ind(2),cur_ind(3)) .* (0.001 .* part_si(tt)).^2 .* part_vel(tt);  %Collection Kernel assuming continuous collection of LWC
                            Kerr = pi ./ 4 .* Ecr(cur_ind(1),cur_ind(2),cur_ind(3)) .* (0.001 .* part_si(tt)).^2 .* ...
                                abs(part_vel(tt)-velrain_Q(cur_ind(1),cur_ind(2),cur_ind(3)).*...
                                sqrt(sfcdens./comp_rhoa(cur_ind(1),cur_ind(2),cur_ind(3))));     %collection kernel for raindrops, uses mass weighted fallspeed 
                                                                                                 %and density correction for rain (hail density correction done above). Ignores rain sizes for now.    
                            mass_accret = Kerc .* LWC_cloud(cur_ind(1),cur_ind(2),cur_ind(3)) +...
                                Kerr .* LWC_rain(cur_ind(1),cur_ind(2),cur_ind(3));  %mass accreted from cloud and rain liquid water kg/s
                        
                            ice_accret = 0;  %Ice accreted from cloud and snow ---> 0 Unless wet growth
                                                    
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % COMPUTE TEMPERATURE of Hailstone Surface Ts to see if there is Wet Growth     %     
                            %
                            % The following takes the T_sfc equation and subtracts it from Tarray. 
                            % make_it_hail should be minimal where the new T_sfc is. Nelson iterated,
                            % Here we just find the zero root                                               %
                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                                        
                         %USING RH87 heat coefficients
                           if N_Re <= 6000
                               make_it_hail = abs(Tarray - ...
                                (2.0.*pi.*(0.001.*part_si(tt)).*ventH.*kT(cur_ind(1),cur_ind(2),cur_ind(3)).*T(cur_ind(1),cur_ind(2),cur_ind(3))-...  % Conduction of thermal energy
                                2.0 .* pi .* (0.001.*part_si(tt)) .* ventV .* ls .* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc1-rhov(cur_ind(1),cur_ind(2),cur_ind(3))) + ...% away from hailstone term
                                ice_accret .* cpi .* T(cur_ind(1),cur_ind(2),cur_ind(3)) + ...         %term associated with heating up collected ice
                                mass_accret .* (lf + cpw.*T(cur_ind(1),cur_ind(2),cur_ind(3))) + ...   %terms associated with freezing liquid and heating up accreted liquid
                                available_mass./delt .* lf)./... + soaked.*lf./delt)./...                              %Freezing the remaining unfrozen mass on outside
                                (2.0.*pi.*(0.001.*part_si(tt)).*ventH.*kT(cur_ind(1),cur_ind(2),cur_ind(3)) + ...  %thermal conduction
                                ice_accret .* cpi + mass_accret .* cpw));
                           
                            %diffusion growth coefficient defined here based on Reynolds # regime
                                  diffcoef = 2.0 .* ventV;
                                                              
                           elseif N_Re >= 6000 && N_Re < 2e4
                               chiX = 0.76;  %Note: for spheres only. Increases with increasing oblatness 
                               make_it_hail = abs(Tarray - ...
                                    (chiX.*pi.*(0.001.*part_si(tt)).*sqrt(N_Re).*N_Pr.^(1.0./3.0).*kT(cur_ind(1),cur_ind(2),cur_ind(3)).*T(cur_ind(1),cur_ind(2),cur_ind(3))-...  % Conduction of thermal energy
                                    chiX.*pi.*(0.001.*part_si(tt)) .*sqrt(N_Re).*N_Sc.^(1.0./3.0).*ls.* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc1-rhov(cur_ind(1),cur_ind(2),cur_ind(3))) + ...% away from hailstone term
                                    ice_accret .* cpi .* T(cur_ind(1),cur_ind(2),cur_ind(3)) + ...         %term associated with heating up collected ice
                                    mass_accret .* (lf + cpw.*T(cur_ind(1),cur_ind(2),cur_ind(3))) + ...   %terms associated with freezing liquid and heating up accreted liquid
                                    available_mass./delt .* lf)./...  + soaked.*lf./delt)./...                              %Freezing the remaining unfrozen mass on outside
                                    (chiX.*pi.*(0.001.*part_si(tt)).*sqrt(N_Re).*N_Pr.^(1.0./3.0).*kT(cur_ind(1),cur_ind(2),cur_ind(3)) + ...  %thermal conduction
                                    ice_accret .* cpi + mass_accret .* cpw)); 
                           
                            %diffusion growth coefficient defined here based on Reynolds # regime
                                   diffcoef = chiX .* sqrt(N_Re).*N_Sc.^(1.0./3.0);
                                
                          elseif N_Re >= 2e4
                                chiX = 0.57 + 9.0e-6 .* N_Re;
                                make_it_hail = abs(Tarray - ...
                                    (chiX.*pi.*(0.001.*part_si(tt)).*sqrt(N_Re).*N_Pr.^(1.0./3.0).*kT(cur_ind(1),cur_ind(2),cur_ind(3)).*T(cur_ind(1),cur_ind(2),cur_ind(3))-...  % Conduction of thermal energy
                                    chiX.*pi.*(0.001.*part_si(tt)) .*sqrt(N_Re).*N_Sc.^(1.0./3.0).*ls.* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc1-rhov(cur_ind(1),cur_ind(2),cur_ind(3))) + ...% away from hailstone term
                                    ice_accret .* cpi .* T(cur_ind(1),cur_ind(2),cur_ind(3)) + ...         %term associated with heating up collected ice
                                    mass_accret .* (lf + cpw.*T(cur_ind(1),cur_ind(2),cur_ind(3))) + ...   %terms associated with freezing liquid and heating up accreted liquid
                                    available_mass./delt .* lf)./... .* lf + soaked.*lf./delt)./...                              %Freezing the remaining unfrozen mass on outside
                                    (chiX.*pi.*(0.001.*part_si(tt)).*sqrt(N_Re).*N_Pr.^(1.0./3.0).*kT(cur_ind(1),cur_ind(2),cur_ind(3)) + ...  %thermal conduction
                                    ice_accret .* cpi + mass_accret .* cpw)); 
                          
                              %diffusion growth coefficient defined here based on Reynolds # regime
                                diffcoef = chiX .* sqrt(N_Re).*N_Sc.^(1.0./3.0);
                           
                           end
                           
                           %Find Index of minimum in T_sfc eqn above, assign that to new particle T
                            min_ind = find(make_it_hail == min(make_it_hail));
                            part_Ts(tt) = Tarray(min_ind(1));   %Added (1) in case there is T perfectly between two Ts (3/18/19)
                            part_Ts(tt) = min(part_Ts(tt),273.15);  %make sure it doesnt go over 0 celsius --> if 0, triggers wet growth below
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % DETERMINE THE DENSITY of RIME ACCRETION
                            % 
                            % Uses the Heymsfield and Pflaum (1985) parameterization of rime density, 
                            % but with 500 kg/m3 as min (otherwise might get too low
                            %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            fRIME = (0.5.*D_c(cur_ind(1),cur_ind(2),cur_ind(3)) .* (part_vel(tt).*vIMP) ./  (T0 - part_Ts(tt)));
                                if fRIME < 1.60 && T(cur_ind(1),cur_ind(2),cur_ind(3)) <= 268.15;
                                    rhoh_f = 0.30 .* fRIME.^0.44 .* 1000; %conversion to kg/m3
                                elseif fRIME >= 1.60;
                                    rhoh_f = 0.30 .* fRIME.^0.44 .* 1000; %conversion to kg/m3
                                else
                                    rhoh_f = 1000.0.*exp(-0.03115-1.7030.*fRIME+0.9116.*fRIME.^2-0.1224.*fRIME.^3);
                                end

                                    rhoh_f = min(rhosolid, rhoh_f);  %Max allowed is solid ice
                                    rhoh_f = max(rhomin, rhoh_f);    %Min allowed is min density threshold, default 500
                                    rhoh(tt) = rhoh_f;               %Dump it into hail density for this timestep
                           
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           %  Trigger for Wet Growth. Here we allow ice collection, and appropriate enthalpy %
                           %   becomes lv (rather than ls). Also, Tsfc = 0. Explicitly predicts Frozen       %
                           %   Fraction following Nelson (1983, JAS)                                         %
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                           rhossfc0=es0./Rv./T0;   %Equilibrium vapor pressure is at Tsfc=0
                            if part_Ts(tt) == 273.15;% && mass_accret > 0.0;    
                                which_l = lv;           %enthalpy is for vaporization
                                rhossfc = rhossfc0;   %Force surface rho to be rhossfc0
                                Eci = 1.0;    %switch on collection of ice
                                %Collection of ice crystals only if wet growth
                                Keri = pi ./ 4 .* Eci .* (0.001 .* part_si(tt)).^2 .* part_vel(tt);  %Kernel for ice collection
                                ice_accret =  IWC_total(cur_ind(1),cur_ind(2),cur_ind(3)) .* Keri;  %Collection Kernel assuming continuous collection of IWC, with eff =1  only for wet growth!
                                
                                                             
                                if N_Re < 6000
                                    FrozenFrac = (-2.0.*pi .* (part_si(tt) .* 0.001) .* ventH .* kT(cur_ind(1),cur_ind(2),cur_ind(3)) .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) +...
                                    2.0.* pi .* (part_si(tt) .* 0.001) .* ventV .* lv .* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc0-rhov(cur_ind(1),cur_ind(2),cur_ind(3)))-... - rhossfc0) - ...
                                    ice_accret .* cpi .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) - mass_accret .* cpw .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0))./...
                                    (mass_accret .* lf + available_mass.*lf./delt );%+ soaked.*lf./delt);   %freeze the remaining unfrozen mass-- added to liquid from accretion
                                
                                %diffusion growth coefficient defined here based on Reynolds # regime
                                   diffcoef = 2.0 .* ventV;
                                
                                elseif N_Re >= 6000 && N_Re < 2e4
                                    chiX = 0.76;  %Note: this value is for spheres. Increases with increasing oblateness
                                    FrozenFrac = (-chiX.*pi.*(part_si(tt).*0.001).*sqrt(N_Re).*(N_Pr.^(1.0./3.0)).* kT(cur_ind(1),cur_ind(2),cur_ind(3)) .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) +...
                                    chiX.*pi.*(part_si(tt).*0.001).*sqrt(N_Re).*(N_Sc.^(1.0./3.0)).*lv.* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc0-rhov(cur_ind(1),cur_ind(2),cur_ind(3)))-... - rhossfc0) - ...
                                    ice_accret .* cpi .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) - mass_accret .* cpw .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0))./...
                                    (mass_accret .* lf + available_mass.*lf./delt );% soaked.*lf./delt);   %freeze the remaining unfrozen mass-- added to liquid from accretion
                              
                                %diffusion growth coefficient defined here based on Reynolds # regime
                                   diffcoef = chiX .* sqrt(N_Re).*N_Sc.^(1.0./3.0);
                                
                                elseif N_Re >= 2e4
                                    chiX = 0.57 + 9e-6 .* N_Re;
                                    FrozenFrac = (-chiX.*pi.*(part_si(tt).*0.001).*sqrt(N_Re).*(N_Pr.^(1.0./3.0)).* kT(cur_ind(1),cur_ind(2),cur_ind(3)) .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) +...
                                    chiX.*pi.*(part_si(tt).*0.001).*sqrt(N_Re).*(N_Sc.^(1.0./3.0)).*lv.* Dv(cur_ind(1),cur_ind(2),cur_ind(3)) .* (rhossfc0-rhov(cur_ind(1),cur_ind(2),cur_ind(3)))-... - rhossfc0) - ...
                                    ice_accret .* cpi .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0) - mass_accret .* cpw .* (T(cur_ind(1),cur_ind(2),cur_ind(3))-T0))./...
                                    (mass_accret .* lf + available_mass.*lf./delt);% + soaked.*lf./delt);   %freeze the remaining unfrozen mass-- added to liquid from accretion
                                
                                %diffusion growth coefficient defined here based on Reynolds # regime
                                    diffcoef = chiX .* sqrt(N_Re).*N_Sc.^(1.0./3.0);   
                                
                                end
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % Rasmussen and Heymsfield 1987 Ff Parameterization -- uncomment   %
                                %   -- Does not explicitly predict Ff, but makes bigger hail       %
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               %  if (LWC_cloud(cur_ind(1),cur_ind(2),cur_ind(3))+LWC_rain(cur_ind(1),cur_ind(2),cur_ind(3))) >= 0.002
                               %     FrozenFrac2 = 0.25 + (0.75./(1.0 + 0.1798.*(record_liqr(tt) + record_liqc(tt) - 0.002)))-0.2;  %RH87 FF parameterization
                               %  else
                               %      FrozenFrac2 = 1.0;
                               %  end
                               %
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
                               %Flags for where diffusion mass change goes
                                 diff_to_ice = 0;  %Diffusion growth/loss not to ice if wet growth
                                 diff_to_liq = 1;  %Diffusion growth/loss to liquid if wet growth 
                            
                            else   %SFC TEMP NOT 273.15 -- i.e., not wet growth
                                FrozenFrac = 1.0;       %Frozen fractions should be unity
                                FrozenFrac2 = 1.0;      %
                                which_l = ls;           %enthalpy is sublimation
                                Eci = 0.0;              %turn off ice collection
                                soaked = 0.0;           %clear out unfrozen liquid if not in wet growth Apr 1 *NOTE added soaked into heat balance Apr1
                                available_mass = 0.0;   %clear out unfrozen liquid if not in wet growth Apr 1 
                                esTs = es0 .* exp(which_l./Rv .* (1./T0-1./part_Ts(tt)));   %equilibrium vapor pressure as function of Ts (Pa)
                                rhossfc=esTs./Rv./part_Ts(tt);                            %equilibrium vapor density over hailstone sfc (kg/m3)
                                
                                %Flags for where diffusion mass change goes
                                 diff_to_ice = 1;  %Diffusion growth/loss to ice if dry growth
                                 diff_to_liq = 0;  %Diffusion growth/loss not to liquid if dry growth 
                            end
                                    
                               %Clean up FrozenFrac in case something bad happens     
                                FrozenFrac = min(1.0, FrozenFrac);  %Ensure Ff not > 1.0
                                FrozenFrac = max(0.0, FrozenFrac);  %Ensure Ff not < 0.0
                                FF(tt) = FrozenFrac;  %Dump into FrozenFrac for this timestep
                                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % MASS GROWTH FROM DIFFUSION  dm/dt|diff%
                            %
                            % Uses Diffusion coefficients determined from the Reynolds number regime above 
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
                            mass_diff = -1.0.*pi.*Dv(cur_ind(1),cur_ind(2),cur_ind(3)).*(0.001 .* part_si(tt)) .* ...
                                diffcoef .* ...   %Mass growth rate from diffusion, k
                                (rhossfc-rhov(cur_ind(1),cur_ind(2),cur_ind(3)));  
                            
                            %Growth rates (kg/s) for determining increase in size and mass
                            diffgrow(tt) = mass_diff;
                            massgrow(tt) = (mass_accret + ice_accret + available_mass./delt) .*FrozenFrac;  %Liquid on outside retained from previous calc; this is total mass frozen in timestep
                                     %includes outside and inside% 
                            accretgrow(tt) = (mass_accret + ice_accret).*FrozenFrac;  %acrreted mass frozen (not counting retained liquid       
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % COMBINE GROWTH - accreted liquid and ice, and diffusion, if applicable   %
                            %     - diffusion change in mass may be from liquid -- uses flag 0/1       %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            GROWTH_MICRO(tt) = (accretgrow(tt) + diffgrow(tt).*diff_to_ice) .* delt; %growth increment in the time step (kg), called delta m in notes
                                
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          % SOAKING MODULE - allows unfrozen water to soak into hailstone if density is %
                          %less than solid ice. Using average density for this calculation.             % 
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if FrozenFrac < 1.0   %Only if part of accreted mass is unfrozen (otherwise freezes on contact)
                                %Compute how much mass is needed to fill in to get solid density of ice (air pockets)
                                needed_to_soak = (rhosolid .* pi./6 .* (part_si(tt).*0.001).^3.0) - tot_m; %mass needed to soak and densify (kg)
                                needed_to_soak = max(needed_to_soak, 0.0);  %make sure it doesnt drop below 0 kg!
                                
                                %Compute the mass of unfrozen liquid available to soak (kg) 
                                available_mass = (1-FrozenFrac)*(mass_accret + ...    %liquid mass accreted available to soak (kg); 
                                    ice_accret + available_mass./delt)*delt + diff_to_liq.*diffgrow(tt).*delt; %includes unfrozen from previous time step
                                                                                    %as well as any liquid added/lost to diffusion
                              
                                if available_mass >= needed_to_soak   %IF more is available than needed to soak
                                 
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %Shedding condition, following RH87   %
                                    %$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   soaked = needed_to_soak;                                         %soaked mass of liquid yet to be frozen APR1
                                   tot_m = tot_m + needed_to_soak + GROWTH_MICRO(tt);               %This accounts for ice growth, plus soaked. 
                                   rhoh(tt) = ((1.0 - 0.08.*FrozenFrac).*FrozenFrac) .* 1000;       %RH87 spongy growth parameterization
                                   rhoh(tt) = min(rhosolid, rhoh(tt));                              %making sure we dont get super-dense ice
                                   rhoh(tt) = max(rhomin, rhoh(tt));                                %making sure we dont get super-fluff ice. CHANGED TO 0 APR 1 $$$$
                                   remaining_unfrozen_mass = available_mass - needed_to_soak;       %Available for shedding on outside (kg)
                                   mwcrit = (2.68e-4 + 0.1389 .* tot_m);                            %Critical mass of liquid retainable following RH87; includes soaked and ice mass added in this timestep
                                   remaining_unfrozen_mass = min(mwcrit, remaining_unfrozen_mass);  %Shed anything in excess of mwcrit; this should be passed to next timestep
                                   available_mass = remaining_unfrozen_mass;                        %no longer available, its passed to remaining unfrozen 
                                    
                                else
                                    
                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   %Soaked but no shedding -- still wet  %
                                   %$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   soaked = available_mass;                                        %soaked but not frozen
                                   tot_m = tot_m + available_mass + GROWTH_MICRO(tt);              %this accounts for mass frozen inside plus diffusion 
                                   rhoh(tt) = ((1.0 - 0.08.*FrozenFrac).*FrozenFrac) .* 1000;      %RH87 spongy growth parameterization
                                   rhoh(tt) = min(rhosolid, rhoh(tt));                             %making sure we dont get super-dense ice.
                                   rhoh(tt) = max(rhomin, rhoh(tt));                               %making sure we dont get super-fluff ice. 
                                   available_mass = 0.0;                                           %All of the available mass was soaked 
                                
                                end
                            else   %If not in wet growth, no soaking, shedding, and Ff=1.
                                available_mass = 0;                 %none available for soaking or shedding if dry growth
                                soaked = 0;
                                tot_m = tot_m + GROWTH_MICRO(tt);   %if no wet growth/soaking, then total mass increase is just what grew in timestep.
                                remaining_unfrozen_mass = 0;        %ensure that everything was frozen, nothing carried into next time step.
                            end
                            %%%%%%%%% End soaking module%%%%%%%%%% 
               
                    
                      end   %%%END OF IF STATEMENT ABOUT GROWTH ABOVE/BELOW 0 Celsius
             
                     %cleaning up (make sure no NaNs)
                        cleanup = find(isnan(GROWTH_MICRO) == 1);
                        GROWTH_MICRO(cleanup) = 0;
                    
                   %UPDATE PARTICLE SIZE, VEL, XYZ POSITIONS, and TEMPERATURE
                   % Note: accretion/frozen growth assumes rhoh at that
                   % timestep, which may be < solid. Deposition growth
                   % assumed to be solid, but could be changed
                    if GROWTH_MICRO(tt) == 0            %No growth = no size change!
                        part_si(tt+1) = part_si(tt);
                    else
                        %If there is growth (of ice only); uses diffusion flag to put into ice if applicable
                        part_si(tt+1) = (2 .* ((0.001.*part_si(tt)).^3 ./ 8 + (3.* massgrow(tt)./4./pi./rhoh(tt) + 3.* diff_to_ice.*diffgrow(tt)./4./pi./rhosolid)).^(1.0/3.0)) .* 1000.0;
                    end
                    
                    avg_rhoh(tt+1) = tot_m ./ (pi./6 .* (part_si(tt+1).*0.001).^3)  ;%average particle density
                    ah = sqrt(1.333333 .* avg_rhoh(tt+1)./ Cdrag .* g);   %Update ah parameter based on density
                    part_vel(tt+1) = (ah .* (part_si(tt+1)./1000).^bh).*sqrt(sfcdens./comp_rhoa(cur_ind(1),cur_ind(2),cur_ind(3)));  %update fall speed
                    part_posx(tt+1) = part_posx(tt) + comp_u(cur_ind(1),cur_ind(2),cur_ind(3)).*delt.*0.001;                %update x position in km
                    part_posy(tt+1) = part_posy(tt) + comp_v(cur_ind(1),cur_ind(2),cur_ind(3)).*delt.*0.001;                %udpate y position in km
                    part_posz(tt+1) = part_posz(tt) + (comp_w(cur_ind(1),cur_ind(2),cur_ind(3))-part_vel(tt+1)).*delt.*0.001; %update z position in km

                %Find nearest index to where hailstone has ended up
                temp = abs(part_posx(tt+1)-xgrid);      %
                mintemp = find(temp == min(temp));
                cur_ind(1) = mintemp(1);   %update 8/29
                %cur_ind(1) = find(temp == min(temp));   %finding closest grid point
                temp = abs(part_posy(tt+1)-ygrid);      %so that we can use those
                mintemp = find(temp == min(temp));
                cur_ind(2) = mintemp(1);   %update 8/29
                %cur_ind(2) = find(temp == min(temp));   %values of u,v,w, and process
                tempz = abs(part_posz(tt+1)-zgrid);     %rates for next set of calculations
                mintempz = find(tempz == min(tempz));
                cur_ind(3) = mintempz(1);
                %cur_ind(3) = find(tempz == min(tempz)); % 

                %Record useful storm/hailstone quantities for tracking 
                record_Theta(tt+1)=comp_theta(cur_ind(1),cur_ind(2),cur_ind(3));    %theta at hailstone level
                record_liqc(tt+1) =LWC_cloud(cur_ind(1),cur_ind(2),cur_ind(3));     %Cloud liquid water content kg/m3
                record_liqr(tt+1) =LWC_rain(cur_ind(1),cur_ind(2),cur_ind(3));      %Rain iquid water content kg/m3
                record_ice(tt+1) = IWC_total(cur_ind(1),cur_ind(2),cur_ind(3));     %Ice water content kg/m3 (includes snow and cloud ice)
                record_qv(tt+1) = comp_qv(cur_ind(1),cur_ind(2),cur_ind(3));        %ambient vapor density (kg/m3)
                record_rhoa(tt+1) = comp_rhoa(cur_ind(1),cur_ind(2),cur_ind(3));    %ambient air density (kg/m3)
                record_u(tt+1) = comp_u(cur_ind(1),cur_ind(2),cur_ind(3));          %u component of wind (m/s)
                record_v(tt+1) = comp_v(cur_ind(1),cur_ind(2),cur_ind(3));          %v component of wind (m/s)
                record_w(tt+1) = comp_w(cur_ind(1),cur_ind(2),cur_ind(3));          %w component of wind (m/s)
                record_T(tt+1) = (T(cur_ind(1),cur_ind(2),cur_ind(3)));             %ambient air Temperature (K)
                record_sfcliq(tt+1) = available_mass;                               %mass of liquid on hailstone surface (kg)
                record_soaked(tt+1) = soaked;                                       %mass of liquid soaked into hailstone per timestep (kg)              
                BUGGYSMALLZ(tt) =soaked;                                            %Useful for debugging things
                
                tt=tt+1;  %Advance time by one timestep
               
                end  %-----------END GROWTH LOOP

             %Check to make sure particle is nonexistent if it hits ground
                if part_posz(tt+1) <= 0
                     part_posz(tt+1:end) = NaN;
                     part_posx(tt+1:end) = NaN;
                     part_posy(tt+1:end) = NaN;
                     GROWTH_MICRO(tt+1:end) = NaN;
                     part_si(tt+1:end) = NaN;
                     part_Ts(tt+1:end) = NaN;
                     record_T(tt+1:end) = NaN;
                     record_Theta(tt+1:end) = NaN;
                     record_liqr(tt+1:end) = NaN;
                     record_liqc(tt+1:end) = NaN;
                     record_ice(tt+1:end) = NaN;
                     record_qv(tt+1:end) = NaN;
                     record_rhoa(tt+1:end) = NaN;
                     record_u(tt+1:end) = NaN;
                     record_v(tt+1:end) = NaN;
                     record_w(tt+1:end) = NaN;
                     avg_rhoh(tt+1:end) = NaN;
                     diffgrow(tt+1:end) = NaN;
                     massgrow(tt+1:end) = NaN;
                     accretgrow(tt+1:end) = NaN;
                     FF(tt+1:end) = NaN;
                     BUGGYSMALLZ(tt+1:end) = NaN;
                end
                
            %DETERMINE THE MAXIMUM SIZE ATTAINED BY EACH HAILSTONE
                part_max_si(pp)=max(part_si);
           
            %Construct the grid showing the final size of a hailstone initiated from an embryo at that point.
                tstgr(pp,yi,zz) = part_max_si(pp);

            %Construct the final location of an embryo starting at given index at that point.
                final_posx(pp,yi,zz) = part_posx(tt);   %final x position
                final_posy(pp,yi,zz) = part_posy(tt);   %final y position
                final_posz(pp,yi,zz) = part_posz(tt);   %final z position
                total_traj_time(pp,yi,zz) = tt;         %total trajectory time
   
            end  %LOOP OVER PARTICLES ----- END
         
         end  %%%LOOP OVER Y POSITIONS -------- END
 end  %%% LOOP OVER Z POSITIONS -----------------------------END

