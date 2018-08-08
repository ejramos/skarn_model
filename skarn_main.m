%file: skarn_main.m

%author:  Evan J. Ramos

%Description: Model for the Empire Mountain skarn found in Ramos et al.,
%2018 G3.

close all
clear
clc

%% Simulation parameters

filename                   = 'mineralisotopeconstant.txt';

%% Dimensional constants

%Space
%%%%
Depth    = 5000;        %[m]
Length   = 7880;        %[m]
Width    = 4000;        %[m] z-direction width of Empire Mountain pluton

%Time
%%%%
yrs2s    = 3600*24*365; %conversion factor to seconds from years
dt       = 1e8;         %[s] time step
endtime  = 20e3*yrs2s;
t_totals = 200;         %Total number of time steps saved
times    = (1:t_totals)*endtime/t_totals;


%Darcy
%%%%
Mu       = 8.9e-4;      %[Pa*s] Dynamic viscosity of water at 20 deg C %%%%%%%%%%
rho_f    = 1000;        %[kg/m^3] Density of water
rho_s    = 2800;        %[kg/m^3] Density of solid
rho_m    = 1500;        %[kg/m^3] Density of magmatic fluid
grav     = 9.81;        %[m/s^2] Acceleration due to gravity on Earth
tensile  = 20e6;        %[Pa] average tensile strength of limestone (Mark Zoback)
phi_min  = 1e-4;        %[m^3/m^3] minimum imposed porosity
phi_0    = 5e-3;        %[m^3/m^3] reference porosity
k_0      = 1e-17;       %[m^2] reference permeability
n        = 3;           %integer exponent for power law calc
kvkh     = 1;


%Temperature
%%%%
T_s      = 293;         %[K]
qt_s     = .06;         %[W/m^2] surface heat flux 
A        = 2e-6;        %[W/m^3] radiogenic heat production: average
kf       = 0.65;        %[W/m/K] thermal conductivity of water
ks       = 2.08;        %[W/m/K] thermal conductivity of rock
cp_f     = 1046;        %[J/(kg K)] specific heat of water
cp_s     = 908;         %[J/(kg K)] specific heat of solid

%pluton
phi_p    = .0025;             %[m^3/m^3] porosity of pluton
rho_p    = 2800;              %[kg/m^3] pluton density
cp_p     = 1450;              %[J/kg/K] specific heat 
Tp       = 1355;              %[K] pluton temperature
Tliq     = 1355;              %[K] liquidus temperature
Tsol     = 983;               %[K] solidus temperature
t_inj    = dt;                %[s] timescale over which pluton injects its fluids
wt_perc  = .0089;             %[1] weight percent volatiles in magma %%%%%%%%%%%%%%%%%%%%%
Q_gamma  = wt_perc*rho_m/(rho_f*t_inj); %[kg/m^3/s] fluid source injection rate 
Q_latent = 190e3;             %[J/kg] total latent heat of crystallization
dT       = Tliq - Tsol;       %[K] temperature range of crystallization (Hanson, 1995)

dHdT_l   = Q_latent/dT;       %[J/kg/K] latent heat of crystallization pluton
melt     = @(T) (T - Tsol)/dT;%[1] function that computes melt fraction
T_crit   = Tsol;%773;               %[K] temperature at which pluton becomes permeable %%%%REVIEWER suggestion
ho_f     = 117390;            %[J/kg] specific enthalpy of water at STP

%Stable isotopes
%%%%
tau      = sqrt(2);        %[m/m] tortuosity
vphi_m   = 1;              %[m^3/m^3] volume fraction of mineral in solid part
Do       = 1e-9;           %[m^2/s] solute diffusivity in water (Marc 2016)
delta_s  = -10;            %[per mil] surface delta 18O value
delta_p  = 6.4;              %[per mil] plutonic delta 18O value (D'Errico 2010)
R_std    = .0020052;       %[1] 18O/16O of VSMOW
R_gamma  = delta_p/1000*R_std + R_std;
xw       = .89;            %[kg/kg] mass ratio of oxygen to total molar mass

%Following Clechenko and Valley (2003)
D_andr = -3.03 - .57*.6 + 2.51; %Grs-qtz + Grs-andr correction + qtz-h2o
E_andr = 0;
F_andr = 0 + -1.46; %Grs-qtz + qtz-h2o %1.96

alpha_andr = @(T) exp(D_andr*1000./(T.^2) + E_andr./T + F_andr/1000);
delta_andr = @(T,delta_w) alpha_andr(T).*(1000 + delta_w) - 1000;
%convert deltas to ratios
Rwi      = (delta_s/1000)*R_std + R_std;

%% Build Grid and operators

Nx = 400; Ny = 100;
% Defining the computational grid and conductivities
Grid.xmin = 0; Grid.xmax = Length;  Grid.Nx = Nx;
Grid.ymin = 0; Grid.ymax = Depth;   Grid.Ny = Ny; 
Grid.Nz = 1;
Grid.dz = Width;
Grid.psi_dir = 'xy';
Grid      = build_grid(Grid);
[Xc,Yc]   = meshgrid(Grid.xc,Grid.yc);
[Xf,Yf]   = meshgrid(Grid.xf,Grid.yf);
[Xfx,Yfx] = meshgrid(Grid.xf,Grid.yc);
[Xfy,Yfy] = meshgrid(Grid.xc,Grid.yf);

%Operators
[D,G,I]           = build_ops(Grid);

%% Create Medium
Emp_coords = make_emp_structure;
[Mineral, phi, Kd, k, Pluton] = make_EMP_medium(filename,Grid,Emp_coords);
phi_m                         = vphi_m*(1-phi);
phi_lim_d = phi(Grid.dof_xmin);

%% Characteristic scaling terms

%Topography driven flow
H      = Depth;         %[m]
dh     = .1*H;          %[m] reviewer suggestion

%Heat transport
kbar   = phi/tau*kf + (1-phi)*ks;
cp_s_i = cp_s + dHdT_l;
rhobar = phi*rho_f*cp_f + (1-phi)*rho_s*cp_s;
rhobar(Pluton.dof) = phi(Pluton.dof)*rho_f*cp_f + ...
                     (1-phi(Pluton.dof))*rho_p*cp_s_i;
qt_b   = qt_s - A*H;                             %[W/m^2] basal heat flux

Teq    = @(z) (A/kbar(Grid.dof_ymin(1)))*(H*z - .5*(H^2 + z.^2)) + ...
              qt_s/kbar(Grid.dof_ymin(1))*(H-z) + T_s; %[K]

%% Boundary conditions

% Parameters: pressure
Param.p.dof_dir   = Grid.dof_ymax;
Param.p.dof_f_dir = Grid.dof_f_ymax;
Param.p.g         = rho_f*grav*(linspace(0,dh,Grid.Nx)');
Param.dof_out     = Param.p.dof_dir;
Param.p.dof_neu   = [];
Param.p.dof_f_neu = [];
Param.p.qb        = [];

%% Build boundaries

%Pressure
[Bp,Np,fn_p] = build_bnd(Param.p,Grid);

%% Initialize pressure, temperature, and fluid stable isotope values

%%%% steady state initial conditions %%%%%
load EMP_ss_NEW
%pressure
Ps               = zeros(Grid.N,t_totals);
p_litho          = reshape(rho_s*grav*(H-Yc),Grid.N,1);
Ps(:,1)          = p_litho/rho_s*rho_f;
fs_ps            = zeros(t_totals,1);
%fluxes
Qx               = zeros(Grid.N,t_totals);
Qy               = zeros(Grid.N,t_totals);
Qx_null          = zeros(Grid.Nfx,t_totals);
Qy_null          = zeros(Grid.Nfy,t_totals);
Psi              = zeros((Grid.Nx+1)*(Grid.Ny+1),t_totals);

%Heat transport operators
theta_a          = 1;                          %advection solved explicitly
theta_d          = 0;                          %diffusion solved implicitly
Ts               = zeros(Grid.N,t_totals+1);
Ts(:,1)          = T_i;
Ts(Pluton.dof,1) = Tp;
fs_T             = A*ones(Grid.N,1);
fs_T(Pluton.dof) = 0;


%Stable isotope transport operators, initialized
theta              = 0;
[~,~,~,del_w,del_m]= make_avg_rxn_terms(Grid,Mineral,Ts(:,1));
Rs                 = zeros(2*Grid.N,t_totals+1);
Rs(1:Grid.N,1)     = del_w/1000*R_std + R_std;       %Rw
Rs(Grid.dof(~ismember(Grid.dof,Pluton.dof)),1) = delta_s/1000*R_std + R_std;
Rs(Grid.N+1:end,1) = del_m/1000*R_std + R_std;       %Rm
delta_ws           = zeros(Grid.N,t_totals+1);
delta_ms           = zeros(Grid.N,t_totals+1);
delta_gt           = zeros(Grid.N,t_totals+1);
delta_ws(:,1)      = del_w;
delta_ms(:,1)      = del_m;
fs_O               = zeros(2*Grid.N,1);

%diagonalized matrices
Phi                = comp_mean(reshape(phi,Grid.Ny,Grid.Nx),-1,kvkh,Grid);
Phi_diag           = spdiags(phi,0,Grid.N,Grid.N);
Phi_m              = spdiags(phi_m,0,Grid.N,Grid.N);
Rhobar             = spdiags(rhobar,0,Grid.N,Grid.N);
Kbar               = comp_mean(reshape(kbar,Grid.Ny,Grid.Nx),1,kvkh,Grid);
Xw                 = spdiags(xw*ones(Grid.N,1),0,Grid.N,Grid.N);
Rho_f              = rho_f;
KdMu               = Kd/Mu;

%GW Age
Ages_surf          = zeros(Grid.N,t_totals+1); %surface fluid tracer
Ages_volc          = zeros(Grid.N,t_totals+1); %metavolcanic fluid tracer
Ages_plut          = zeros(Grid.N,t_totals+1); %plutonic fluid tracer
Ages_lime          = zeros(Grid.N,t_totals+1); %limestone fluid tracer

volc_dof = Grid.dof(Mineral.number == 5 | Mineral.number == 6);
lime_dof = Grid.dof(Mineral.number == 4);
Ages_volc(volc_dof,1) = 1;
Ages_lime(lime_dof,1) = 1;
Ages_plut(Pluton.dof,1) = 1;

%Initialize loop variables
t = 0;
count = 1;
%temperature
T_i  = Ts(:,1);
T_ii = Ts(:,1);
%tracers
Ages_surf_i = Ages_surf(:,1);
Ages_volc_i = Ages_volc(:,1);
Ages_plut_i = Ages_plut(:,1);
Ages_lime_i = Ages_lime(:,1);

%Oxygen isotopes
R_i = Rs(:,1);
while t < endtime
    
    t = t + dt;
    %% Solve steady state pressure --> Darcy Flux
    %pressure
    Lp             = -D*KdMu*G;
    fs_grav        = D*(KdMu*Rho_f*grav)*G*Yc(:);
    fs_p           = zeros(Grid.N,1);

    if t > dt
        %tracking how much crystallization has occurred between the
        %subsequent time steps so that the appropriate scale to the source
        %term can be applied to the magmatic fluid input.
        diff_melt_i        = melt(T_i(Pluton.dof));
        diff_melt_ii       = melt(T_ii(Pluton.dof));
        diff_melt_i(diff_melt_i < 0) = 0;
        diff_melt_ii(diff_melt_ii < 0) = 0;
        diff_melt = diff_melt_ii - diff_melt_i;
        diff_melt(diff_melt<0) = 0;
        fs_p(Pluton.dof) = fs_p(Pluton.dof) + Q_gamma*diff_melt;
        fs_T(Pluton.dof) = fs_p(Pluton.dof).*new_rho(Pluton.dof)*ho_f;
        fs_O(Pluton.dof) = fs_p(Pluton.dof)*R_gamma*xw;                
    end

    fp             = fs_p + fn_p + fs_grav;
    P              = solve_lbvp(Lp,fp,Bp,Param.p.g,Np);
    %STORE AT APPROPRIATE TIME STEP
    if xor(t == times(count),t > times(count))
        Ps(:,count+1) = P;
        fs_ps(count)  = sum(fs_p);
    end
    q              = comp_flux_p(D,KdMu,Rho_f*grav*G*Yc(:),G,P,fs_p,Grid,Param.p);
    
    qtop           = q(Grid.dof_f_ymax);
    qdown_dof      = Grid.dof_ymax(qtop < 0);
    qdown_dof_f    = Grid.dof_f_ymax(qtop < 0);
    
    [qx,qy,qmags]  = center_flux(q,Grid);
    [psi,~,~]      = comp_streamfun(q,Grid);
    %STORE AT APPROPRIATE TIME STEP
    if xor(t == times(count), t > times(count))
        Qx(:,count)    = reshape(qx,Grid.N,1);
        Qy(:,count)    = reshape(qy,Grid.N,1);
        Qx_null(:,count) = q(1:Grid.Nfx);
        Qy_null(:,count) = q(Grid.Nfx+1:end);
    end
    Psi(:,count)       = reshape(psi,numel(psi),1);
    % 1st order upwind of fluid fluxes
    Aq = flux_upwind(q,Grid);
    
    %% Solve GW tracer
    molten_dof        = find(T_i > Tsol);
    molten_dof        = intersect(molten_dof,Pluton.dof);
    Da                = D(molten_dof,:);
    molten_dof_f      = (Grid.dof_f(abs(sum(Da,1)) > eps))';

    %Surface tracer boundary conditions
    Param.as.dof_dir   = qdown_dof;
    Param.as.dof_f_dir = qdown_dof_f;
    Param.as.g         = ones(length(qdown_dof),1);
    Param.as.dof_neu   = [];
    Param.as.dof_f_neu = [];
    Param.as.qb        = [];
    [BA.s, NA.s, fn_A.s] = build_bnd(Param.as,Grid);

    %Metavolcanic tracer boundary conditions
    Param.av.dof_dir   = [];
    Param.av.dof_f_dir = [];
    Param.av.g         = [];
    Param.av.dof_neu   = [];
    Param.av.dof_f_neu = [];
    Param.av.qb        = [];
    [BA.v, NA.v, fn_A.v] = build_bnd(Param.av,Grid);

    %Plutonic tracer boundary conditions
    Param.ap.dof_dir   = molten_dof;
    Param.ap.dof_f_dir = molten_dof_f;
    Param.ap.g         = ones(length(molten_dof),1);
    Param.ap.dof_neu   = [];
    Param.ap.dof_f_neu = [];
    Param.ap.qb        = [];
    [BA.p, NA.p, fn_A.p] = build_bnd(Param.ap,Grid);

    %Limestone tracer boundary conditions
    Param.al.dof_dir   = [];
    Param.al.dof_f_dir = [];
    Param.al.g         = [];
    Param.al.dof_neu   = [];
    Param.al.dof_f_neu = [];
    Param.al.qb        = [];
    [BA.l, NA.l, fn_A.l] = build_bnd(Param.al,Grid);

    alpha             = 0; %explicit-implicit integer computation
    fs_A              = zeros(Grid.N,1);
    ImA               = Phi_diag + dt*D*Aq*(1-alpha);
    ExA               = Phi_diag - dt*D*Aq*alpha;
    RHS.s             = (fs_A + fn_A.s)*dt + ExA*Ages_surf_i;
    RHS.v             = (fs_A + fn_A.v)*dt + ExA*Ages_volc_i;
    RHS.p             = (fs_A + fn_A.p)*dt + ExA*Ages_plut_i;
    RHS.l             = (fs_A + fn_A.l)*dt + ExA*Ages_lime_i;

    Ages_surf_new = solve_lbvp(ImA,RHS.s,BA.s,Param.as.g,NA.s);
    Ages_volc_new = solve_lbvp(ImA,RHS.v,BA.v,Param.av.g,NA.v);
    Ages_plut_new = solve_lbvp(ImA,RHS.p,BA.p,Param.ap.g,NA.p);
    Ages_lime_new = solve_lbvp(ImA,RHS.l,BA.l,Param.al.g,NA.l);
    %store for next time step
    Ages_surf_i   = Ages_surf_new;
    Ages_volc_i   = Ages_volc_new;
    Ages_plut_i   = Ages_plut_new;
    Ages_lime_i   = Ages_lime_new;
    %STORE AT APPROPRIATE TIME STEP
    if xor(t == times(count), t > times(count))
        Ages_surf(:,count+1) = Ages_surf_new;
        Ages_volc(:,count+1) = Ages_volc_new;
        Ages_plut(:,count+1) = Ages_plut_new;
        Ages_lime(:,count+1) = Ages_lime_new;
    end
        
    
    %% Build parameters
    % Parameters: temperature
    Param.T.dof_dir   = qdown_dof;
    Param.T.dof_f_dir = qdown_dof_f;
    Param.T.g         = T_s*ones(length(qdown_dof),1);
    Param.T.dof_neu   = Grid.dof_ymin(~ismember(Grid.dof_ymin,Pluton.dof));
    Param.T.dof_f_neu = Grid.dof_f_ymin(~ismember(Grid.dof_f_ymin,Pluton.dof_f));
    Param.T.qb        = qt_b;
    
    %Parameters: stable isotope
    Param.O_w.dof_dir   = [qdown_dof; molten_dof];
    Param.O_w.dof_f_dir = [qdown_dof_f; molten_dof_f];
    Param.O_w.g         = [Rwi*ones(length(qdown_dof),1); R_gamma*ones(length(molten_dof),1)];
    Param.O_w.dof_neu   = [];
    Param.O_w.dof_f_neu = [];
    Param.O_w.qb        = [];
    
    %Stable isotope
    [BO,NO,fn_O] = build_bnd(Param.O_w,Grid,'coupled');
    %Temperature
    [BT,NT,fn_T] = build_bnd(Param.T,Grid);
    
    %% Solve for temperature
    % Build heat operators
    Im_T             = @(theta_a,theta_d) ...
                     Rhobar + dt*(1-theta_a)*D*(rho_f*cp_f*Aq) ...
                            - dt*(1-theta_d)*D*(Kbar*G); %%%%%%
    Ex_T             = @(theta_a,theta_d) ...
                     Rhobar - dt*theta_a*D*(rho_f*cp_f*Aq)...
                            + dt*theta_d*D*(Kbar*G);     %%%%%%
    T_new      = solve_lbvp(Im_T(theta_a,theta_d),...
                            Ex_T(theta_a,theta_d)*T_i + ...
                            dt*(fs_T + fn_T),BT,Param.T.g,NT);
    %store for next time step
    T_ii = T_i;
    T_i  = T_new;
    %STORE AT APPROPRIATE TIME STEP
    if xor(t == times(count), t > times(count))
        Ts(:,count+1) = T_new;
    end
    
    %% Solve for stable isotopes
    %update diagonalized matrices
    [Alpha,Rk,Xm]  = make_avg_rxn_terms(Grid,Mineral,T_new);
     
    %create stable isotope operators
    %dRw/dt PDE
    Im_ww          = @(theta) Phi_diag*Xw...
                     + (1-theta)*dt*...
                     (D*(Aq*xw - Phi/tau*xw*Do*G) + ...
                     (Phi_m*Xm*Rk*Alpha));
    Im_wm          = @(theta) -dt*(1-theta)*(Phi_m*Xm*Rk);
    Ex_ww          = @(theta) Phi_diag.*Xw...
                     - theta*dt*...
                     (D*(Aq*xw - Phi/tau*xw*Do*G) - ...
                     (Phi_m*Xm*Rk*Alpha));
    Ex_wm          = @(theta) theta*dt*(Phi_m*Xm*Rk);
    
    %dRm/dt ODE
    Im_mw          = @(theta) -dt*(1-theta)*Rk*Alpha;
    Im_mm          = @(theta) I + dt*(1-theta)*Rk;
    Ex_mw          = @(theta) dt*theta*Rk*Alpha;
    Ex_mm          = @(theta) I - dt*theta*Rk;
    
    %coupled
    Im             = [Im_ww(theta) Im_wm(theta); Im_mw(theta) Im_mm(theta)];
    Ex             = [Ex_ww(theta) Ex_wm(theta); Ex_mw(theta) Ex_mm(theta)];
    
    %solve stable isotope values
    R_new         = solve_lbvp(Im,Ex*R_i + dt*(fs_O),...
                                BO,Param.O_w.g,NO);
    %store for next time step;                        
    R_i           = R_new;
    if xor(t == times(count),t > times(count))
        Rs(:,count+1)        = R_new;
        delta_ws(:,count+1)  = 1000*(R_new(1:Grid.N)-R_std)/R_std;
        delta_ms(:,count+1)  = 1000*(R_new(Grid.N+1:end)-R_std)/R_std;
        delta_gt(:,count+1)  = delta_andr(T_new,delta_ws(:,count+1));
    end
    
    
    %% Update T-dependent variables
    %update fluid density
    [new_rho,~]   = steam_table(T_new,P,Grid);
    Rho_f         = comp_mean(reshape(new_rho,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    
    %update source terms, diagonalized matrices
    solid            = find(T_new < T_crit);
    solid            = intersect(solid,Pluton.dof);
    fs_frac          = pluton_source(T_new(Pluton.dof));
    solidus          = find(T_new < Tsol);
    rhobar(solidus)  = phi(solidus)*rho_f*cp_f + (1-phi(solidus))*rho_s*cp_s;
    Rhobar           = spdiags(rhobar,0,Grid.N,Grid.N);

    [~,solid_loc]    = ismember(Yc(solid),Grid.yc);
    phi(solid)       = phi_lim_d(solid_loc);
    k(solid)         = perm_from_poro(1e-17,phi(solid),1e-4,5e-3,3);
    Kd               = comp_mean(k,-1,kvkh,Grid);
    KdMu             = Kd/Mu;
    Phi              = comp_mean(reshape(phi,Grid.Ny,Grid.Nx),-1,1,Grid);
    Phi_diag         = spdiags(phi,0,Grid.N,Grid.N);

    if xor(t == times(count),t > times(count))
        fprintf('\n%d of %d iterations\n',count,t_totals)
        count = count + 1;
    end
end

times = [0 times];