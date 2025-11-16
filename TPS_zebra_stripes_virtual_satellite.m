function TPS_zebra_stripes_virtual_satellite(Time_obs,L_obs,MLT_obs,alpha_obs,Flux_obs,Tsim_initial,Tsim_end,E_phi,E_R)
% Function Name: TPS_zebra_stripes_virtual_satellite.m
% Description:
%   This function is an extension to "TPS_zebra_stripes_90deg_e.m" for 
%   simulating electron zebra stripes "observed" by a virtual satellite
%   in the presence of an ad-hoc electric field. This function can be used
%   to reproduce simulation results shown in Figure 8 of Mei et al. (2025)
%
%   Main features:
%     - Plot zebra stripes in electron flux spectrogram based on satellite
%       observation (CIRBE CubeSat)
%     - Simulate zebra stripes that would be observed by a virtual
%     satellite following the same track as the real spacecraft for
%     apples-to-apples comparisons with observation
%
%   Input variables:
%     - Information of a segment of satellite measurement: Time_obs,L_obs,
%   MLT_obs,alpha_obs,Flux_obs
%     - Simulation starting and ending time in the MATLAB datenum format:
%   Tsim_initial,Tsim_end
%     - Ad-hoc electric field set up: E_phi,E_R
%
% Author: Yang Mei
% Affiliation: Laboratory for Atmospheric and Space Physics, University of
% Colorado, Boulder
% Contact: Yang.Mei@colorado.edu
% GitHub: https://github.com/yangmei-spacephysics/Test-Particle-Simulation_90deg-electron_zebra-stripes
%
% Version: 1.0
% Date: 2025-11-16
%
% License: This code is released under the MIT License (see LICENSE file).
%
% Notes:
%   Assumptions: Dipole magnetic field
%                An ad-hoc electrostatic field pointing dawn to dusk
%                For equatorially trapped electrons
%                The first adiabatic invariant is conserved
%   Reference: Mei, Y., Li, X., O Brien, D., Xiang, Z., Zhao, H., Sarris, T., 
% et al. (2025). Characteristics of “zebra stripes” of relativistic electrons 
% unveiled by CIRBE/REPTile-2 measurements and test particle simulations. 
% Journal of Geophysical Research: Space Physics, 130, e2024JA033187. 
% https://doi.org/10.1029/2024JA033187
%
%% Simulation set up
% Simulation time duration and time step
time_step=10; % unit: seconds
tsim_end=10*3600; % sec
Time_sim=0:time_step:tsim_end;
nT=length(Time_sim)-1;
len_T=nT+1;

% Specify the time range in which simulated particles are recorded
rec_time_1=(Time_obs(1)-Tsim_initial)*24*3600;
rec_time_2=(Time_obs(end)-Tsim_initial)*24*3600;
time_sim=24*3600*(Time_obs-Tsim_initial); % time of measurements shifted to the simulation time (5hr+1000sec);

%% particle set up
dR=0.1; % Unit: Earth radii (R_E)
R_aim=[1.14:0.02:1.9,1.95:0.05:2.2 2.3:dR:4.5]'; % The radial distance of the electrons at the beginning of simulation
e_num0=1000; % numbers of the electrons at a L for one energy

phi_step=2*pi/e_num0;
phi_ini0=[0:phi_step:2*pi-phi_step]';
phi_ini=phi_ini0;
dEk=0.0125;  % Unit: MeV
Ek_TPS=[0.1:dEk:1 1.025:0.0125:1.5 1.55:0.025:2 2.1 2.2]; % Initial energy
n_Ek0=length(Ek_TPS);
e_num_Ek=length(R_aim)*(e_num0); % how many particles

Posit_phi_ini=repmat(phi_ini,1,length(R_aim));
Posit_R_ini=repmat(R_aim',length(phi_ini),1);
Ek_ini=ones(e_num_Ek,1)*Ek_TPS;

%% Virtual measurements set up
dPhi_m=pi/5; % delta azimuthal angle where measurements are made


%% Energy channel set up based on REPTile-2 channel energy range (based on Khoo et al.(2022))
RNGe.Ebeg=[0.24 0.24 0.25 0.25 0.26 0.27 0.28 .29 0.30 0.32 0.33 0.35 ...
    0.37 0.39 0.41 0.43 0.46 0.49 0.52 0.56 0.60 0.64 0.69 0.74 0.79 ...
    0.85 0.92 0.99 1.07 1.15 1.24 1.34 1.45 1.57 1.69 1.83 1.98 2.14 ...
    2.32 2.51 2.72 2.94 3.19 3.46 3.74 4.06 4.41 4.77 5.21 5.62];
RNGe.Eend=[0.30 0.35 0.32 0.31 0.32 0.33 0.33 0.34 0.35 0.36 0.37 0.39 ...
    0.41 0.43 0.45 0.47 0.50 0.53 0.57 0.61 0.65 0.69 0.74 0.80 0.86 ...
    0.92 0.99 1.07 1.15 1.25 1.34 1.45 1.57 1.69 1.83 1.98 2.15 2.32 ...
    2.52 2.72 2.95 3.20 3.46 3.75 4.07 4.46 4.84 5.20 5.61 6.08];
eleCh=(RNGe.Ebeg+RNGe.Eend)/2;
%% Background field set up
phi_step_B=phi_step/5;
phi_B(:,1)=0:phi_step_B:2*pi;
R_step=0.1;
R_B=1.1:R_step:8;

B0=Geomagnetic_model_dipole_eq(R_B);
B_total=repmat(B0,length(phi_B),1);

%% Set up an ad-hoc 'penetration electric field'
t_PEF_start=0;
t_PEF_end=400;

%% Set up the data saving nc file
SaveDir=pwd;
SimID=1;
fileDir=dir([SaveDir,'\*TPS_VirtualSatellite_',num2str(SimID),'*.nc']);
Datafilename=[SaveDir,'\Data_TPS_VirtualSatellite_',num2str(SimID),'.nc'];
if isempty(fileDir)
    nccreate(Datafilename,'ID','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'T','Dimensions',{'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'R0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'phi0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'Ek0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'P_R','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'P_phi','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'Ek','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
else
    delete ([SaveDir,'\',fileDir(1).name])
    nccreate(Datafilename,'ID','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'T','Dimensions',{'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'R0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'phi0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'Ek0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'P_R','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'P_phi','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
    nccreate(Datafilename,'Ek','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable','DeflateLevel',0)
end
ncwrite(Datafilename,'ID',[1:e_num_Ek*n_Ek0])
ncwrite(Datafilename,'T',Time_sim)
ncwrite(Datafilename,'R0',repmat(reshape(Posit_R_ini,[],1),n_Ek0,1))
ncwrite(Datafilename,'phi0',repmat(reshape(Posit_phi_ini,[],1),n_Ek0,1))
ncwrite(Datafilename,'Ek0',reshape(Ek_ini,[],1))

%% Motion of particles

temp=find(Time_sim>rec_time_1 & Time_sim<rec_time_2);
ind_rec1=temp(1);
ind_rec2=temp(end);

phi_B_temp=repmat(phi_B(1:end),1,length(R_B));
R_B_temp=repmat(R_B,length(phi_B),1);
B_mat=B_total;
[~,Vradial_E]=V_Drift_E(E_phi,B_mat);
[~,Vphi_E]=V_Drift_E(E_R,B_mat);
for m_E=1:length(Ek_TPS)
    tic
    Posit_R=zeros(e_num_Ek,nT+1);
    Posit_phi=zeros(e_num_Ek,nT+1);

    Posit_R(:,1)=Posit_R_ini(:);
    Posit_phi(:,1)=repmat(phi_ini,length(R_aim),1);
    Energy_MeV=ones(e_num_Ek,1)*Ek_TPS(m_E);
    tao_count=zeros(size(Posit_R(:,1)));

    for m1=1:nT
        alpha_temp=interp1(L_obs,alpha_obs,Posit_R(:,m1));
        [tao_D_temp,V_phi1,~]=V_Drift_GradCurv_Relativity_Roederer(alpha_temp,Posit_R(:,m1),Energy_MeV(:,m1)); % V_phi here is omega_phi
        if Time_sim(m1)>=t_PEF_start && Time_sim(m1)<=t_PEF_end
            % when the penetration E-field is on
            V_R0=interp2(R_B_temp,phi_B_temp,Vradial_E,Posit_R(:,m1),Posit_phi(:,m1));
            V_R=V_R0;
            V_phi_E2=interp2(R_B_temp,phi_B_temp,Vphi_E,Posit_R(:,m1),Posit_phi(:,m1));
            V_phi2=V_phi_E2./Posit_R(:,m1);
            V_phi=V_phi1+V_phi2;
        else
            % when the E-field is off
            V_R=0;
            V_phi=V_phi1;
        end

        Posit_phi(:,m1+1)=Posit_phi(:,m1)+V_phi*time_step;
        Posit_phi(Posit_phi(:,m1+1)>2*pi,m1+1)=Posit_phi(Posit_phi(:,m1+1)>2*pi,m1+1)-2*pi;
        Posit_phi(Posit_phi(:,m1+1)<0,m1+1)=Posit_phi(Posit_phi(:,m1+1)<0,m1+1)+2*pi;
        Posit_R(:,m1+1)=Posit_R(:,m1)+V_R*time_step;
        tao_count=tao_count+time_step;
        tao_count_temp=mod(tao_count,tao_D_temp);

        Energy_MeV(:,m1+1)=Energy_MeV(:,m1).*((Posit_R(:,m1)./Posit_R(:,m1+1)).^3);
        tao_count=tao_count_temp;
    end

    ncwrite(Datafilename,'P_R',Posit_R(:,ind_rec1:ind_rec2),[e_num_Ek*(m_E-1)+1 ind_rec1])
    ncwrite(Datafilename,'P_phi',Posit_phi(:,ind_rec1:ind_rec2),[e_num_Ek*(m_E-1)+1 ind_rec1])
    ncwrite(Datafilename,'Ek',Energy_MeV(:,ind_rec1:ind_rec2),[e_num_Ek*(m_E-1)+1 ind_rec1])

    timer_sim=toc;
    ['Simulation progress: ',num2str(m_E/length(Ek_TPS)*100),'%. Elapsed time: ',num2str(timer_sim),' seconds']
    clear Posit_phi Posit_R Ek tao_D_temp
end

%% Plot the observed electron spectrogram
fig=figure;
fig.Position=[346 75 1400 900];
Flux_plot=Flux_obs;
Flux_plot(Flux_plot<10)=nan;
h1=pcolor(Time_obs,eleCh,log10(Flux_plot));
set(h1,'linestyle','none')
set(gca, 'YScale', 'log')
hold on
temp1=Time_obs(1):1/24/60*3:Time_obs(end)+1/24/60*3;
xticks(temp1)
xticklabels(cellstr(datestr(temp1,'HH:MM')))
xlim([Time_obs(1) Tsim_end])
ylim([eleCh(1) eleCh(end)])
ylim([eleCh(1) 2])
colormap('jet')
c_CB=colorbar;
c_CB.Label.String = 'log_{10}(flux)(cm^{-2}sr^{-1}s^{-1}MeV^{-1})';
c_CB.Label.FontSize=30;
c_CB.TickDirection='out';
c_CB.TickLength=0.018;
c_CB.FontSize=25;
clim([2 6])
yticks([0.2 0.5 1 1.5 2 3 4 5])
grid on

hAX = gca();
ylabel(hAX,'Energy (MeV)')
xlabel(hAX,'Time')
title('Observation')
hAX.FontSize=27;
hAX.LineWidth=2;
set(hAX,'FontWeight','bold')
set(hAX,'XMinorTick','on','YMinorTick','on')
set(hAX,'TickLength',[0.01, 0.005])
set(hAX,'TickDir','out')


%% get smoothed flux spectrogram as initial condition
temp=log10(Flux_obs);
temp(isinf(temp))=nan;
temp2=movmean(temp,7); % moving average over energy (+/- 3 channels)
Flux_ini_smooth=10.^temp2;
%% set up virtual satellite and plot the simulation results
Phi_1=MLT_obs/24*2*pi;
ind_plot=ind_rec1;
temp=find(Time_sim>time_sim(end));
ind_T_plot=ind_plot:1:temp(1)-1;
T_plot=Time_sim(ind_T_plot); % sec
J_Ek_plot=zeros(length(T_plot),length(eleCh));
L_plot= interp1(time_sim,L_obs,T_plot,'linear','extrap');
Phi_plot=interp1(time_sim,Phi_1,T_plot);
Ek_ini=repmat(Ek_TPS,e_num_Ek,1);
Ek_ini=Ek_ini(1:length(Ek_TPS)*e_num_Ek)';
Flux_smooth_ext=Flux_ini_smooth;
Flux_smooth_ext(isnan(Flux_ini_smooth)&(~isnan(Flux_obs)))=Flux_obs(isnan(Flux_ini_smooth)&(~isnan(Flux_obs)));
Read_count=[inf 1];
Read_stride=[1 1];
Q_ini=interp2(sort(eleCh),sort(L_obs),Flux_smooth_ext',Ek_TPS,R_aim); % initial weight of test particles
ind_Eklow=find(Ek_TPS<eleCh(1));
ind_Ekhigh=find(Ek_TPS>=eleCh(1));
% truncate particle energies by the instrument energy range
for m1=1:length(R_aim)
    Qtemp=interp1(Ek_TPS(ind_Ekhigh),log10(Q_ini(m1,ind_Ekhigh)),Ek_TPS(ind_Eklow),'linear','extrap');
    Q_ini(m1,ind_Eklow)=10.^(Qtemp);
end
for n1=1:length(T_plot)
    Read_start=[1 ind_plot+n1-1];
    P_R=ncread(Datafilename,'P_R',Read_start,Read_count,Read_stride);
    P_phi=ncread(Datafilename,'P_phi',Read_start,Read_count,Read_stride);
    Ek_sim=ncread(Datafilename,'Ek',Read_start,Read_count,Read_stride);
    Meas_phi_edge=[Phi_plot(n1)-1/2*dPhi_m Phi_plot(n1)+1/2*dPhi_m];
    Meas_R_edge=[L_plot(n1)-dR/2 L_plot(n1)+dR/2];
    ind_R=find(P_R(:,1)>Meas_R_edge(1) & P_R(:,1)<=Meas_R_edge(2));
    Meas_temp=zeros(size(eleCh));
    for n2=1:50
        Meas_Ek_edge=[RNGe.Ebeg(n2) RNGe.Eend(n2)];
        Meas_temp(n2)=histcounts2(P_phi(ind_R,1),Ek_sim(ind_R,1),Meas_phi_edge,Meas_Ek_edge);
    end
    ind_phi=find(P_phi(ind_R,1)>Meas_phi_edge(1) & P_phi(ind_R,1)<=Meas_phi_edge(2));
    ind_combined=ind_R(ind_phi);
    temp=repmat(Posit_R_ini(:),length(Ek_TPS),1);
    R_ini=temp(ind_combined);
    R_ini_edge=1/2*([R_aim(1)-0.05; R_aim]+[R_aim; R_aim(end)+0.1]);
    Q_edge=[Ek_TPS-1/2*dEk Ek_TPS(end)+1/2*dEk];
    Q_temp=zeros(length(Ek_TPS),length(R_aim),50);
    for n2=1:50
        Meas_Ek_edge=[RNGe.Ebeg(n2) RNGe.Eend(n2)];
        Q_temp0=histcnd(Ek_ini(ind_combined),R_ini,Ek_sim(ind_combined,1),Q_edge,R_ini_edge,Meas_Ek_edge);
        Q_temp(:,:,n2)=Q_temp0(1:end-1,1:end-1,1);
    end
    clear Q_temp0
    temp=Q_ini'.*Q_temp;
    temp(isnan(temp))=0;
    J_Ek_plot(n1,:)=sum(sum(temp,1),2)./reshape(Meas_temp,1,1,[]);
end

T_plot2=T_plot/24/3600+Tsim_initial;

fig=figure;
fig.Position=[346 75 1400 900];
Flux_plot=J_Ek_plot;
Flux_plot(Flux_plot<10)=nan;
h1=pcolor(T_plot2,eleCh,log10(Flux_plot'));
set(h1,'linestyle','none')
set(gca, 'YScale', 'log')
hold on
ylim([eleCh(1) 2])
temp1=Time_obs(1):1/24/60*3:Time_obs(end)+1/24/60*3;
xticks(temp1)
xticklabels(cellstr(datestr(temp1,'HH:MM')))
xlim([Time_obs(1) Tsim_end])
colormap('jet')
c_CB=colorbar;
c_CB.Label.String = 'log_{10}(flux)(cm^{-2}sr^{-1}s^{-1}MeV^{-1})';
c_CB.Label.FontSize=30;
c_CB.TickDirection='out';
c_CB.TickLength=0.018;
c_CB.FontSize=25;
clim([2 6])

temp0=get(get(gca,'title'),'Position');
temp0(2)=111.5;
set(get(gca,'title'),'Position',temp0)
yticks([0.2 0.5 1 1.5 2 3 4 5])

grid on
hAX = gca();
ylabel(hAX,'Energy (MeV)')
xlabel(hAX,'Time')
title('Simulation')
hAX.FontSize=27;
hAX.LineWidth=2;
set(hAX,'FontWeight','bold')
set(hAX,'XMinorTick','on','YMinorTick','on')
set(hAX,'TickLength',[0.01, 0.005])
set(hAX,'TickDir','out')

end


%% functions
function [V_m,V_L]=V_Drift_E(E,B)
% E:mV/m; B:nT

R_E=6371e3; % meter
V_m=E./B*1e6; % m/s
V_L=V_m/R_E; % size: phi*L
end

function [tao_D,omega_D,V_d]=V_Drift_GradCurv_Relativity_Roederer(alpha_e,L,Ek_MeV)
% Based on (Roederer 1970; Schulz and Lauzerottim 1974)
R_E=6371e3; % meter
B0=3.11e4;  % unit: nT
E0=0.51099895; % MeV
gama=Ek_MeV/E0+1;
beta_square=1-1./gama.^2;
T_ae=1.30-0.56*sin(alpha_e);
D_ae=1/12*(5.520692-2.357194*sin(alpha_e)+1.279385*(sin(alpha_e)).^(3/4));
Rc=R_E*L;

omega_D=3*E0.*gama.*beta_square.*L/(B0*R_E^2)./T_ae.*D_ae*1e15;
tao_D=2*pi./omega_D; % unit sec
V_d=omega_D.*Rc; % unit: m/s
end

function B=Geomagnetic_model_dipole_eq(L)
B0=3.11e4; % equatorial magnetic field on the Earth's surface, unit: nT
lamda_pi=0;
a1=(1+3*sin(lamda_pi).^2)^0.5;
b1=cos(lamda_pi).^6;
B=B0*a1./(L.^3*b1);
end
