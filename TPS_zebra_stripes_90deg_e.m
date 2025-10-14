%% ========================================================================
% Script Name: TPS_zebra_stripes_90deg_e.m
% Description:
%   This script can show how electron zebra stripes will show up in
%   electron flux spectrogram in the presence of an ad-hoc electric field
%
%   Main features:
%     - Computes electron drift trajectories in the presence of static
%     electric fields
%     - Test particle trajectories are saved as .nc files
%     - Generate plots to illustrate the formation of zebra stripes in
%     electron flux spectrogram
%
% Author: Yang Mei
% Affiliation: Laboratory for Atmospheric and Space Physics, University of
% Colorado, Boulder
% Contact: Yang.Mei@colorado.edu
% GitHub: https://github.com/yangmei-spacephysics/Test-Particle-Simulation_90deg-electron_zebra-stripes
%
% Version: 1.0
% Date: 2025-10-13
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
%
% ============================================================


clear
clc
% addpath('D:\Research\Function by Mei')
%% Simulation set up
% Simulation time duration and time step
time_step=50; % unit: seconds
tsim_end=10000; % sec
Time_sim=0:time_step:tsim_end;
nT=length(Time_sim)-1;
len_T=nT+1;


%% particle set up
dR=0.1; % Unit: Earth radii (R_E)
R_aim=[1.1:dR:5]'; % The radial distance of the electrons at the beginning
e_num0=1000; % numbers of the electrons at a L for one energy

phi_step=2*pi/e_num0;
phi_ini0=[0:phi_step:2*pi-phi_step]';
phi_ini=phi_ini0;
dEk=0.0125; % Unit: MeV
Ek_TPS=[0.05:dEk:2 2.05:0.1:2.3]; % Initial energy
n_Ek0=length(Ek_TPS);
e_num_Ek=length(R_aim)*(e_num0); % how many particles

Posit_phi_ini=repmat(phi_ini,1,length(R_aim));
Posit_R_ini=repmat(R_aim',length(phi_ini),1);
Ek_ini=ones(e_num_Ek,1)*Ek_TPS;

%% Virtual measurements set up
dPhi_m=pi/5; % delta azimuthal angle where measurements are made
Phi_m_mean=pi/2; % centered azimuthal angle of the virtual measurements
Phi_m_lo=Phi_m_mean-dPhi_m/2;
Phi_m_hi=Phi_m_mean+dPhi_m/2;
if Phi_m_lo<0
    Phi_m_lo=Phi_m_lo+2*pi;
end
R_m_mean=R_aim;
dR_m=0.1;
N_m=zeros(length(R_m_mean),length(Time_sim));

%% Weighting of the particle
% no weighting factor over L
fr=ones(size(R_aim));
% Assuming a kappa distribution as the initial energy spectrum
E0=0.5*1e-3; % keV
k_Ek=2.9;
f0=1e16;
f_Ek0=f0*(1+Ek_TPS/((k_Ek-3/2)*E0)).^(-k_Ek-1); % Unit: (cm^2*s*sr*MeV)^-1

%% Energy channel set up
% setting up energy channel
dEk2=0.05;
Ek_Meas=0.2:dEk2:2; % MeV
%% Field set up
phi_step_B=phi_step/5;
phi_B(:,1)=0:phi_step_B:2*pi;
phi_E=phi_B-phi_step_B/2;
R_step=dR;
R_B=1.1:R_step:8;

B0=Geomagnetic_model_dipole_eq(R_B);
B_total=repmat(B0,length(phi_B),1);

%% Set up an ad-hoc 'penetration electric field'
t_PEF_start=0;
t_PEF_end=400;

% add radial dependence
Ey0=zeros(size(R_B));
Ey0(R_B<=4)=2*R_B(R_B<=4)-10;
Ey0(R_B>4)=-2;
Ex0=0;
E_phi=cos(phi_B)*Ey0-sin(phi_B)*Ex0*ones(size(R_B));
E_R=sin(phi_B)*Ey0+cos(phi_B)*Ex0*ones(size(R_B));

E_y=E_phi.*cos(phi_B)+E_R.*sin(phi_B);
E_x=-E_phi.*sin(phi_B)+E_R.*cos(phi_B);
%% Plotting the ad-hoc electric field
fig=figure;
fig.Position=[299 165 1226 781];
plot(R_B,abs(E_y),'linewidth',2)
xlim([1 6])
ylim([0 10])
xlabel('L')
ylabel('E (mV/m)')
title(['Electric field from 0-',num2str(t_PEF_end),' sec'])
set(gca,'LineWidth',2);
set(gca,'FontSize',25);
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLength',[0.01, 0.005])
set(gca,'FontWeight','bold')
%% Set up the data saving nc file
SaveDir=pwd;
SimID=1;
fileDir=dir([SaveDir,'\*TPS_',num2str(SimID),'*.nc']);
Datafilename=[SaveDir,'\Data_TPS_',num2str(SimID),'.nc'];
if isempty(fileDir)
    nccreate(Datafilename,'ID','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'T','Dimensions',{'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'R0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'phi0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'Ek0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'P_R','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'P_phi','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'Ek','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
else
    delete ([SaveDir,'\',fileDir(1).name])
    nccreate(Datafilename,'ID','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'T','Dimensions',{'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'R0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'phi0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'Ek0','Dimensions',{'ID',e_num_Ek*n_Ek0},...
        'FillValue','disable')
    nccreate(Datafilename,'P_R','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'P_phi','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
    nccreate(Datafilename,'Ek','Dimensions',{'ID',e_num_Ek*n_Ek0,'T',len_T},...
        'FillValue','disable')
end
ncwrite(Datafilename,'ID',[1:e_num_Ek*n_Ek0])
ncwrite(Datafilename,'T',Time_sim)
ncwrite(Datafilename,'R0',repmat(reshape(Posit_R_ini,[],1),n_Ek0,1))
ncwrite(Datafilename,'phi0',repmat(reshape(Posit_phi_ini,[],1),n_Ek0,1))
ncwrite(Datafilename,'Ek0',reshape(Ek_ini,[],1))

%% Motion of particles
Ek_all=cell(length(Ek_TPS),1);
for m_E=1:length(Ek_TPS)
    tic
    Posit_R=zeros(e_num_Ek,nT+1);
    Posit_phi=zeros(e_num_Ek,nT+1);
    
    Posit_R(:,1)=Posit_R_ini(:);
    Posit_phi(:,1)=repmat(phi_ini,length(R_aim),1);
    Energy_MeV0=ones(e_num_Ek,1)*Ek_TPS(m_E);
    Energy_MeV=zeros(length(Energy_MeV0),nT+1);
    Energy_MeV(:,1)=Energy_MeV0;

    phi_B_temp=repmat(phi_B(1:end),1,length(R_B));
    R_B_temp=repmat(R_B,length(phi_B),1);
    tao_count=zeros(size(Posit_R(:,1)));
    k_0=0;
    B_mat=B_total;
    [~,Vradial_E]=V_Drift_E(E_phi,B_mat);
    [~,Vphi_E]=V_Drift_E(E_R,B_mat);
    for m1=1:nT
        % penetration electric field
        [tao_D_temp,V_phi1,~]=V_Drift_GradCurv_Relativity_Roederer(pi/2,Posit_R(:,m1),Energy_MeV(:,m1)); % V_phi here is actually omega_phi
       
        if Time_sim(m1)>=t_PEF_start && Time_sim(m1)<=t_PEF_end
            V_R0=interp2(R_B_temp,phi_B_temp,Vradial_E,Posit_R(:,m1),Posit_phi(:,m1));
            V_R=V_R0;
            V_phi_E2=interp2(R_B_temp,phi_B_temp,Vphi_E,Posit_R(:,m1),Posit_phi(:,m1));
            V_phi2=V_phi_E2./Posit_R(:,m1);
            V_phi=V_phi1+V_phi2;
        else
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
    
    ncwrite(Datafilename,'P_R',Posit_R,[e_num_Ek*(m_E-1)+1 1])
    ncwrite(Datafilename,'P_phi',Posit_phi,[e_num_Ek*(m_E-1)+1 1])
    ncwrite(Datafilename,'Ek',Energy_MeV,[e_num_Ek*(m_E-1)+1 1])

    timer_sim=toc;
    ['Simulation progress: ',num2str(m_E/length(Ek_TPS)*100),'%. Elapsed time: ',num2str(timer_sim),' seconds']
    clear Posit_phi Posit_R Ek tao_D_temp
end
%% plot the particle spectrogram
d_T_plot=1; % time resolution
ind_T_plot=1:d_T_plot:nT+1;
T_plot=Time_sim(ind_T_plot); % sec

J_Ek_plot=zeros(length(T_plot),length(Ek_Meas));
L_aim=[1.5 4];
dL=0.1;
Ek_ini=repmat(Ek_TPS,e_num_Ek,1);
Ek_ini=Ek_ini(1:length(Ek_TPS)*e_num_Ek)';

for n1=1:length(L_aim)
    for n3=1:length(ind_T_plot)
        Read_start=[1 n3];
        Read_count=[inf 1];
        Read_stride=[1 1];
        P_R=ncread(Datafilename,'P_R',Read_start,Read_count,Read_stride);
        P_phi=ncread(Datafilename,'P_phi',Read_start,Read_count,Read_stride);
        Ek_sim=ncread(Datafilename,'Ek',Read_start,Read_count,Read_stride);
        Meas_phi_edge=[Phi_m_lo Phi_m_hi];
        Meas_R_edge=[L_aim(n1)-dR/2 L_aim(n1)+dR/2];
        ind_R=find(P_R(:,1)>Meas_R_edge(1) & P_R(:,1)<=Meas_R_edge(2));
        Meas_Ek_edge=[Ek_Meas(1)-1/2*(Ek_Meas(2)-Ek_Meas(1)) (Ek_Meas(1:end-1)+Ek_Meas(2:end))/2  Ek_Meas(end)+1/2*(Ek_Meas(end)-Ek_Meas(end-1))];
        Meas_temp=histcounts2(P_phi(ind_R,1),Ek_sim(ind_R,1),Meas_phi_edge,Meas_Ek_edge);
        ind_phi=find(P_phi(ind_R,1)>Meas_phi_edge(1) & P_phi(ind_R,1)<=Meas_phi_edge(2));
        ind_combined=ind_R(ind_phi);
        R_ini=P_R(ind_combined,1);
        R_ini_edge=[R_aim-1/2*dR;R_aim(end)+1/2*dR];
        Q_edge=[Ek_TPS-1/2*dEk Ek_TPS(end)+1/2*dEk];
        Q_temp0=histcnd(Ek_ini(ind_combined),R_ini,Ek_sim(ind_combined,1),Q_edge,R_ini_edge,Meas_Ek_edge);
        Q_temp=Q_temp0(1:end-1,1:end-1,1:end-1);
        clear Q_temp0
        J_Ek_plot(n3,:)=sum(fr'.*sum(f_Ek0'.*Q_temp,1),2)./reshape(Meas_temp,1,1,[]);
    end
    
    fig=figure;
    fig.Position=[670 250 1100 600];
    h1=pcolor(T_plot,Ek_Meas,log10(J_Ek_plot'));
    set(h1,'linestyle','none')
    colormap('jet')
    c2=colorbar;
    c2.Label.String = 'log_{10}(flux)(cm^{-2}sr^{-1}s^{-1}MeV^{-1})';
    c2.Label.FontSize=18;
    clim([1 7])
    xlim([0 10000])
    ylim([0.2 1.6])
    title(['L=',num2str(L_aim(n1)),' MLT=',num2str(Phi_m_mean/2/pi*24)])
    xlabel('Time(sec)')
    ylabel('E_k (MeV)')
    set(gca,'FontSize',18)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'TickLength',[0.008, 0.0016])
    set(gca,'TickDir','out')
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