clear
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%"!" means the default value
eco_pars.Nsupply = 5;    %%VALUES: !3                                           %Nutrient concentration
eco_pars.tau     = 1;    %%VALUES: 0.2 0.4 0.6 0.8 !1 1/0.8 1/0.6 1/0.4 1/0.2   %Trophic trade-off parameter
eco_pars.T       = 20;   %%VALUES: 0 5 10 15 !20 25 30 35 40 45                 %Temperature
eco_pars.gamma   = 0.1;  %%VALUES: 0.2:0.15:1                                   %Growth limitation term
eco_pars.tempfuA = 0.05; %%VALUES: 0.01 !0.05 0.1 0.5 1                         %Temperature Afunction parameter
eco_pars.tempfuZ = 0.05; %%VALUES: 0.01 !0.05 0.1 0.5 1                         %Temperature Zfunction parameter
eco_pars.kappa   = 0.01; %%VALUES: !0.01                                        %Nutrient Flow

%% MEGALOOP
%% PARAMETER LOOP
par.Nsupply = [5];
par.tau     = [1];
par.T       = [20];
par.gamma   = [0.1];
par.tempfuA = [0.05];
par.tempfuZ = [0.05];
par.kappa   = [0.01];

c=0;
fig = 2;
for taucount = 1:length(par.tau)
    for tempcount = 1:length(par.T)
        for gammacount = 1:length(par.gamma)
            clearvars -except par eco_pars c Tcount gammacount taucount tempcount fig ext
            
            % count iterations
            if c == length(par.T)*length(par.gamma)
                c = 1;
            else
                c = c+1;
            end
            eco_pars.tau = par.tau(taucount);
            eco_pars.T = par.T(tempcount);
            eco_pars.gamma = par.gamma(gammacount);
            
MCM_initialise           %%INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eco_pars.nsize      = 11;
%eco_pars.ntroph     = 11;
%eco_pars.jpmax=eco_pars.nsize.*eco_pars.ntroph;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define simulation conditions
save_iterations = true;
years_per_iteration = 1;
eco_pars.seasonalcycle = false;  %Seasonal Cycle
EiE = true;                      %"Everything is Everywhere" option
SS = false;                      %Load Steady State
nyears = 1000;                   %50 YEARS FOR EiE

%% Extinction parameters (works onyl when SS is true)
%KILL 1=AUTOTROPHS, 2=MIXOTROPHS, 3=HETEROPHS, 4=AUT-MIX, 5=HETER-MIX, 6=AUT-HETER, 7=NONE, 8=AUT-MIX-HETER
ext = 7;
%Phenotypic Space indexes to kill
% Try to avoid entering numbers when the same job could be done more flexibly with parameters
PSindex = 1:eco_pars.jpmax; % PSindex = 1:441;
PSindex = reshape(PSindex,[eco_pars.nsize,eco_pars.ntroph]); % PSindex = reshape(PSindex,[21,21]);
PSindex = PSindex';

%% Seed Change for Everything is Everywhere option
if EiE == true
    iseed = [];
else [iseed] = Seedchange(eco_pars,PSindex);
end

%% Initial conditions
N_0 = eco_pars.Nsupply;
P_0 = zeros(eco_pars.jpmax,1);
ndays_per_iteration=365*years_per_iteration; % needs to be an exact multiple of 365
n_ode_iterations=nyears;%./ndays_per_iteration*365;
tf=ndays_per_iteration;
eco_pars.t_res=7;
t0=0;
ndata_yr=numel([t0:eco_pars.t_res:tf]);
data=zeros(ndata_yr.*n_ode_iterations,eco_pars.jpmax+1);
data=sparse(data);
if isempty(iseed)
    P_0(:) = eco_pars.seed_val; % seed  everything
    N_0    = N_0 - sum(P_0); 
    eco_pars.Pmut=speye(size(eco_pars.Pmut));
    minphy=eco_pars.extnct;
    eco_pars.minphy = minphy; % set minimum 'everything is everywhere' threshold
    eco_pars.extnct = eco_pars.extnct ./ 1000; % decrease 'extinction' threshold
else
    P_0(iseed) = eco_pars.seed_val;
    minphy=0;
    eco_pars.minphy = minphy;
end
v0=[N_0;P_0];
dead=[];

%% Load the Steady State
if SS == true
    [v0] = Extinction(PSindex,ext);
end

%% DIRECTORY NAMES
% CA = COMMUNITY ASSEMBLY
% CP = COMMUNITY PERTURBATION
if SS == false
    state='_CA';
    else state='_CP';
end
% "Everything is Everywhere" option
if EiE == true
    EiE = '_EiE';
    else EiE = '';
end
% Seasonalcycle, Y = Yes, N = No
if eco_pars.seasonalcycle == true
    sc = '_Seasons_';
    else sc = '';
end
% Extinction
d = ["_AUT","_MIX","_HETER","_AUT-MIX","_HETER-MIX","_AUT-HETER","","_AUT-MIX-HETER"];
if ext == 1
    ex=d(1);
end
if ext == 2
    ex=d(2);
end
if ext == 3
    ex=d(3);
end
if ext == 4
    ex=d(4);
end
if ext == 5
    ex=d(5);
end
if ext == 6
    ex=d(6);
end
if ext == 7
    ex=d(7);
end
if ext == 8
    ex=d(8);
end
ex = convertStringsToChars(ex);
% Name directory of results
% Example: "Run_CA_N3_Tau1_T20_2021_02_29"
%
% '_' datestr(now,'yyyy_mm_dd')
Output_fdir = ['Run' state '_N' num2str(eco_pars.Nsupply) '_Tau' num2str(eco_pars.tau) ...
               EiE sc ex '_T' num2str(eco_pars.T) ...  
               '_tfuA' num2str(eco_pars.tempfuA) '_tfuZ' num2str(eco_pars.tempfuZ) ...
               '_gamma' num2str(eco_pars.gamma) ];
if exist(['Output/' Output_fdir])==0
    mkdir(['Output/' Output_fdir]);
    mkdir(['Output/' Output_fdir '/Figures']);
end

%% Solve!
tic
for k=1:n_ode_iterations
    eco_pars.k = k;
    eco_pars.nyears = nyears;
    eco_pars.path = ['Output/' Output_fdir '/Network1.mat'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%ODE%%%%%%%%%%%%%%%%%%%%%%%%%  
    % call ecosystem function
    [tout,yout] = ode45(@(t,y)  ecosystem(t,y,eco_pars,dead,minphy) ,  [t0:eco_pars.t_res:tf],        v0);
                        %       control on the biomass                 control on the time(375/53)   initial state

    data((k-1).*ndata_yr+1:k.*ndata_yr,:)=yout;
    
    t0=tout(end);
    tf=t0+ndays_per_iteration;
    v0=yout(end,:);
    

    %%%%%%%%%%%%%%%%%NETWORK DATA%%%%%%%%%%%%%%%%%%%%
     if ismember(k,[1:nyears+1])
        path = ['Output/' Output_fdir '/steadystate.mat'];
        seesteady(k,:) = yout(53,:);
        save(path, 'seesteady')
     end

    %%%%%%%%%%%%%%%MODEL PIC YEAR %%%%%%%%%%%%%%%%%%%
    if ismember(k,[10 20 100 200 500 700 1000 1500 2000 2500 3000 3500 4000 4500 5000 6000 7000 8000 9000 10000])
        fh=figure(1);
        fh.Position=[103 1309 846 902];
        disp_indx=numel(tout);
        plot_output(tout,yout,eco_pars,disp_indx)    %shows the final result of the model
        fname=['Output/' Output_fdir '/Figures/FitLand_' num2str(k,'%05i') '.png'];
        saveas(fh,fname);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
    %%%%%%%%%%%Display every 10 years%%%%%%%%%%%%%%%%%
    if rem(k,10)==0
        disp(['Year = ' num2str(k.*ndays_per_iteration./365) '; ' num2str(toc) ' seconds.'])
        tic
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

save(path, 'seesteady')
save(['Output/' Output_fdir '/Workspace_dump.mat'],'data','eco_pars','nyears','years_per_iteration','iseed','dead','minphy') %MODIFY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEADY STATE PIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg=figure(2);
loglog(seesteady)
if save_iterations
    fname=['Output/' Output_fdir '/Figures/steadystate_' num2str(k,'%05i') '.png'];             %MODIFY
    saveas(gcf,fname);
end

%% END OF THE PARAMETER LOOP
        end
    end
end
