%% This main file runs the ODE Herrera_2023_model_SK_Expression
clear
close all
clc

%% Preconditions
AF = 0; % 0 = nSR; 1 = AF
SK_expression = 1; % Change this to either [0.5 1 2] = 50%, 100%, 200% SK Expression
prot_rate = 1; % prot_rate = pacing frequency

%% Load Initial Conditions 
if AF == 0
    load yfin_nSR_1Hz_1SK.mat % nSR 100% GSK (1 Hz)
elseif AF == 1
    load yfin_AF_1Hz_1SK.mat % cAF 100% GSK (1 Hz)
end

%% Simulation parameters
% Stimulation protocol parameters
prot_index = 1; % Selection of the protocol (-)

prot_interval = 1000;

% Drug parameters
drug_index = 0; drug_conc = 0; % Drug Free

% Other experimental conditions
exp_Temp = 310; % [K]
exp_Nao = 140; % [Na]o
exp_ISO = 0; %; % (boolean)
exp_Ach = 0; %; % (boolean)   % no IKach

% IKur properties
% 0 drug free; 1 drug-CLOSED; 2 drug-OPEN; 3 drug-INACTIVATED; 4 drug-OPEN
% & INACTIVE
drug_index_ikur = 0; % (-)
Constant_Block = 0; % 1 if consant 50% block, 0 if not
Complete_Block = 0;
drug_conc_kur = 0; % mM
kon = 0.00001; % 1/ms
koff = 0.00001;  % 1/ms

%% Run Simulation
duration = 3e3;  % (ms)
period = 1000/prot_rate;
num_beats = floor(duration/period);
duration = (num_beats*period);
tspan = [0; duration];
options = odeset('RelTol',1e-6,'MaxStep',1);

y0 = yfinal;
prot_vm = yfinal(39); % (mV)
prot_par = [prot_index prot_rate prot_interval prot_vm];    % 1 2 3 4
drug_par = [drug_index drug_conc];                          % 5 6
exp_par = [exp_Temp exp_Nao exp_ISO exp_Ach];               % 7 8 9 10
ikur_par = [drug_index_ikur Constant_Block drug_conc_kur kon koff Complete_Block]; % 11 12 13 14 15 16
par_SA = ones(1, 45);
p = [prot_par, drug_par, exp_par, ikur_par, AF, SK_expression par_SA]; %1 p.option

[t,y] = ode15s(@Herrera_2023_model_SK_expression,tspan,y0,options,p);

time = t; % (ms)
Vm = y(:,39); % (mV)
Ca = y(:,38); % (mM)
Na = y(:,34); % (mM)

yfinal = y(end,:);
currents = calcCurrents(t,y,p);

%% All Figures
% AP
figure,set(gcf,'color','w')
subplot(3,1,1)
hold on,
plot(t,y(:,39));
ax = gca;
ax.TickDir = 'out';
ylabel(('Em (mV)'))

% CaT total
subplot(3,1,2)
hold on,
plot(t,y(:,38)*1e6);
ylabel(('[Ca2+]i (nM)'))
ax = gca;
ax.TickDir = 'out';
box off

% ISK total
subplot(3,1,3)
hold on,
plot(t,currents(:,6));
ylabel(('ISK (pA/pF)'))
ax = gca;
ax.TickDir = 'out';
xlabel('Time (ms)')

%% Biomarker Analysis
dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
dVm = [dVm_array(1); dVm_array];
AP_index = 2; % w/ 1 first AP, otherwise last AP
%     outputs = [dVm_max Vm_max -Vm_min AP_amp APD90 APD70 APD50 APD30...
%         Ca_max Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min];
period = 1000/prot_rate;

outputs = function_beat_analysis(time,Vm,Ca,Na,dVm,period,AP_index);
dVm_max = outputs(1);
Vm_max = outputs(2);
RMP = outputs(3);
AP_amp = outputs(4);
APD90 = outputs(5);
APD70 = outputs(6);
APD50 = outputs(7);
APD30 = outputs(8);
Ca_max = outputs(9);
Ca_min = outputs(10);
CaT_amp = outputs(11);
CaT_rise = outputs(12);
CaT_decay_50 = outputs(13);
CaT_decay_63 = outputs(14);
Na_min = outputs(15);

output = [dVm_max Vm_max RMP AP_amp APD90 APD70 APD50 APD30 Ca_max Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min];
