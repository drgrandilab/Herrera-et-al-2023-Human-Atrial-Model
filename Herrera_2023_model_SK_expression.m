function output = Herrera_2023_model_SK_expression(t,y,param,runType)

%% Parameters
p = param;

ydot = zeros(size(y));

par_SA = p(19:end);

% Perturbations used for populations
S_GNa   = par_SA(1);
S_GClCa = par_SA(2);
S_GCaL  = par_SA(3);
S_Gto   = par_SA(4);
S_GKur  = par_SA(5);
S_GKr   = par_SA(6);
S_GKs   = par_SA(7);
S_GK1   = par_SA(8);
S_GNCX  = par_SA(9);
S_GNaK  = par_SA(10);
S_GK2P  = par_SA(11);
S_GNaB  = par_SA(12);
S_GCaB  = par_SA(13);
S_GClB  = par_SA(14);
S_GKAch = par_SA(15);
S_GNaL = par_SA(16);
S_GKp = par_SA(17);
S_GSK = par_SA(18);
S_GCaP  = par_SA(19);
S_SERCA  = par_SA(20);
S_RyR  = par_SA(21);
S_RyR_leak  = par_SA(22);

%% IKur properties
IKur_M_flag = 1; % if 1 Markov model, 0 otherwise
IKur_ur_50 = p(12); % 1 if constant 50% IKur block, 0 if not
Complete_Block = p(16);
drug_kur_index = p(11); % (-)
drug_kur_conc = p(13); % (mM)
kon_kur =  p(14); % (1/(mM*ms))
koff_kur = p(15); % (1/ms)
% Kd_kur = koff/kon; % (mM)

%% Ca j and sl clamp
cajl_clamp=0;

if t>=2e3 && cajl_clamp==1
    tt=mod(t,1/(p(2)*1e-3));
    t_in_roi=find(p.altTime>=tt); t_in_index=t_in_roi(1);
    y(36) = p.altCaj(t_in_index);
    y(37) = p.altCasl(t_in_index);
end

%% Simulation Parameters
protocol_index = p(1);

% 1) 'pace_cc'; 2) 'pace_cc_erp'; 3) 'dte_finder'; 4) 'ap_clamp';
% 5) 'v_hold'; 6) 'ssa_prot'; 7) 'ssi_prot'; 8) 'rec_prot';
% 9) 'v_step_tb'; 10) 'v_step'; 11) 'v_step_interval'; 12) 'rec_im_prot';
if protocol_index == 1, protocol = 'pace_cc';
    elseif protocol_index == 2, protocol = 'pace_cc_erp';
    elseif protocol_index == 3, protocol = 'dte_finder';
    elseif protocol_index == 4, protocol = 'ap_clamp';
    elseif protocol_index == 5, protocol = 'v_hold';
    elseif protocol_index == 6, protocol = 'ssa_prot';
    elseif protocol_index == 7, protocol = 'ssi_prot';
    elseif protocol_index == 8, protocol = 'rec_prot';
    elseif protocol_index == 9, protocol = 'v_step_tb';
    elseif protocol_index == 10, protocol = 'v_step';
    elseif protocol_index == 11, protocol = 'v_step_interval';
    elseif protocol_index == 12, protocol = 'rec_im_prot';
    elseif protocol_index == 13, protocol = 'pace_cc_ead';
    elseif protocol_index == 14, protocol = 'pace_cc_dad';
    elseif protocol_index == 15, protocol = 'pace_cc_ead_sa';
    elseif protocol_index == 16, protocol = 'pace_cc_ead_2beat';
    elseif protocol_index == 17, protocol = 'ap_clamp_w_step';
    elseif protocol_index == 18, protocol = 'ap_clamp_ramp';
    elseif protocol_index == 19, protocol = 'tri_ap_clamp_wc';
    elseif protocol_index == 20, protocol = 'pace_hr_control';
end

pacing_rate = p(2); % (Hz) used only with 'pace_cc' and 'v_step'
rec_interval = p(3); % (ms) used only for "recovery" protocols (and ERP) (and 17)
vm_test = p(4);% (mV)

% Drug concentrations
drug_index = p(5); drug = p(6);

% Ach = p(10) * 0.1; % Acetylcholine (uM)
Ach = p(10) * 1; % Acetylcholine (uM)

%% Model Flags
% ISO
ISO = p(9);
if protocol_index == 15
    first_interval_duration = 30e3;
    if t < first_interval_duration
        ISO=0;
    else
        ISO=1;
    end
end

% AF
AF = p(17);

% SK Expression
SK_flag = p(18);

% Right ATRIUM
RA = 0;

% other
epi = 1; % EPI or ENDO?
bGal = 1;
EAD = 0;

% Model for INa
flagMina = 0; % flagMina = 1 (Markov INa) if flagMina = 0 (CM INa) 

%% Model Parameters
% Constants
Temp = p(7); % 310;     % [K]
R = 8314;       % [J/kmol*K]
Frdy = 96485;   % [C/mol]
FoRT = Frdy/R/Temp;
Cmem = 1.1e-10; % [F] membrane capacitance 1.3810e-10;
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um] 113;%100
cellRadius = 10.25;   % cell radius [um] 12;%10.25
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec]
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 1*0.0539*.01*Vcell;
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 9.9664e-014
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.7460e-012
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 6.6240e-013
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.5618e-011
% tau's from c-code]
J_ca_juncsl = 1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11;    Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = p(8); %140; % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT )*log(Cli/Clo);            % [mV]

%% Na transport parameters
GNa = 10*(1-0.1*AF);%*(1+0.2*EAD);  % [mS/uF] %Used for Markov; -Using CM instead
GNaB = S_GNaB*1*0.597e-3;    % [mS/uF] 
IbarNaK = S_GNaK*1*1.26;     % [uA/uF]
KmNaip = 11*(1-0.25*ISO);         % [mM]11
KmKo = 1.5;         % [mM]1.5
Q10NaK = 1.63;
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;
gkp = 0.002*S_GKp;

%% Cl current parameters
GClCa = S_GClCa* 0.0548;   % [mS/uF]
GClB = S_GClB *9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]
GClCFTR = 0;%4.9e-3*ISO;     % [mS/uF]

%% Ca transport parameters
% I_ca parameteres
pNa =S_GCaL*(1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       % [cm/sec]
pCa =S_GCaL*(1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       % [cm/sec]
pK =S_GCaL*(1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        % [cm/sec]
Q10CaL = 1.8;
% I_cabk parameteres
GCaB = S_GCaB*6.0643e-4;   % [uA/uF] 3
% NCX parameteres
IbarNCX = S_GNCX*(1+0.4*AF)*3.15;      % [uA/uF]5.5 before - 9 in rabbit

KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]
nu = 0.35;          % [none]
Kdact = 0.384e-3;   % [mM]
Q10NCX = 1.57;      % [none]

% I_pca parameteres
IbarSLCaP = S_GCaP*0.0471; % IbarSLCaP FEI changed [uA/uF] (2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa = 0.5e-3;     % [mM]
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3*S_SERCA;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = (2.5-1.5*ISO)*0.246e-3;          % [mM] default 2.5-1.25*ISO
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25*S_RyR;                 % [1/ms]
koCa = 10+20*AF+10*ISO*(1-AF);               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45*(1-0.5*AF);           % [mM] % NH added

%% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = (1+0.5*ISO)*19.6e-3;    % [1/ms]
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity
koff_tnchca = 0.032e-3; % [1/ms]
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms]
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM]  % CaM buffering
koff_cam = 238e-3;      % [1/ms]
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM]
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    % Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM]
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM]  %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms]
kon_csqn = 100;         % [1/mM/ms]

%% Membrane Currents
% Fast I_Na

% Courtemanche, Am J Physiol. 1998
scale = 5; % scaled to match Markov INaL
GNa_model_adjusment = 12; % original HH value is 23
ah = (y(39) >= -40.0) * ( 0.0 ) + (y(39) < -40.0) * ( 0.135 * exp( -( y(39) + 80.0 ) / 6.8 ) );
bh = (y(39) >= -40.0) * ( 1.0 / ( 0.13 * ( 1.0 + exp( -(y(39) + 10.66) / 11.1 ) ) ) ) + (y(39) < -40.0) * ((3.56 * exp( 0.079 * y(39)) + 3.1e5 * exp(0.35 * y(39))));
tauh = 1 / (ah + bh);
aj = (y(39) >= -40.0) * (0.0) +(y(39) < -40.0) * (( ( -127140 * exp(0.2444*y(39)) - 3.474e-5 * exp(-0.04391 * y(39))) * (y(39) + 37.78)) / (1.0 + exp( 0.311 * (y(39) + 79.23) ) ));
bj = (y(39) >= -40.0) * ((0.3 * exp(-2.535e-7*y(39))) / (1.0 + exp( -0.1 * (y(39) + 32.0) ))) + (y(39) < -40.0) * ((0.1212 * exp( -0.01052 * y(39) )) / (1.0 + exp( -0.1378 * (y(39) + 40.14) )));
tauj = 1 / (aj + bj);
am = (y(39) == -47.13) * ( 3.2 ) + (y(39) ~= -47.13) * ( 0.32 * (y(39) + 47.13) / (1.0 - exp( -0.1 * (y(39) + 47.13) ) ) );
bm = 0.08 * exp( -y(39) / 11.0 );
taum = 1 / (am + bm); % I added this__NH
mss = ( am/(am + bm) );
hss = ( ah/(ah + bh) );
jss = ( aj/(aj + bj) );

ydot(1) = (mss - y(1)) / taum;
ydot(2) = (hss - y(2)) / tauh;
ydot(3) = (jss - y(3)) / tauj;

% I_Na Conductance
GNa_CM = S_GNa*GNa_model_adjusment*(1-0.1*AF);  % [mS/uF] %CM- Adjusted maximal conductance 12/15/22

I_Na_junc1 = Fjunc*285/421*GNa_CM*y(1)^3*y(2)*y(3)*(y(39)-ena_junc); 
I_Na_sl1 = Fsl*285/421*GNa_CM*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);

% Late I_Na Conductance
GNaL = S_GNaL*0.0025*scale; % Adjusted to reproduce the same current density during AP -NH 12/15/22

% Late Activation
ydot(60) = 0; % Using y(1) from fast activation-gate

% Haibo Formulation Adjusted
hlinf = 1/(1+exp((y(39)+70)/6.1)); %(22 -> 37) %15 deg %V_shift -21 1/26/23 NH
tauhl = 600/2.83; %ms (divide by Q10 to go to 37deg) which means 600/~2.83
ydot(61) = (hlinf-y(61))/tauhl;


I_NaL_junc = Fjunc*GNaL*y(1)^3*y(61)*(y(39)-ena_junc); % Courtmanche Sodium Activation Gate 12/15/22
I_NaL_sl = Fsl*GNaL*y(1)^3*y(61)*(y(39)-ena_sl); % Courtmanche Sodium Activation Gate 12/15/22
I_NaL = I_NaL_junc + I_NaL_sl;
if t<9050
    ydot(62)=0;
else
    ydot(62)=I_NaL;
end

%% I_Na: Voltage-Gated Na Current (NEW MARKOV MODEL)
% Parameters
P1a1=3.802;
P2a1=0.1027;
P3a1=2.5;
P4a1=17;
P5a1=0.20;
P6a1=150;
P4a2=15;
P5a2=0.23;
P4a3=12;
P5a3=0.25;
P1b1=0.1917;
P2b1=20.3;
P1b2=0.2;
P2b2=2.5;
P1b3=0.22;
P2b3=7.5;
P1a4=0.188495;
P2a4=16.6;
P3a4=0.393956;
P4a4=7;
P1a5=7e-7;
P2a5=7.2; % TG 7.8
P1b5=0.0044; % TG 0.00467
P2b5=2e-5;
P1a6=100;
P1b6=8.9554e-7; % TG 6.5449e-7
P2b6=11.3944;
P1a7=0.487e-4; % TG 0.3377e-3
P2a7=23.2696;
P1b7=0.2868e-3; % TG 1.868e-4
P2b7=35.9898;
P1a8=0.1e-7; % TG 6.5e-6
P1b8=9.8e-3; % TG 3.8e-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (drug_index == 0 || drug_index == 1) % RAN (or Drug Free)
    diffusion = 5500;       % Ranolazine        %drug = 1 * (1E-6); % (M)
    pH = 7.4;
    pKa = 7.2;
    dd = -0.7;
    kd0 = 100.5 * (1e-6);
    kd0_b = 1.5012 * (1e-6); % bursting
    k_off_0 = 400 * (1e-6);
    ki_off_0 = 5.4 * (1e-6);
    kc_off_0 = 800 * (1e-6);
    Pa3_c = 3.6811;
    Pa4_c = 6.8705e+04;
    Pa5_c = 4.0832e-02;
    Pb5_c = 1.7561e-01;
    Pa6_c = 1*8;
    Pb6_c = 1/4;
    Pa7_c = 1;
    Pb7_c = 1;
    Pa3_n = 2.3570e+02;
    Pa4_n = 2.1182e+02;
    Pb5_n = 1.2197e-03;
    Pa6_n = 1;
    Pa7_n = 1;
elseif drug_index == 2     % Lidocaine         %drug = 1 * (1E-6); % (M)
    diffusion = 500;
    pH = 7.6;
    pKa = 7.2;
    dd = -0.7;
    kd0 = 318 * (1e-6);
    kd0_b = 1.5012 * (1e-6); % as in ranolazine model
    k_off_0 = 400 * (1e-6);
    ki_off_0 = 3.4 * (1e-6);
    kc_off_0 = 900 * (1e-6);
    Pa3_c = 5.6974e-03;
    Pa4_c = 6.7067e-06;
    Pa5_c = 3.2976;
    Pb5_c = 1.9698e-05;
    Pa6_c = 1;
    Pb6_c = 1;
    Pa7_c = 1;
    Pb7_c = 1;
    Pa3_n = 8.4559e+01;
    Pa4_n = 1.7084e-05;
    Pb5_n = 4.8477;
    Pa6_n = 1;
    Pa7_n = 1;
elseif drug_index == 3     % Flecainide        %drug = 1 * (1E-6); % (M)
    diffusion = 5500;
    pH = 7.4;
    pKa = 9.3;
    dd = -0.7;
    kd0 = 11.2 * (1e-6);
    kd0_b = 1.5012 * (1e-6); % as in ranolazine model
    k_off_0 = 400 * (1e-6);
    ki_off_0 = 5.4 * (1e-6);
    kc_off_0 = 800 * (1e-6);
    Pa3_c = 3.6324e-03;
    Pa4_c = 1.4847e+03;
    Pa5_c = 6.7505e-05;
    Pb5_c = 1.7352e-06;
    Pa6_c = 1;
    Pb6_c = 1;
    Pa7_c = 1;
    Pb7_c = 1;
    Pa3_n = 2.6452;
    Pa4_n = 4.2385e+01;
    Pb5_n = 2.1181;
    Pa6_n = 1;
    Pa7_n = 1;
elseif drug_index == 4     % GS967             %drug = 1 * (1E-6); % (M)
    % parameters -> only neutral layer, no charged layer!
    diffusion = 27314.13689124808;
    pH = 1000;%7.4; % to keep portion=0
    pKa = 0;%7.2;   % to keep portion=0
    dd = -0.7;
    kd0 = 0;%100.5 * (1e-6);
    kd0_b = 0;%1.5012 * (1e-6);
    k_off_0 = 582.095765800779 * (1e-6);		%400 * (1e-6);
    ki_off_0 = 5.000980473143621e-02 * (1e-6);	%5.4 * (1e-6);
    kc_off_0 = 884.1020871925977 * (1e-6);      %800 * (1e-6);
    Pa3_c = 1;%3.6811;
    Pa4_c = 1;%6.8705e+04;
    Pa5_c = 1;%4.0832e-02;
    Pb5_c = 1;%1.7561e-01;
    Pa6_c = 1;
    Pb6_c = 1;
    Pa7_c = 1;
    Pb7_c = 1;
    Pa3_n = 6392.988579704396;%2.3570e+02;
    Pa4_n = 3808.596738506219;%2.1182e+02;
    Pb5_n = 0.02248682933320539;%1.2197e-03;
    Pa6_n = 1;
    Pa7_n = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition rates
Q10_INa = 2.1; % Tfactor = 1 @ T = 300 K, Tfactor = Q10 @ T = 310 K
Tfactor_INa = Q10_INa^((Temp-300)/10);

% drug
portion = 1/(1+(10^(pH-pKa)));
drug_charged = drug * portion;
drug_neutral = drug * (1-portion);
kd_open = kd0 * exp(dd*y(39)*Frdy/(R*Temp));
kd_open_b = kd0_b * exp(dd*y(39)*Frdy/(R*Temp));

% charged drug
kon = drug_charged * diffusion;
koff = kd_open * diffusion;
kcon = kon;
kcoff = koff;
kbon = kon; % bursting
kboff = kd_open_b * diffusion;
kcbon = kbon;
kcboff = kboff;

% neutral drug
k_on = drug_neutral * diffusion;
k_off = k_off_0 * diffusion;
ki_on = k_on/1;
ki_off = ki_off_0 * diffusion;
kc_on = k_on/1;
kc_off = kc_off_0 * diffusion;
% kb_on = k_on; % bursting
% kb_off = k_off;
% kbc_on = kc_on;
% kbc_off = kc_off;

% Drug Free
alphaNa1 = Tfactor_INa * P1a1/(P2a1*exp(-(y(39)+P3a1)/P4a1)+P5a1*exp(-(y(39)+P3a1)/P6a1));
alphaNa2 = Tfactor_INa * P1a1/(P2a1*exp(-(y(39)+P3a1)/P4a2)+P5a2*exp(-(y(39)+P3a1)/P6a1));
alphaNa3 = Tfactor_INa * P1a1/(P2a1*exp(-(y(39)+P3a1)/P4a3)+P5a3*exp(-(y(39)+P3a1)/P6a1));
betaNa1 = Tfactor_INa * P1b1*exp(-(y(39)+P3a1)/P2b1); % shift
betaNa2 = Tfactor_INa * P1b2*exp(-(y(39)-P2b2)/P2b1);
betaNa3 = Tfactor_INa * P1b3*exp(-(y(39)-P2b3)/P2b1);
alphaNa4 = Tfactor_INa * 1/(P1a4*exp(-(y(39)+P4a4)/P2a4)+P3a4);
alphaNa5 = Tfactor_INa * P1a5*exp(-(y(39)+P4a4)/P2a5);
betaNa5 = Tfactor_INa * (P1b5+P2b5*(y(39)+P4a4));
betaNa6 = Tfactor_INa * P1b6*exp(-(y(39))/P2b6);
alphaNa7 = Tfactor_INa * P1a7*exp((y(39))/P2a7);
betaNa7 = Tfactor_INa * P1b7*exp(-(y(39))/P2b7);
alphaNa8 = Tfactor_INa * P1a8;
betaNa8 = Tfactor_INa * P1b8;
alphaNa6 = alphaNa4/P1a6;
betaNa4 = (alphaNa3*alphaNa4*alphaNa5)/(betaNa3*betaNa5); % REV

% Charged Drug
alphaNa1_c = alphaNa1; % constrained
alphaNa2_c = alphaNa2; % constrained
alphaNa3_c = Pa3_c * alphaNa3;                      % can be changed
betaNa1_c = betaNa1; % constrained
betaNa2_c = betaNa2; % constrained
%betaNa3_c = betaNa3; % constrained (REV)
alphaNa4_c = Pa4_c * alphaNa4;                      % can be changed
alphaNa5_c = Pa5_c * alphaNa5;                      % can be changed
betaNa5_c = Pb5_c * betaNa5;                        % can be changed
alphaNa6_c = Pa6_c * alphaNa6;                      % can be changed
betaNa6_c = Pb6_c * betaNa6;                        % can be changed
alphaNa7_c = Pa7_c * alphaNa7;                      % can be changed
betaNa7_c = Pb7_c * betaNa7;                        % can be changed
alphaNa8_c = alphaNa8; % constrained
betaNa8_c = betaNa8; % constrained
%betaNa4_c = betaNa4; % constrained (REV)

% Neutral Drug
alphaNa1_n = alphaNa1; % constrained
alphaNa2_n = alphaNa2; % constrained
alphaNa3_n = Pa3_n * alphaNa3;                      % can be changed
betaNa1_n = betaNa1; % constrained
betaNa2_n = betaNa2; % constrained
%betaNa3_n = betaNa3; % constrained (REV)
alphaNa4_n = Pa4_n * alphaNa4;                      % can be changed
%alphaNa5_n = alphaNa5; % constrained (REV)
betaNa5_n = Pb5_n * betaNa5;                        % can be changed
alphaNa6_n = Pa6_n * alphaNa6;                      % can be changed
alphaNa7_n = Pa7_n * alphaNa7;                      % can be changed
betaNa6_n = alphaNa6_n*betaNa6/alphaNa6; % constrained (REV)
betaNa7_n = alphaNa7_n*betaNa7/alphaNa7; % constrained (REV)
alphaNa8_n = alphaNa8; % constrained
betaNa8_n = betaNa8; % constrained
%betaNa4_n = betaNa4; % constrained (REV)

% Microscopic reversibility (REV)
if ( drug == 0 || drug_charged == 0 )
    betaNa3_c = 0;
else
    betaNa3_c = ( betaNa3 * kcon * koff * alphaNa3_c ) / ( kon * kcoff * alphaNa3);
end
if ( betaNa3_c == 0 )
    betaNa4_c = 0;
else
    betaNa4_c = ( alphaNa3_c * alphaNa4_c * alphaNa5_c ) / ( betaNa3_c * betaNa5_c );
end
if ( drug == 0 || drug_neutral == 0 )
    alphaNa5_n = 0;
else
    alphaNa5_n = ( ki_off * alphaNa5 * kc_on * betaNa5_n ) / ( ki_on * kc_off * betaNa5 );
end
if ( drug == 0 || drug_neutral == 0 )
    betaNa3_n = 0;
else
    betaNa3_n = ( betaNa3 * kc_on * alphaNa3_n * k_off ) / ( kc_off * alphaNa3 * k_on );
end
if ( betaNa3_n == 0 )
    betaNa4_n = 0;
else
    betaNa4_n = ( alphaNa5_n * alphaNa3_n * alphaNa4_n ) / ( betaNa5_n * betaNa3_n );
end

% State variables
CNa3 = y(63); CNa2 = y(64); CNa1 = y(65); ONa = y(66);
LCNa3 = y(67); LCNa2 = y(68); LCNa1 = y(69); LONa = y(70);
ICNa3 = y(71); ICNa2 = y(72); IFNa = y(73); I1Na = y(74);
%I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3));

CNa3_c = y(75); CNa2_c = y(76); CNa1_c = y(77); ONa_c = y(78);
LCNa3_c = y(79); LCNa2_c = y(80); LCNa1_c = y(81); LONa_c = y(82);
ICNa3_c = y(83); ICNa2_c = y(84); IFNa_c = y(85); I1Na_c = y(86); I2Na_c = y(87);

CNa3_n = y(88); CNa2_n = y(89); CNa1_n = y(90); ONa_n = y(91);
LCNa3_n = y(92); LCNa2_n = y(93); LCNa1_n = y(94); LONa_n = y(95);
ICNa3_n = y(96); ICNa2_n = y(97); IFNa_n = y(98); I1Na_n = y(99); I2Na_n = y(100);

I2Na = ( 1 - (ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3+...
    +ONa_c+CNa1_c+CNa2_c+CNa3_c+IFNa_c+I1Na_c+I2Na_c+ICNa2_c+ICNa3_c+LONa_c+LCNa1_c+LCNa2_c+LCNa3_c+...
    +ONa_n+CNa1_n+CNa2_n+CNa3_n+IFNa_n+I1Na_n+I2Na_n+ICNa2_n+ICNa3_n+LONa_n+LCNa1_n+LCNa2_n+LCNa3_n) );

% ODEs - Drug Free
coeff_CNa2  = (betaNa1+alphaNa2+betaNa5+alphaNa8 +kcon+kc_on);
coeff_CNa1  = (betaNa2+alphaNa3+betaNa5+alphaNa8 +kcon+kc_on);
coeff_ONa   = (betaNa3+alphaNa4+alphaNa8 +kon+k_on);
coeff_IFNa  = (betaNa4+alphaNa5+alphaNa6+betaNa2 +ki_on);
coeff_I1Na  = (betaNa6+alphaNa7 +ki_on);
coeff_CNa3  = (alphaNa1+betaNa5+alphaNa8 +kcon+kc_on);
coeff_ICNa2 = (betaNa1+alphaNa2+alphaNa5 +ki_on);
coeff_ICNa3 = (alphaNa1+alphaNa5 +ki_on);
coeff_LONa  = (betaNa8+betaNa3 +kbon+k_on);
coeff_LCNa1 = (betaNa8+betaNa2+alphaNa3 +kcbon+kc_on);
coeff_LCNa2 = (betaNa8+betaNa1+alphaNa2 +kcbon+kc_on);
coeff_LCNa3 = (betaNa8+alphaNa1 +kcbon+kc_on);
%coeff_I2Na  = (betaNa7 +ki_on);

dCNa2  = kcoff*CNa2_c+kc_off*CNa2_n+ betaNa8*LCNa2+alphaNa1*CNa3+betaNa2*CNa1+alphaNa5*ICNa2-(coeff_CNa2)*CNa2;
dCNa1  = kcoff*CNa1_c+kc_off*CNa1_n+ betaNa8*LCNa1+alphaNa2*CNa2+betaNa3*ONa+alphaNa5*IFNa-(coeff_CNa1)*CNa1;
dONa   = koff*ONa_c+k_off*ONa_n+ betaNa8*LONa+alphaNa3*CNa1+betaNa4*IFNa-(coeff_ONa)*ONa;
dIFNa  = ki_off*IFNa_n+ alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa2*ICNa2-(coeff_IFNa)*IFNa;
dI1Na  = ki_off*I1Na_n+ alphaNa6*IFNa+betaNa7*I2Na-(coeff_I1Na)*I1Na;
dCNa3  = kcoff*CNa3_c+kc_off*CNa3_n+ betaNa8*LCNa3+betaNa1*CNa2+alphaNa5*ICNa3-(coeff_CNa3)*CNa3;
dICNa2 = ki_off*ICNa2_n+ alphaNa1*ICNa3+betaNa2*IFNa+betaNa5*CNa2-(coeff_ICNa2)*ICNa2;
dICNa3 = ki_off*ICNa3_n+ betaNa1*ICNa2+betaNa5*CNa3-(coeff_ICNa3)*ICNa3;
dLONa  = kboff*LONa_c+k_off*LONa_n+ alphaNa3*LCNa1+alphaNa8*ONa-(coeff_LONa)*LONa;
dLCNa1 = kcboff*LCNa1_c+kc_off*LCNa1_n+ alphaNa8*CNa1+alphaNa2*LCNa2+betaNa3*LONa-(coeff_LCNa1)*LCNa1;
dLCNa2 = kcboff*LCNa2_c+kc_off*LCNa2_n+ betaNa2*LCNa1+alphaNa8*CNa2+alphaNa1*LCNa3-(coeff_LCNa2)*LCNa2;
dLCNa3 = kcboff*LCNa3_c+kc_off*LCNa3_n+ alphaNa8*CNa3+betaNa1*LCNa2-(coeff_LCNa3)*LCNa3;
%dI2Na  = ki_off*I2Na_n+ alphaNa7*I1Na-(coeff_I2Na)*I2Na;

% ODEs - Charged Drug
coeff_CNa2_c  = (betaNa1_c+alphaNa2_c+betaNa5_c+alphaNa8_c +kcoff);
coeff_CNa1_c  = (betaNa2_c+alphaNa3_c+betaNa5_c+alphaNa8_c +kcoff);
coeff_ONa_c   = (betaNa3_c+alphaNa4_c+alphaNa8_c +koff);
coeff_IFNa_c  = (betaNa4_c+alphaNa5_c+alphaNa6_c+betaNa2_c);
coeff_I1Na_c  = (betaNa6_c+alphaNa7_c);
coeff_CNa3_c  = (alphaNa1_c+betaNa5_c+alphaNa8_c +kcoff);
coeff_ICNa2_c = (betaNa1_c+alphaNa2_c+alphaNa5_c);
coeff_ICNa3_c = (alphaNa1_c+alphaNa5_c);
coeff_LONa_c  = (betaNa8_c+betaNa3_c +kboff);
coeff_LCNa1_c = (betaNa8_c+betaNa2_c+alphaNa3_c +kcboff);
coeff_LCNa2_c = (betaNa8_c+betaNa1_c+alphaNa2_c +kcboff);
coeff_LCNa3_c = (betaNa8_c+alphaNa1_c +kcboff);
coeff_I2Na_c  = (betaNa7_c);

dCNa2_c  = kcon*CNa2+ betaNa8_c*LCNa2_c+alphaNa1_c*CNa3_c+betaNa2_c*CNa1_c+alphaNa5_c*ICNa2_c-(coeff_CNa2_c)*CNa2_c;
dCNa1_c  = kcon*CNa1+ betaNa8_c*LCNa1_c+alphaNa2_c*CNa2_c+betaNa3_c*ONa_c+alphaNa5_c*IFNa_c-(coeff_CNa1_c)*CNa1_c;
dONa_c   = kon*ONa+ betaNa8_c*LONa_c+alphaNa3_c*CNa1_c+betaNa4_c*IFNa_c-(coeff_ONa_c)*ONa_c;
dIFNa_c  = alphaNa4_c*ONa_c+betaNa5_c*CNa1_c+betaNa6_c*I1Na_c+alphaNa2_c*ICNa2_c-(coeff_IFNa_c)*IFNa_c;
dI1Na_c  = alphaNa6_c*IFNa_c+betaNa7_c*I2Na_c-(coeff_I1Na_c)*I1Na_c;
dCNa3_c  = kcon*CNa3+ betaNa8_c*LCNa3_c+betaNa1_c*CNa2_c+alphaNa5_c*ICNa3_c-(coeff_CNa3_c)*CNa3_c;
dICNa2_c = alphaNa1_c*ICNa3_c+betaNa2_c*IFNa_c+betaNa5_c*CNa2_c-(coeff_ICNa2_c)*ICNa2_c;
dICNa3_c = betaNa1_c*ICNa2_c+betaNa5_c*CNa3_c-(coeff_ICNa3_c)*ICNa3_c;
dLONa_c  = kbon*LONa+ alphaNa3_c*LCNa1_c+alphaNa8_c*ONa_c-(coeff_LONa_c)*LONa_c;
dLCNa1_c = kcbon*LCNa1+ alphaNa8_c*CNa1_c+alphaNa2_c*LCNa2_c+betaNa3_c*LONa_c-(coeff_LCNa1_c)*LCNa1_c;
dLCNa2_c = kcbon*LCNa2+ betaNa2_c*LCNa1_c+alphaNa8_c*CNa2_c+alphaNa1_c*LCNa3_c-(coeff_LCNa2_c)*LCNa2_c;
dLCNa3_c = kcbon*LCNa3+ alphaNa8_c*CNa3_c+betaNa1_c*LCNa2_c-(coeff_LCNa3_c)*LCNa3_c;
dI2Na_c  = alphaNa7_c*I1Na_c-(coeff_I2Na_c)*I2Na_c;

% ODEs - Neutral Drug
coeff_CNa2_n  = (betaNa1_n+alphaNa2_n+betaNa5_n+alphaNa8_n +kc_off);
coeff_CNa1_n  = (betaNa2_n+alphaNa3_n+betaNa5_n+alphaNa8_n +kc_off);
coeff_ONa_n   = (betaNa3_n+alphaNa4_n+alphaNa8_n +k_off);
coeff_IFNa_n  = (betaNa4_n+alphaNa5_n+alphaNa6_n+betaNa2_n +ki_off);
coeff_I1Na_n  = (betaNa6_n+alphaNa7_n +ki_off);
coeff_CNa3_n  = (alphaNa1_n+betaNa5_n+alphaNa8_n +kc_off);
coeff_ICNa2_n = (betaNa1_n+alphaNa2_n+alphaNa5_n +ki_off);
coeff_ICNa3_n = (alphaNa1_n+alphaNa5_n +ki_off);
coeff_LONa_n  = (betaNa8_n+betaNa3_n +k_off);
coeff_LCNa1_n = (betaNa8_n+betaNa2_n+alphaNa3_n +kc_off);
coeff_LCNa2_n = (betaNa8_n+betaNa1_n+alphaNa2_n +kc_off);
coeff_LCNa3_n = (betaNa8_n+alphaNa1_n +kc_off);
coeff_I2Na_n  = (betaNa7_n +ki_off);

dCNa2_n  = kc_on*CNa2+ betaNa8_n*LCNa2_n+alphaNa1_n*CNa3_n+betaNa2_n*CNa1_n+alphaNa5_n*ICNa2_n-(coeff_CNa2_n)*CNa2_n;
dCNa1_n  = kc_on*CNa1+ betaNa8_n*LCNa1_n+alphaNa2_n*CNa2_n+betaNa3_n*ONa_n+alphaNa5_n*IFNa_n-(coeff_CNa1_n)*CNa1_n;
dONa_n   = k_on*ONa+ betaNa8_n*LONa_n+alphaNa3_n*CNa1_n+betaNa4_n*IFNa_n-(coeff_ONa_n)*ONa_n;
dIFNa_n  = ki_on*IFNa+ alphaNa4_n*ONa_n+betaNa5_n*CNa1_n+betaNa6_n*I1Na_n+alphaNa2_n*ICNa2_n-(coeff_IFNa_n)*IFNa_n;
dI1Na_n  = ki_on*I1Na+ alphaNa6_n*IFNa_n+betaNa7_n*I2Na_n-(coeff_I1Na_n)*I1Na_n;
dCNa3_n  = kc_on*CNa3+ betaNa8_n*LCNa3_n+betaNa1_n*CNa2_n+alphaNa5_n*ICNa3_n-(coeff_CNa3_n)*CNa3_n;
dICNa2_n = ki_on*ICNa2+ alphaNa1_n*ICNa3_n+betaNa2_n*IFNa_n+betaNa5_n*CNa2_n-(coeff_ICNa2_n)*ICNa2_n;
dICNa3_n = ki_on*ICNa3+ betaNa1_n*ICNa2_n+betaNa5_n*CNa3_n-(coeff_ICNa3_n)*ICNa3_n;
dLONa_n  = k_on*LONa+ alphaNa3_n*LCNa1_n+alphaNa8_n*ONa_n-(coeff_LONa_n)*LONa_n;
dLCNa1_n = kc_on*LCNa1+ alphaNa8_n*CNa1_n+alphaNa2_n*LCNa2_n+betaNa3_n*LONa_n-(coeff_LCNa1_n)*LCNa1_n;
dLCNa2_n = kc_on*LCNa2+ betaNa2_n*LCNa1_n+alphaNa8_n*CNa2_n+alphaNa1_n*LCNa3_n-(coeff_LCNa2_n)*LCNa2_n;
dLCNa3_n = kc_on*LCNa3+ alphaNa8_n*CNa3_n+betaNa1_n*LCNa2_n-(coeff_LCNa3_n)*LCNa3_n;
dI2Na_n  = ki_on*I2Na+ alphaNa7_n*I1Na_n-(coeff_I2Na_n)*I2Na_n;

% ydot(63) = dCNa3; ydot(64) = dCNa2; ydot(65) = dCNa1; ydot(66) = dONa;
% ydot(67) = dLCNa3; ydot(68) = dLCNa2; ydot(69) = dLCNa1; ydot(70) = dLONa;
% ydot(71) = dICNa3; ydot(72) = dICNa2; ydot(73) = dIFNa; ydot(74) = dI1Na; % ydot(75) = dI2Na;
%
%     ydot(75) = dCNa3_c; ydot(76) = dCNa2_c; ydot(77) = dCNa1_c; ydot(78) = dONa_c;
%     ydot(79) = dLCNa3_c; ydot(80) = dLCNa2_c; ydot(81) = dLCNa1_c; ydot(82) = dLONa_c;
%     ydot(83) = dICNa3_c; ydot(84) = dICNa2_c; ydot(85) = dIFNa_c; ydot(86) = dI1Na_c; ydot(87) = dI2Na_c;
%
%     ydot(88) = dCNa3_n; ydot(89) = dCNa2_n; ydot(90) = dCNa1_n; ydot(91) = dONa_n;
%     ydot(92) = dLCNa3_n; ydot(93) = dLCNa2_n; ydot(94) = dLCNa1_n; ydot(95) = dLONa_n;
%     ydot(96) = dICNa3_n; ydot(97) = dICNa2_n; ydot(98) = dIFNa_n; ydot(99) = dI1Na_n; ydot(100) = dI2Na_n;

ydot(63:74) = [dCNa3 dCNa2 dCNa1 dONa dLCNa3 dLCNa2 dLCNa1 dLONa dICNa3 dICNa2 dIFNa dI1Na]; % ydot(75) = dI2Na;
ydot(75:87) = [dCNa3_c dCNa2_c dCNa1_c dONa_c dLCNa3_c dLCNa2_c dLCNa1_c dLONa_c dICNa3_c dICNa2_c dIFNa_c dI1Na_c dI2Na_c];
ydot(88:100) = [dCNa3_n dCNa2_n dCNa1_n dONa_n dLCNa3_n dLCNa2_n dLCNa1_n dLONa_n dICNa3_n dICNa2_n dIFNa_n dI1Na_n dI2Na_n];

I_Na_junc2 = Fjunc*GNa*(ONa+LONa)*(y(39)-ena_junc); % +(y(39)>-82.5)*0.2e-6 %Markov
I_Na_sl2 = Fsl*GNa*(ONa+LONa)*(y(39)-ena_sl); % +(y(39)>-82.5)*0.2e-6 %Markov
% I_Na2=I_Na_junc2+I_Na_sl2;

% Compute total INa (HH or Markov)
I_Na_junc = (I_Na_junc1+I_NaL_junc)*(1-flagMina)+I_Na_junc2*flagMina; %-I_Na_junc1*(t>21010);
I_Na_sl = (I_Na_sl1+I_NaL_sl)*(1-flagMina)+I_Na_sl2*flagMina; %-I_Na_sl1*(t>21010);
I_Na = I_Na_junc+I_Na_sl;

%% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

%% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = 1*Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = 1*Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
if (flagMina == 1) && (drug_index == 1) % ranolazine
    IC50_kr = 35*(1e-6);
    factor_rano_kr = 1/(1+(drug/IC50_kr));
else
    factor_rano_kr = 1;
end
gkr = S_GKr*0.035*sqrt(Ko/5.4)*factor_rano_kr;

%gkr =0.035*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+10)/5));
tauxr = 550/(1+exp((-22-y(39))/9))*6/(1+exp((y(39)-(-11))/9))+230/(1+exp((y(39)-(-40))/20));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+74)/24));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

%% I_ks: Slowly Activating K Current
% pcaks_junc = -log10(y(36))+3.0;
% pcaks_sl = -log10(y(37))+3.0;
% gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
% gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));
eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));

markov_iks = 0;
if markov_iks == 0
    gks_junc=S_GKs*(1+1*AF+2*ISO)*0.0035*1;
    gks_sl=S_GKs*(1+1*AF+2*ISO)*0.0035*1; %FRA
    xsss = 1 / (1+exp(-(y(39)+40*ISO + 3.8)/14.25)); % fitting Fra
    tauxs=990.1/(1+exp(-(y(39)+40*ISO+2.436)/14.12));
    ydot(13) = (xsss-y(13))/tauxs;
    I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
    I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);
    I_ks = I_ks_junc+I_ks_sl;
else
    gks_junc=S_GKs*1*0.0065;
    gks_sl=S_GKs*1*0.0065; %FRA
    alpha=3.98e-4*exp(3.61e-1*y(39)*FoRT);
    beta=5.74e-5*exp(-9.23e-2*y(39)*FoRT);
    gamma=3.41e-3*exp(8.68e-1*y(39)*FoRT);
    delta=1.2e-3*exp(-3.3e-1*y(39)*FoRT);
    teta=6.47e-3;
    eta=1.25e-2*exp(-4.81e-1*y(39)*FoRT);
    psi=6.33e-3*exp(1.27*y(39)*FoRT);
    omega=4.91e-3*exp(-6.79e-1*y(39)*FoRT);

    ydot(42)=-4*alpha*y(42)+beta*y(43);
    ydot(43)=4*alpha*y(42)-(beta+gamma+3*alpha)*y(43)+2*beta*y(44);
    ydot(44)=3*alpha*y(43)-(2*beta+2*gamma+2*alpha)*y(44)+3*beta*y(45);
    ydot(45)=2*alpha*y(44)-(3*beta+3*gamma+alpha)*y(45)+4*beta*y(46);
    ydot(46)=1*alpha*y(44)-(4*beta+4*gamma)*y(46)+delta*y(50);
    ydot(47)=gamma*y(43)-(delta+3*alpha)*y(47)+beta*y(48);
    ydot(48)=2*gamma*y(44)+3*alpha*y(47)-(delta+beta+2*alpha+gamma)*y(48)+2*beta*y(49)+2*delta*y(51);
    ydot(49)=3*gamma*y(45)+2*alpha*y(48)-(delta+2*beta+1*alpha+2*gamma)*y(49)+3*beta*y(50)+2*delta*y(52);
    ydot(50)=4*gamma*y(46)+1*alpha*y(49)-(delta+3*beta+0*alpha+3*gamma)*y(50)+2*delta*y(53);
    ydot(51)=1*gamma*y(48)-(2*delta+2*alpha)*y(51)+beta*y(52);
    ydot(52)=2*gamma*y(49)+2*alpha*y(51)-(2*delta+beta+1*alpha+gamma)*y(52)+2*beta*y(53)+3*delta*y(54);
    ydot(53)=3*gamma*y(50)+1*alpha*y(52)-(2*delta+2*beta+2*gamma)*y(53)+3*delta*y(55);
    ydot(54)=1*gamma*y(52)-(3*delta+1*alpha)*y(54)+beta*y(55);
    ydot(55)=2*gamma*y(53)+1*alpha*y(54)-(3*delta+1*beta+1*gamma)*y(55)+4*delta*y(56);
    ydot(56)=1*gamma*y(55)-(4*delta+teta)*y(56)+eta*y(57);
    O2=1-(y(42)+y(43)+y(44)+y(45)+y(46)+y(47)+y(49)+y(48)+y(50)+y(51)+y(52)+y(53)+y(54)+y(55)+y(56)+y(57));
    ydot(57)=1*teta*y(56)-(eta+psi)*y(57)+omega*O2;
    I_ks_junc = Fjunc*gks_junc*(y(57)+O2)*(y(39)-eks);
    I_ks_sl = Fsl*gks_sl*(y(57)+O2)*(y(39)-eks);
    I_ks = I_ks_junc+I_ks_sl;
end

%% I_kp: Plateau K Current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_k,ach: Muscarinic-Receptor-Activated K Current
%Ach=1*0.1; % [Ach in uM]
I_KAch =S_GKAch* 1/(1+(0.03/Ach)^2.1)*(0.08+0.4./(1+exp((y(39)+91)/12))).*(y(39)-ek);

%% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes
GtoFast=S_Gto*(1.0-0.7*AF)*0.165*1.0; %nS/pF maleckar; %human atrium

% equations for activation;
xtoss = ( (1)./ ( 1 + exp( -(y(39)+1.0)/11.0 ) ) );
tauxtof = 3.5*exp(-((y(39)/30.0)^2.0))+1.5;

% equations for inactivation;
ytoss = ( (1.0)./ ( 1 + exp( (y(39)+40.5)/11.5) ) ) ;
tauytof = 25.635*exp(-(((y(39)+52.45)/15.8827)^2.0))+24.14;%14.14

ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = 1.0*GtoFast*y(10)*y(11)*(y(39)-ek);

I_to = I_tof;

%% I_kur: Ultra rapid delayed rectifier Outward K Current
% Equation for IKur; from Maleckar et al. 2009 - EG
% atrium

%nS/pF maleckar 0.045
gkur_hh = 0.85*S_GKur*(1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA);   %after incorporation of ISK and IK2P, IKur reduced by 15% first.
% Chaos/Frontiers paper version
% gkur_hh = (1-Complete_Block)*(1.0-0.5*IKur_ur_50)*(1.0-0.5*AF)*(1+2*ISO)*0.045*1.36*(1+0.2*RA); %nS/pF maleckar 0.045*1.36 based on Feng physiological data (Figure 4)
% testing 1
% gkur_hh = (1-Complete_Block)*(1.0-0.5*IKur_ur_50)*(1.0-0.5*AF)*(1+2*ISO)*0.045*1.45*(1+0.2*RA); %nS/pF maleckar 0.045*1.36 based on Feng physiological data (Figure 4)
% testing 2, scaled conductance
% gkur_hh = (1-Complete_Block)*(1.0-0.5*IKur_ur_50)*(1.0-0.5*AF)*(1+2*ISO)*0.045*1.06*(1+0.2*RA); %nS/pF maleckar 0.045*1.36 based on Feng physiological data (Figure 4)

% equations for activation;
xkurss = ( (1)./ ( 1 + exp( (y(39)+6)/-8.6 ) ) );
tauxkur = 9/(1+exp((y(39)+5)/12.0))+0.5;
% equations for inactivation;
ykurss = ( (1)./ ( 1 + exp( (y(39)+7.5)/10 ) ) );
tauykur = 590/(1+exp((y(39)+60)/10.0))+3050;

ydot(58) = (xkurss-y(58))/tauxkur;
ydot(59) = (ykurss-y(59))/tauykur;
I_kur_hh = gkur_hh*y(58)*y(59)*(y(39)-ek);

%% I_kur: Markov model (w/ drugs)
y_kur = y(100+1:100+11);
ydot_kur = ydot(100+1:100+11);

% Conductance
%gkur_m = 0.100; %0.0975*(0.35)*(1.05); % modify also gtof = 0.165*(0.65);
% gkur_m = (1-Complete_Block)*(1.0-0.5*IKur_ur_50)*(1.0-0.5*AF)*(1+2*ISO)*0.045*1.06*(1+0.2*RA); %nS/pF maleckar 0.045*1.36 based on Feng physiological data (Figure 4)
gkur_m = S_GKur*(1-Complete_Block)*(1.0-0.5*IKur_ur_50)*(1.0-0.5*AF)*(1+2*ISO)*0.045*1.36*(1+0.2*RA); %nS/pF maleckar 0.045*1.36 based on Feng physiological data (Figure 4)
% gkur_m = 0.9*gkur_hh;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhou et al base model %
%%%%%%%%%%%%%%%%%%%%%%%%%
%shift = 0; ktau = 1; kinact = 1;

% shift = -15; ktau = 3; kinact = 1/40 + (39/40./(1 + exp((y(39)+5)./10)));
shift = -27.5; kinact = 1./(1 + exp((y(39)+10)./10)); ktau =  5; kall = 2.65;

% TRANSITION RATES (a1_ur, b1_ur, kfur, kbur)
% C4ur to O in ms^-1
% a1_ur = exp((y(39) - 48.4 + shift)/49.2); % modified model original
a1_ur = kall*(exp((y(39) - 48.4 + shift)./70)); %modified model new
% a1_ur = 8*(exp((y(39) - 48.4 + shift)./27)); %modified model new, EADS
% a1_ur = 9*(exp((y(39) - 76)./25.4)); %modified model new, v2
% a1_ur = 3*(exp((y(39) - 48.4 + shift)./70)); %modified model new
% C2ur to C1ur in ms^-1
b1_ur = kall*(exp((y(39) - 48.4 + shift)./70) .* (exp(-(y(39)+ 48.4 + shift)./10) ./ ...
    (1 + 0.3.*exp(-(y(39) + 48.4 + shift)./10)))); %modified new
% b1_ur = 0.04*(exp((y(39) - 48.4 + shift)./27) .* (exp(-4*(y(39))./27) ./ ...
%     (1 + 0.023.*exp(-1.4*(y(39) + 48.4 + shift)./10)))); %EADS
% b1_ur = 0.04*(exp((y(39) - 76)./27) .* (exp(-4*(y(39))./24) ./ ...
%     (1 + 0.023.*exp(-1.7*(y(39) + 20.9)./10)))); %modified model changed beta, v2
% b1_ur = kall*b1_ur;
% b1_ur = 0.75*(exp((y(39) - 48.4 + shift)./70) .* (exp(-(y(39) + 48.4 + shift)./10) ./ ...
%     (1 + 0.3.*exp(-(y(39) + 48.4 + shift)./10))));
% O to I in ms^-1
kf_ur = ktau*0.0001;
% I to O in ms^-1
kb_ur = ktau*0.0001*kinact;
% C1ur to C2ur in ms^-1
a4_ur = 4*a1_ur;
% C2ur to C3ur in ms^-1
a3_ur = 3*a1_ur;
% C3 to C2 in ms^-1
b2_ur = 2*b1_ur;
% C3ur to C4ur in ms^-1
a2_ur = 2*a1_ur;
% C4ur to C3ur in ms^-1
b3_ur = 3*b1_ur;
% O to C4ur in ms^-1
b4_ur = 4*b1_ur;
% Drug-free to Drug-bound in ms^-1
db_ur = drug_kur_conc*kon_kur;
% Drug-bound to Drug-free in ms^-1
dr_ur = koff_kur;

I_ur = 1 - sum(y_kur(1:11));

if drug_kur_index == 0 % DRUG-FREE
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur - y_kur(5)*(b4_ur + kf_ur);
    % ydot_kur(6-11) still set to 0

elseif drug_kur_index == 1 % CLOSED STATE
    ydot_kur(1) = b1_ur*y_kur(2) + dr_ur*y_kur(6) - y_kur(1)*(a4_ur + db_ur);
    ydot_kur(2) = a4_ur*y_kur(1) + dr_ur*y_kur(7) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur + db_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + dr_ur*y_kur(8) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur + db_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + dr_ur*y_kur(9) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur + db_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur - y_kur(5)*(b4_ur + kf_ur);
    % I_ur
    ydot_kur(6) = b1_ur*y_kur(7) + db_ur*y_kur(1) - y_kur(6)*(a4_ur + dr_ur);
    ydot_kur(7) = a4_ur*y_kur(6) + db_ur*y_kur(2) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur + dr_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + db_ur*y_kur(3) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur + dr_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + db_ur*y_kur(4) + b4_ur*y_kur(10) - y_kur(9)*(b3_ur + a1_ur + dr_ur);
    ydot_kur(10) = a1_ur*y_kur(9) + kb_ur*y_kur(11) - y_kur(10)*(b4_ur + kf_ur);
    ydot_kur(11) =  kf_ur*y_kur(10) - kb_ur*y_kur(11);

elseif drug_kur_index == 2 % OPEN STATE
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(10) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    % I_ur
    ydot_kur(6) = b1_ur*y_kur(7) - a4_ur*y_kur(6);
    ydot_kur(7) = a4_ur*y_kur(6) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + b4_ur*y_kur(10) - y_kur(9)*(b3_ur + a1_ur);
    ydot_kur(10) = a1_ur*y_kur(9) + kb_ur*y_kur(11) + db_ur*y_kur(5) - y_kur(10)*(b4_ur + kf_ur + dr_ur);
    ydot_kur(11) =  kf_ur*y_kur(10) - kb_ur*y_kur(11);

elseif drug_kur_index == 3 % INACTIVATED STATE
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur - y_kur(5)*(b4_ur + kf_ur);
    % I_ur
    ydot_kur(6) = b1_ur*y_kur(7) - a4_ur*y_kur(6);
    ydot_kur(7) = a4_ur*y_kur(6) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + b4_ur*y_kur(10) - y_kur(9)*(b3_ur + a1_ur);
    ydot_kur(10) = a1_ur*y_kur(9) + kb_ur*y_kur(11) - y_kur(10)*(b4_ur + kf_ur);
    ydot_kur(11) = kf_ur*y_kur(10) + db_ur*I_ur - y_kur(11)*(kb_ur + dr_ur);

elseif drug_kur_index == 4 % OPEN and INACTIVATED STATE
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(10) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    %I_ur
    ydot_kur(6) = b1_ur*y_kur(7) - a4_ur*y_kur(6);
    ydot_kur(7) = a4_ur*y_kur(6) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + b4_ur*y_kur(10) - y_kur(9)*(b3_ur + a1_ur);
    ydot_kur(10) = a1_ur*y_kur(9) + kb_ur*y_kur(11) + db_ur*y_kur(5) - y_kur(10)*(b4_ur + kf_ur + dr_ur);
    ydot_kur(11) = kf_ur*y_kur(10) + db_ur*I_ur - y_kur(11)*(kb_ur + dr_ur);

elseif drug_kur_index ==5 % OPEN and INACTIVATED and CLOSED
    ydot_kur(1) = b1_ur*y_kur(2) + dr_ur*y_kur(6) - y_kur(1)*(a4_ur + db_ur);
    ydot_kur(2) = a4_ur*y_kur(1) + dr_ur*y_kur(7) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur + db_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + dr_ur*y_kur(8) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur + db_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + dr_ur*y_kur(9) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur + db_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(10) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    % I_ur
    ydot_kur(6) = b1_ur*y_kur(7) + db_ur*y_kur(1) - y_kur(6)*(a4_ur + dr_ur);
    ydot_kur(7) = a4_ur*y_kur(6) + db_ur*y_kur(2) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur + dr_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + db_ur*y_kur(3) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur + dr_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + db_ur*y_kur(4) + b4_ur*y_kur(10) - y_kur(9)*(b3_ur + a1_ur + dr_ur);
    ydot_kur(10) = a1_ur*y_kur(9) + kb_ur*y_kur(11) + db_ur*y_kur(5) - y_kur(10)*(b4_ur + kf_ur + dr_ur);
    ydot_kur(11) = kf_ur*y_kur(10) + db_ur*I_ur - y_kur(11)*(kb_ur + dr_ur);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non-Reciprocal Models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif drug_kur_index == 6 % OPEN ONLY
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(6) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    % I_ur
    ydot_kur(6) = y_kur(5)*db_ur - y_kur(6)*dr_ur; %drug bound OPEN state

elseif drug_kur_index == 7 % INACTIVE ONLY
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur - y_kur(5)*(b4_ur + kf_ur);
    % I_ur
    ydot_kur(6) = db_ur*I_ur - y_kur(6)*dr_ur; %drug bound to inactive state

elseif drug_kur_index == 8 % ClOSED ONLY
    ydot_kur(1) = b1_ur*y_kur(2) + dr_ur*y_kur(6) - y_kur(1)*(a4_ur + db_ur);
    ydot_kur(2) = a4_ur*y_kur(1) + dr_ur*y_kur(7) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur + db_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + dr_ur*y_kur(8) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur + db_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + dr_ur*y_kur(9) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur + db_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur - y_kur(5)*(b4_ur + kf_ur);
    % I_ur
    ydot_kur(6) = b1_ur*y_kur(7) + db_ur*y_kur(1) - y_kur(6)*(a4_ur + dr_ur);
    ydot_kur(7) = a4_ur*y_kur(6) + db_ur*y_kur(2) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur + dr_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + db_ur*y_kur(3) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur + dr_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + db_ur*y_kur(4) - y_kur(9)*(b3_ur + dr_ur);

elseif drug_kur_index == 9 %Open and Closed (like 4-AP)
    ydot_kur(1) = b1_ur*y_kur(2) + dr_ur*y_kur(6) - y_kur(1)*(a4_ur + db_ur);
    ydot_kur(2) = a4_ur*y_kur(1) + dr_ur*y_kur(7) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur + db_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + dr_ur*y_kur(8) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur + db_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + dr_ur*y_kur(9) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur + db_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(10) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    %I_ur
    ydot_kur(6) = b1_ur*y_kur(7) + db_ur*y_kur(1) - y_kur(6)*(a4_ur + dr_ur);
    ydot_kur(7) = a4_ur*y_kur(6) + db_ur*y_kur(2) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur + dr_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + db_ur*y_kur(3) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur + dr_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + db_ur*y_kur(4) - y_kur(9)*(b3_ur + dr_ur);
    ydot_kur(10) = y_kur(5)*db_ur - y_kur(10)*dr_ur; %drug bound OPEN state

elseif drug_kur_index == 10 %Open and Inactive (like Lee/Hill)
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(6) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    % I_ur
    ydot_kur(6) = y_kur(5)*db_ur - y_kur(6)*dr_ur; %drug bound OPEN state
    ydot_kur(7) = db_ur*I_ur - y_kur(7)*dr_ur; %drug bound to INACTIVE state

elseif drug_kur_index == 11 %Open, Inactive, Closed (all non-reciprocal)
    ydot_kur(1) = b1_ur*y_kur(2) + dr_ur*y_kur(6) - y_kur(1)*(a4_ur + db_ur);
    ydot_kur(2) = a4_ur*y_kur(1) + dr_ur*y_kur(7) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur + db_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + dr_ur*y_kur(8) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur + db_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + dr_ur*y_kur(9) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur + db_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(10) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    %I_ur
    ydot_kur(6) = b1_ur*y_kur(7) + db_ur*y_kur(1) - y_kur(6)*(a4_ur + dr_ur);
    ydot_kur(7) = a4_ur*y_kur(6) + db_ur*y_kur(2) + b2_ur*y_kur(8) - y_kur(7)*(b1_ur + a3_ur + dr_ur);
    ydot_kur(8) = a3_ur*y_kur(7) + db_ur*y_kur(3) + b3_ur*y_kur(9) - y_kur(8)*(b2_ur + a2_ur + dr_ur);
    ydot_kur(9) = a2_ur*y_kur(8) + db_ur*y_kur(4) - y_kur(9)*(b3_ur + dr_ur);
    ydot_kur(10) = db_ur*y_kur(5) - y_kur(10)*dr_ur;
    ydot_kur(11) = db_ur*I_ur - y_kur(11)*dr_ur;


elseif drug_kur_index == 12 %Open and Inactive (like Lee/Hill), Variable Affinity
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + p(15)*y_kur(6) - y_kur(5)*(b4_ur + kf_ur + p(14)*drug_kur_conc);
    % I_ur
    ydot_kur(6) = y_kur(5)*p(14)*drug_kur_conc - y_kur(6)*p(15); %drug bound OPEN state
    ydot_kur(7) = p(18)*drug_kur_conc*I_ur - y_kur(7)*p(19); %drug bound to INACTIVE state

elseif drug_kur_index == 13 % (dO --> dI horizontal transition, for reviewer)
    ydot_kur(1) = b1_ur*y_kur(2) - a4_ur*y_kur(1);
    ydot_kur(2) = a4_ur*y_kur(1) + b2_ur*y_kur(3) - y_kur(2)*(b1_ur + a3_ur);
    ydot_kur(3) = a3_ur*y_kur(2) + b3_ur*y_kur(4) - y_kur(3)*(b2_ur + a2_ur);
    ydot_kur(4) = a2_ur*y_kur(3) + b4_ur*y_kur(5) - y_kur(4)*(b3_ur + a1_ur);
    ydot_kur(5) = a1_ur*y_kur(4) + kb_ur*I_ur + dr_ur*y_kur(6) - y_kur(5)*(b4_ur + kf_ur + db_ur);
    % I_ur
    ydot_kur(6) = y_kur(5)*db_ur + y_kur(7)*kb_ur - y_kur(6)*(dr_ur + kf_ur); %drug bound OPEN state
    ydot_kur(7) = db_ur*I_ur + y_kur(6)*kf_ur - y_kur(7)*(dr_ur + kb_ur); %drug bound to INACTIVE state

end
ydot(100+1:100+11) = ydot_kur;

I_kur_m = 0.85*gkur_m*y_kur(5)*(y(39)-ek);   %after incorporation of ISK % and IK2P, IKur reduced by 15% first.

I_kur = I_kur_m*IKur_M_flag + I_kur_hh*(1-IKur_M_flag);

%% I_k1: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);

I_ki = S_GK1 * (1+1*AF)*0.0525*sqrt(Ko/5.4)*kiss*(y(39)-ek);  % after incorporating ISK, IK1 was reduced to keep RMP.

%% I_sk: Small-conductance Ca-activated K Current
% Optimization
par_sk = [0.0506381114404388,0.273335569451572,2.96381060498817,0.199981221802789,0.279328126521496,-86.9289059836381,0.00636311816933264,5.22915055145375];

% Maximal conductance
gsk = S_GSK*par_sk(1);

% Ca-dependency
SK_shift = 0;
kdsk = (10^(SK_shift-3.45)); % (mM)
% Simulations with reduced SK-Kd (Supplemental Figure S3) were done with Kdsk = 230 nM

gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(y(36)))/0.3));
gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(y(37)))/0.3));
% Simulations with Ca-dependent block of ISK (Supplemental Figures S5 & S7) were performed using the following lines:
% kdsk_inact = 0.0193; 
% slope_inact = 0.4;
% gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(y(36)))/0.3)).*(1-1./(1+ exp((log10(kdsk_inact)-log10(y(36)))./slope_inact)));
% gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(y(37)))/0.3)).*(1-1./(1+ exp((log10(kdsk_inact)-log10(y(37)))./slope_inact)));

Fjunc_SK = Fjunc;    Fsl_SK = 1-Fjunc_SK; % Homogeneous distribution

% Simulations with increased SK channel expression in the cleft compartment
% (Supplemental Figures S6 & S7) were performed with Fjunc_SK = 0.50 (50/50 distribution)

gsk_vm = par_sk(2)/(1+exp((y(39)-ek+par_sk(3))*par_sk(4))) + par_sk(5)/(1+exp((-(y(39)-ek+par_sk(6))*par_sk(7))));

I_sk_junc = Fjunc_SK*SK_flag*gsk*gsk_ca_junc*gsk_vm*(y(39)-ek);
I_sk_sl = Fsl_SK*SK_flag*gsk*gsk_ca_sl*gsk_vm*(y(39)-ek);

I_sk = I_sk_junc + I_sk_sl;

%% I_K2P: Two-pore-domain K Current
GK2p = S_GK2P * 0.005*(1+AF*2);
o_k2p = y(111+1);

IK2p = GK2p * o_k2p * (y(39) - ek);
IK2p_Shift = -15;

inf_O_K2p = 0.2 + 0.8 / ( 1 + exp(-(y(39) - (10+AF*IK2p_Shift)) / 14.0));
tau = 2.0 + 40 / (1.0 + exp((y(39) + 25) * (y(39) + 25) / 80.0));
% o_k2p = inf + (o_k2p - inf) * exp(-(dt) / tau);

ydot(111+1) = (inf_O_K2p - y(111+1)) / tau;

%% I_ClCa & I_Clbk: Ca-activated Cl Current and Background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;

I_Clbk = GClB*(y(39)-ecl);

I_ClCFTR = GClCFTR*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+3*ISO+9)/6)); %in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
taud = 1*dss*(1-exp(-(y(39)+3*ISO+9)/6))/(0.035*(y(39)+3*ISO+9));
fss = 1/(1+exp((y(39)+3*ISO+30)/7))+0.2/(1+exp((50-y(39)-3*ISO)/20)); % in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+3*ISO+25))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-1*11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-1*11.9e-3*y(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

%% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^2);
Ka_sl = 1/(1+(Kdact/y(37))^2);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);

I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

%% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

%% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = (1)*koCa/kCaSR;%
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R; closed
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]

J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);

J_SRleak = S_RyR_leak*(1.0+0.25*AF)*5.348e-6*(y(31)-y(36));           %   [mM/ms]

%% Na and Ca Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
    +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec]
% ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur+I_KAch+I_sk + IK2p;     % [uA/uF] %SVP: added IKur
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) = 0; % -I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl

if t>=2e3 && cajl_clamp==1, ydot(36)=0; ydot(37)=0; end

% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38));
% ydot(38)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation type
switch lower(protocol)

    case 'pace_hr_control' %make HR irregular
        ind = find(stim_t>=t);
        I_app = 12.5*stim(ind(1));

    case 'pace_cc' % pace w/ current injection at rate 'rate'
        rate = pacing_rate*1e-3;
        if mod(t,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end
    case 'pace_cc_ead' % pace w/ current injection at rate 'rate'
        fast_pacing_duration = 20e3;
        if t < fast_pacing_duration
            rate = 10*1e-3;
            t0 = 0;
        else
            rate = 1*1e-3;
            t0 = fast_pacing_duration + 2000;
        end
%         if mod(t-t0,1/rate) <= 5,
        if mod(t-t0,1/rate) <= 5 && t > t0
            I_app = 12.5;
        else
            I_app = 0.0;
        end

    case 'pace_cc_ead_2beat' % pace w/ current injection at rate 'rate'
        t2 = 30; % (ms) interval before 2 beat
        rate = 1*1e-3;
        if mod(t,1/rate) <= 5 || mod(t-t2,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end

    case 'pace_cc_ead_sa' % pace w/ current injection at rate 'rate'
        %first_interval_duration = 30e3;
        fast_pacing_duration = 20e3;
        if t < first_interval_duration
            rate = 1*1e-3;
            t0 = 0;
        end
        if t >= first_interval_duration && t < first_interval_duration+fast_pacing_duration
            rate = 10*1e-3;
            t0 = first_interval_duration;
        end
        if t >= first_interval_duration+fast_pacing_duration
            rate = 1*1e-3;
            t0 = first_interval_duration+fast_pacing_duration;
        end
        if mod(t-t0,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end

    case 'pace_cc_dad' % pace w/ current injection at rate 'rate'
        if t < 10e3
            rate = pacing_rate*1e-3;
        else
            rate = 0;
        end
        if mod(t,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end


    case 'pace_cc_erp' % pace ERP protocol
        stim_time = 200;
        if t > 0 && t < stim_time
            I_app = 0;
        elseif t <= stim_time+5
            I_app = 2*DTE; % Need to calculate DTE based on conditions
        elseif t > 5 && t <= p(3)
            I_app = 0.0;
        elseif t > p(3) && t <= p(3)+5
            I_app = DTE*2; % 2*DTE
        else
            I_app = 0.0;
        end

    case 'dte_finder'  % DTE
        % model_flag = pacing_frequncies;
        if t > 1000 && t < 1005
            I_app = pacing_frequncies;
        else
            I_app = 0.0;
        end
        
    case 'tri_ap_clamp_wc'
        restEm_shift = 0*10; % (mV) shift in resting Em, same peak and slope
        rate = 0.1*1e-3;
        T_pre = 10000;
        V_rest= -85; V_peak= 35; T_tri= rec_interval; % slope = (V_rest-V_peak)/T_tri
        V_hold1 = V_rest; T_hold1 = 5;
        T_hold2 = T_tri;
        V_hold3 = V_rest;
        if t <= T_pre
            V_clamp = V_hold1;
        else
            if mod(t-T_pre,1/rate) <= T_hold1
                V_clamp = V_hold1;
            elseif mod(t-T_pre,1/rate) > T_hold1 && mod(t-T_pre,1/rate) <= T_hold1+T_hold2
                V_clamp = V_peak-(V_peak-V_rest+0*10)/T_tri*(mod(t-T_pre,1/rate)-T_hold1);
            else
                V_clamp = V_hold3;
            end
        end
        if V_clamp < V_rest+restEm_shift % add shift in resting Em
            V_clamp = V_rest+restEm_shift;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'ap_clamp' % timeVC p.emVC_1 p.emVC_2
        restEm_shift = 0*10; % (mV) shift in resting Em, same peak and slope
        rate = 1e-3;
        emVC = p.emVC_1; % original AP
        %emVC = p.emVC_2; % modified AP (10mV-shift in resting Em, same peak)
        V_rest = emVC(1);
        V_hold1 = V_rest+restEm_shift; T_hold1 = 10;
        T_hold2 = 250;%timeVC(end);
        V_hold3 = V_rest+restEm_shift;
        if mod(t,1/rate) <= T_hold1
            V_clamp = V_hold1;
        elseif mod(t,1/rate) > T_hold1 && mod(t,1/rate) <= T_hold1+T_hold2
            tt=mod(t-T_hold1,1/rate);
            t_in_roi=find(timeVC>=tt); t_in_index=t_in_roi(1);
            V_clamp = emVC(t_in_index);
            if V_clamp < V_rest+restEm_shift % add shift in resting Em
                V_clamp = V_rest+restEm_shift;
            end
        else
            V_clamp = V_hold3;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'ap_clamp_w_step' % timeVC p.emVC_1 p.emVC_2


        AP_clamp_duration = rec_interval;
        emVC = p.EAD_Vm; % original AP
        V_rest = emVC(1);
        V_hold2 = -40; T_hold2 = 2;
        if rec_interval == 0, T_hold2 = 0; end
        V_hold3 = -10; T_hold3 = 50;
        if t <= AP_clamp_duration
            t_in_roi=find(p.EAD_t>=t); t_in_index=t_in_roi(1);
            V_clamp = emVC(t_in_index);
        elseif t > AP_clamp_duration && t <= AP_clamp_duration+T_hold2
            V_clamp = V_hold2;
        elseif t > AP_clamp_duration+T_hold2 && t <= AP_clamp_duration+T_hold2+T_hold3
            V_clamp = V_hold3;
        else
            V_clamp = V_rest;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'v_hold' % Em-clamp at -140 mV
        V_clamp = -140;
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'ssa_prot' % duration 250 ms % f(vm_test)
        V_hold1 = -140; T_hold1 = 5;
        %V_hold1 = -80; T_hold1 = 5;
        V_hold2 = vm_test; T_hold2 = 200;
        V_hold3 = V_hold1; % T_hold3 = 40;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold3;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'ssi_prot' % duration 600 ms % f(vm_test)
        V_hold1 = -140; T_hold1 = 5;
        V_hold2 = vm_test; T_hold2 = 500;
        V_hold3 = -10; T_hold3 = 25;
        V_hold4 = V_hold1;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        elseif t > T_hold1+T_hold2 && t <= T_hold1+T_hold2+T_hold3
            V_clamp = V_hold3;
        else
            V_clamp = V_hold4;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'rec_prot' % rec from inactivation (P1/P2)
        V_test= -10;
        V_hold1 = -140; T_hold1 = 5;
        V_hold2 = V_test; T_hold2 = 1000;
        V_hold3 = V_hold1; T_hold3 = rec_interval;
        V_hold4 = V_test; T_hold4 = 200;
        V_hold5 = V_hold1;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        elseif t > T_hold1+T_hold2 && t <= T_hold1+T_hold2+T_hold3
            V_clamp = V_hold3;
        elseif t > T_hold1+T_hold2+T_hold3 && t <= T_hold1+T_hold2+T_hold3+T_hold4
            V_clamp = V_hold4;
        else
            V_clamp = V_hold5;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'rec_im_prot' % rec from IM inactivation (P1/P2)
        V_test= -10;
        V_hold1 = -140; T_hold1 = 5;
        V_hold2 = V_test; T_hold2 = rec_interval;
        V_hold3 = V_hold1; T_hold3 = 20;
        V_hold4 = V_test; T_hold4 = 200;
        V_hold5 = V_hold1;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        elseif t > T_hold1+T_hold2 && t <= T_hold1+T_hold2+T_hold3
            V_clamp = V_hold3;
        elseif t > T_hold1+T_hold2+T_hold3 && t <= T_hold1+T_hold2+T_hold3+T_hold4
            V_clamp = V_hold4;
        else
            V_clamp = V_hold5;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'v_step_tb' % tonic block (TB) after 30 s at -100 mV
        V_test= -10;
        %rate = 1/60*1e-3; V_hold1 = -100; T_hold1 = 3*10000+5; % w/out pulse at t=0
        rate = 1/30*1e-3; V_hold1 = -100; T_hold1 = 5; % with pulse at t=0
        V_hold2 = V_test; T_hold2 = 200;
        if mod(t,1/rate) <= T_hold1
            V_clamp = V_hold1;
        elseif mod(t,1/rate) > T_hold1 && mod(t,1/rate) <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold1;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'v_step' % use dependent block (UDB)
        rate = pacing_rate * 1e-3; 
        V_test= -10;
        V_hold1 = -100; T_hold1 = 5;
        V_hold2 = V_test; T_hold2 = 25;
            %V_hold1 = -140; % initial conditions
	    if mod(t,1/rate) <= T_hold1
            V_clamp = V_hold1;
        elseif mod(t,1/rate) > T_hold1 && mod(t,1/rate) <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold1;
        end
		R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;

    case 'v_step_interval' % recovery from UDB - IM protocol
        V_test= -10;
        V_hold1 = -100; T_hold1 = 5;
        V_hold2 = V_test; T_hold2 = 25;
        V_hold3 = V_hold1; T_hold3 = rec_interval;
        V_hold4 = V_test; T_hold4 = T_hold2;
        V_hold5 = V_hold1;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        elseif t > T_hold1+T_hold2 && t <= T_hold1+T_hold2+T_hold3
            V_clamp = V_hold3;
        elseif t > T_hold1+T_hold2+T_hold3 && t <= T_hold1+T_hold2+T_hold3+T_hold4
            V_clamp = V_hold4;
        else
            V_clamp = V_hold5;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

    case 'ap_clamp_ramp' % positive ramp
        %vm_test
        %rec_interval
        V_hold1 = vm_test; T_hold1 = 500;
        V_hold2 = vm_test+150/2000*(t-T_hold1); T_hold2 = 2000;
        V_hold3 = vm_test;
        if t <= T_hold1
            V_clamp = V_hold1;
        elseif t > T_hold1 && t <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold3;
        end
        R_clamp = 0.02;
        I_app = (V_clamp-y(39))/R_clamp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk++I_ClCFTR;              % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%ydot(39) = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);

%% Adjust output depending on the function call
if (nargin == 3)
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) && strcmp(runType,'currents')
    currents = [I_Na I_KAch I_Catot I_ncx vmax I_sk I_ki I_kur I_kur_hh I_kur_m, I_sk, IK2p, I_nak I_sk_junc I_sk_sl];
    output = currents;
end

% End of the function