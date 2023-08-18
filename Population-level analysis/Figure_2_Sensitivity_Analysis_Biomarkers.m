%% Figure 2 Results
close all
clear
clc

addpath Conditions

%% Load Data
load Parameter_matrix.mat % sigma 0.1
load biomarkers_population_1000_AF.mat
all_values_AF = all_values;
load biomarkers_population_1000_nSR.mat
%ERP Data
load ERP_1SK_populaton_1ms_AF.mat
all_values_AF_ERP = cell2mat(ERP_total);
load ERP_1SK_populaton_1ms_nSR.mat
all_values_ERP = cell2mat(ERP_total);

%% Inital variables APD90
pacing_threshold_APD90 = (all_values(:,5));
pacing_threshold_AF_APD90 = (all_values_AF(:,5));

%% Inital variables CaT Amp
pacing_threshold_CaTamp = (all_values(:,11));
pacing_threshold_CaTamp_AF = (all_values_AF(:,11));

%% Inital variables ERP
pacing_threshold_ERP = (all_values_ERP(:,1));
pacing_threshold_AF_ERP = (all_values_AF_ERP(:,1));

%% Linear Regression
parameter_names = [{'GNa'},{'GClCa'},{'GCaL'},{'Gtof'},{'GKur'},{'GKr'},{'GKs'},{'GK1'},{'VNCX'},{'VNKA'}, {'GK2P'},{'GNaB'},{'GCaB'}, ...
{'GClB'},{'GKAch'},{'GNaL'},{'GKp'},{'GSK'},{'VPMCA'},{'VSERCA'},{'VRyR'},{'VRyRLeak'}];

%% Plot 1
rm_IKACH_parameters = all_parameters;  
rm_IKACH_parameters_AF = all_parameters; 

X = log(rm_IKACH_parameters);
Y = log(pacing_threshold_APD90'); 
X_AF = log(rm_IKACH_parameters_AF);
Y_AF = log(pacing_threshold_AF_APD90'); 

sel_index = ~ isnan(Y);
sel_index_AF = ~ isnan(Y_AF);

XX = zscore(X(sel_index,:));
YY = zscore(Y(sel_index));

XX_AF = zscore(X_AF(sel_index_AF,:));
YY_AF = zscore(Y_AF(sel_index_AF));

[nn, mm ] = size(XX);
[nn_AF, mm_AF ] = size(XX_AF);

New_X = [ones(nn, 1), XX];
New_X_AF = [ones(nn_AF, 1), XX_AF];

[b,bint,r,rint,stats] = regress(YY',New_X);
stats;
b_nSR = b;
bint_nSR = bint;
r_nSR = r;
rint_nSR = rint;
stats_nSR = stats;

[b,bint,r,rint,stats] = regress(YY_AF',New_X_AF);
stats;
b_AF = b;
bint_AF = bint;
r_AF = r;
rint_AF = rint;
stats_AF = stats;

n_SR_bar = b_nSR(2:end);
AF_bar = b_AF(2:end);
n_SR_AF = [n_SR_bar AF_bar];
barWidth=0.8;

figure(1)
bh = barh([n_SR_AF(:,1) n_SR_AF(:,2)], 'BarWidth', barWidth); hold on;
bh(1).FaceColor = 'k';
bh(2).FaceColor = 'r';
legend ('nSR', 'AF')
legend boxoff
set(gcf,'color','w')

[~,N_pars] = size(X);

set(gca,'YTick',1:N_pars)
set(gca,'YTickLabel',parameter_names, 'ytick', 1:(length(parameter_names)))

box off

for i = 1:N_pars
	plot ([min(b(i+1),0)-0.02, -2], [i i], ':','Color',[0.5 0.5 0.5],'linewidth',1.5, 'HandleVisibility','off'); 
end

set(gca, 'XTickLabelRotation', 90)
set(gcf, 'units','inch','position', [1,1,5.25,8])
set(gca,'Ydir','reverse')
title('Linear Regression Analysis APD90')
box off
set(gca,'YAxisLocation', 'left', 'XAxisLocation', 'bottom')
xlim([-0.6 0.4])
ax = gca;
ax.TickDir = 'out';

edges = [160:3:320];

figure(2)
hold on
m = histogram(pacing_threshold_APD90,'BinEdges',edges);
m(1).FaceColor = 'k';
l = histogram(pacing_threshold_AF_APD90,'BinEdges',edges);
l(1).FaceColor = 'r';
set(gcf, 'units','pixels','position', [93,97,500,300])
set(gcf,'color','w')

ax = gca;
ax.TickDir = 'out';
title('APD90 Population') 
xlabel('APD90 (ms)')
ylabel('Number of cells')
ylim([0 150])
legend('nSR', 'AF')
legend boxoff

%% Plot 2
rm_IKACH_parameters = all_parameters;  
rm_IKACH_parameters_AF = all_parameters; 

X = log(rm_IKACH_parameters);
Y = log(pacing_threshold_CaTamp'); 
X_AF = log(rm_IKACH_parameters_AF);
Y_AF = log(pacing_threshold_CaTamp_AF'); 

sel_index = ~ isnan(Y);
sel_index_AF = ~ isnan(Y_AF);

XX = zscore(X(sel_index,:));
YY = zscore(Y(sel_index));

XX_AF = zscore(X_AF(sel_index_AF,:));
YY_AF = zscore(Y_AF(sel_index_AF));

[nn, mm ] = size(XX);
[nn_AF, mm_AF ] = size(XX_AF);

New_X = [ones(nn, 1), XX];
New_X_AF = [ones(nn_AF, 1), XX_AF];

[b,bint,r,rint,stats] = regress(YY',New_X);
stats;
b_nSR = b;
bint_nSR = bint;
r_nSR = r;
rint_nSR = rint;
stats_nSR = stats;

[b,bint,r,rint,stats] = regress(YY_AF',New_X_AF);
stats;
b_AF = b;
bint_AF = bint;
r_AF = r;
rint_AF = rint;
stats_AF = stats;

n_SR_bar = b_nSR(2:end);
AF_bar = b_AF(2:end);
n_SR_AF = [n_SR_bar AF_bar];
barWidth=0.8;

figure(3)
bh = barh([n_SR_AF(:,1) n_SR_AF(:,2)], 'BarWidth', barWidth); hold on;
bh(1).FaceColor = 'k';
bh(2).FaceColor = 'r';
legend ('nSR', 'AF')
legend boxoff
set(gcf,'color','w')

[~,N_pars] = size(X);

set(gca,'YTick',1:N_pars)
set(gca,'YTickLabel',parameter_names, 'ytick', 1:(length(parameter_names)))

box off

for i = 1:N_pars
	plot ([min(b(i+1),0)-0.02, -2], [i i], ':','Color',[0.5 0.5 0.5],'linewidth',1.5, 'HandleVisibility','off'); 
end

set(gca, 'XTickLabelRotation', 90)
set(gcf, 'units','inch','position', [7.3854166,1.04166,5.25,8])
set(gca,'Ydir','reverse')
title('Linear Regression Analysis CaT Amp')
box off
set(gca,'YAxisLocation', 'left', 'XAxisLocation', 'bottom')
xlim([-0.5 0.8])
ax = gca;
ax.TickDir = 'out';

edges = [0:10:600];
figure(4)
hold on
m = histogram(pacing_threshold_CaTamp*1000000,'BinEdges',edges);
m(1).FaceColor = 'k';
l = histogram(pacing_threshold_CaTamp_AF*1000000,'BinEdges',edges);
l(1).FaceColor = 'r';
set(gcf, 'units','pixels','position', [710,103,500,300])
set(gcf,'color','w')

ax = gca;
ax.TickDir = 'out';
title('CaT Amp Population')
xlabel('CaT Amp (nM)')
ylabel('Number of cells')
legend('nSR', 'AF')
legend boxoff

%% Plot 3
rm_IKACH_parameters = all_parameters;  
rm_IKACH_parameters_AF = all_parameters; 

X = log(rm_IKACH_parameters);
Y = log(pacing_threshold_ERP'); 
X_AF = log(rm_IKACH_parameters_AF);
Y_AF = log(pacing_threshold_AF_ERP'); 

sel_index = ~ isnan(Y);
sel_index_AF = ~ isnan(Y_AF);

XX = zscore(X(sel_index,:));
YY = zscore(Y(sel_index));

XX_AF = zscore(X_AF(sel_index_AF,:));
YY_AF = zscore(Y_AF(sel_index_AF));

[nn, mm ] = size(XX);
[nn_AF, mm_AF ] = size(XX_AF);

New_X = [ones(nn, 1), XX];
New_X_AF = [ones(nn_AF, 1), XX_AF];

[b,bint,r,rint,stats] = regress(YY',New_X);
stats;
b_nSR = b;
bint_nSR = bint;
r_nSR = r;
rint_nSR = rint;
stats_nSR = stats;

[b,bint,r,rint,stats] = regress(YY_AF',New_X_AF);
stats;
b_AF = b;
bint_AF = bint;
r_AF = r;
rint_AF = rint;
stats_AF = stats;

n_SR_bar = b_nSR(2:end);
AF_bar = b_AF(2:end);
n_SR_AF = [n_SR_bar AF_bar];
barWidth=0.8;

figure(5)
bh = barh([n_SR_AF(:,1) n_SR_AF(:,2)], 'BarWidth', barWidth); hold on;
bh(1).FaceColor = 'k';
bh(2).FaceColor = 'r';
legend ('nSR', 'AF')
legend boxoff
set(gcf,'color','w')

[~,N_pars] = size(X);

set(gca,'YTick',1:N_pars)
set(gca,'YTickLabel',parameter_names, 'ytick', 1:(length(parameter_names)))

box off

for i = 1:N_pars
	plot ([min(b(i+1),0)-0.02, -2], [i i], ':','Color',[0.5 0.5 0.5],'linewidth',1.5, 'HandleVisibility','off'); 
end

set(gca, 'XTickLabelRotation', 90)
set(gcf, 'units','inch','position', [14.6979166,0.979,5.2499,7.9999])
set(gca,'Ydir','reverse')
title('Linear Regression Analysis ERP')
box off
set(gca,'YAxisLocation', 'left', 'XAxisLocation', 'bottom')
xlim([-0.6 0.4])
ax = gca;
ax.TickDir = 'out';

edges = [150:2.5:350];
figure(6)
hold on
m = histogram(pacing_threshold_ERP,'BinEdges',edges);
m(1).FaceColor = 'k';
l = histogram(pacing_threshold_AF_ERP,'BinEdges',edges);
l(1).FaceColor = 'r';
set(gcf, 'units','pixels','position', [1411,97,500,293])
set(gcf,'color','w')

ax = gca;
ax.TickDir = 'out';
title('ERP Population')
xlabel('ERP (ms)')
ylabel('Number of cells')
legend('nSR', 'AF')
legend boxoff
