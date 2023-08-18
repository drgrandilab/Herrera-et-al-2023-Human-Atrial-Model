%% Figure 3 Linear Regression Alternans
close all
clear
clc

addpath Conditions

%% Load Perturbations
load Parameter_matrix.mat % sigma 0.1
Number_of_cells = 1000; % Always 1000

%% Alternans
%AF
load Alternans_Population_AF.mat
pacing_threshold_AF = cell2mat(pacing_threshold)';
logistic_pacing_threshold_AF = pacing_threshold_AF;
logistic_parameters_AF = all_parameters;
all_parameters_AF = all_parameters;

%nSR
load Alternans_Population_nSR.mat  

%% Inital variables
pacing_threshold = cell2mat(pacing_threshold)';
logistic_pacing_threshold = pacing_threshold;
logistic_parameters = all_parameters;

%% Linear Regression
parameter_names = [{'GNa'},{'GClCa'},{'GCaL'},{'Gtof'},{'GKur'},{'GKr'},{'GKs'},{'GK1'},{'VNCX'},{'VNKA'}, {'GK2P'},{'GNaB'},{'GCaB'}, ...
{'GClB'},{'GKAch'},{'GNaL'},{'GKp'},{'GSK'},{'VPMCA'},{'VSERCA'},{'VRyR'},{'VRyRLeak'}];

%nSR
empty_spots = find(pacing_threshold == 2); %Indicates no DAD above 100 BCL
number_of_below = numel(empty_spots);
empty_spots_above = find(pacing_threshold == 4);  % 4 Indicates DAD above 1000 BCL
number_of_above = numel(empty_spots_above);
bar_all_parameters = all_parameters;
bar_pacing_threshold = pacing_threshold;

Array_one = [empty_spots; empty_spots_above]';
Asorted = sort(Array_one,'ascend');
all_parameters(Asorted,:) = []; 
pacing_threshold(Asorted,:) = [];

bar_all_parameters(empty_spots,:) = 1; 
bar_pacing_threshold(empty_spots,:) = 1; 
bar_empty_spots_above = find(bar_pacing_threshold == 4);
bar_all_parameters(empty_spots_above,:) = 1000; 
bar_pacing_threshold(empty_spots_above,:) = 1000; 

% AF
empty_spots_AF = find(pacing_threshold_AF == 2);
number_of_below_AF = numel(empty_spots_AF);
empty_spots_above_AF = find(pacing_threshold_AF == 4);
number_of_above_AF = numel(empty_spots_above_AF);
bar_all_parameters_AF = all_parameters_AF;
bar_pacing_threshold_AF = pacing_threshold_AF;

Array_one_AF = [empty_spots_AF; empty_spots_above_AF]';
Asorted_AF = sort(Array_one_AF,'ascend');
all_parameters_AF(Asorted_AF,:) = []; 
pacing_threshold_AF(Asorted_AF,:) = []; 

bar_all_parameters_AF(empty_spots_AF,:) = 1;
bar_pacing_threshold_AF(empty_spots_AF,:) = 1; 
bar_empty_spots_above_AF = find(bar_pacing_threshold_AF == 4);
bar_all_parameters_AF(empty_spots_above_AF,:) = 1000;
bar_pacing_threshold_AF(empty_spots_above_AF,:) = 1000;

%% Regression
%nSR
rm_IKACH_parameters = all_parameters; 
X = log(rm_IKACH_parameters);
Y = log(pacing_threshold'); 
sel_index = ~ isnan(Y);

XX = zscore(X(sel_index,:));
YY = zscore(Y(sel_index));
[nn, mm ] = size(XX);
New_X = [ones(nn, 1), XX];
[b,bint,r,rint,stats] = regress(YY',New_X);
stats;
n_SR_bar = b(2:end);

%AF
rm_IKACH_parameters_AF = all_parameters_AF; 
X_AF = log(rm_IKACH_parameters_AF);
Y_AF = log(pacing_threshold_AF'); 
sel_index_AF = ~ isnan(Y_AF);

XX_AF = zscore(X_AF(sel_index_AF,:));
YY_AF = zscore(Y_AF(sel_index_AF));
[nn_AF, mm_AF ] = size(XX_AF);
New_X_AF = [ones(nn_AF, 1), XX_AF];
[b_AF,bint_AF,r_AF,rint_AF,stats_AF] = regress(YY_AF',New_X_AF);
stats_AF;
n_SR_bar_AF = b_AF(2:end);
n_SR_plot = [n_SR_bar n_SR_bar_AF];
barWidth=0.8;

str_1 = 'nSR Cells = ' + string(numel(pacing_threshold));
str_2 = 'AF Cells = ' + string(numel(pacing_threshold_AF));
str_3 = 'R^2 nSR = ' + string(stats(1));
str_4 = 'R^2 AF = ' + string(stats_AF(1));

figure(5)
bh = barh([n_SR_plot(:,1) n_SR_plot(:,2)], 'BarWidth', barWidth); hold on;
bh(1).FaceColor = 'k';
bh(2).FaceColor = 'r';
legend ('nSR', 'AF')
legend boxoff
set(gcf,'units','inch', 'position', [0.0625,1.010416666666667,5.25,8])
text(0.25,19,str_1)
text(0.25,18.5,str_2)
text(0.25,18,str_3)
text(0.25,17.5,str_4)
ax = gca;
ax.TickDir = 'out';

[~,N_pars] = size(X);

set(gca,'YTick',1:N_pars)
set(gca,'YTickLabel',parameter_names, 'ytick', 1:(length(parameter_names)))
box off

for i = 1:N_pars
	plot ([min(b(i+1),0)-0.02, -2], [i i], ':','Color',[0.5 0.5 0.5],'linewidth',1.5, 'HandleVisibility','off'); 
end
set(gca, 'XTickLabelRotation', 90)
set(gcf,'color','w')
set(gca,'Ydir','reverse')
title('Linear Regression Analysis-Alternans')
xlabel('Regression Coefficients')
box off
set(gca,'YAxisLocation', 'left', 'XAxisLocation', 'bottom')
xlim([-1 1])
