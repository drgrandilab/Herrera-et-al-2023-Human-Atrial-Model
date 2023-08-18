function outputs = function_logisitic_regression(logistic_parameters, logistic_pacing_threshold, logistic_parameters_AF, logistic_pacing_threshold_AF)

%% Establish Vectors
% nSR
all_parameters = logistic_parameters;
allpars_LOGISTIC = all_parameters;
X_LOG = log(allpars_LOGISTIC) ;
[N_trials N_pars] = size(all_parameters);
temp_logistic_pacing_threshold = logistic_pacing_threshold';
find_twos = find(temp_logistic_pacing_threshold == 2); %DAD
find_ones = find(temp_logistic_pacing_threshold == 1); %no DAD
logistic_pacing_threshold(logistic_pacing_threshold==2) = 3;
logistic_pacing_threshold(logistic_pacing_threshold==1) = 2;
logistic_pacing_threshold(logistic_pacing_threshold==3) = 1;

for ii=1:N_pars % z-score
    X_LOGISTIC(:,ii)=(X_LOG(:,ii)-mean(X_LOG(:,ii)))/std(X_LOG(:,ii));
end

[B_LOGISTIC,dev,stats] = mnrfit (X_LOGISTIC,logistic_pacing_threshold);

%% DAD Probability
%nSR
B0 = B_LOGISTIC(1);

% PEAD in each model of the population
P_DAD_array = zeros(1,N_trials);
for iii=1:N_trials
    P_DAD_array(iii) = 1/(1+exp(-(B0+sum(B_LOGISTIC(2:end).*X_LOGISTIC(iii,:)'))));
end

all_outputs_DAD_sum = (logistic_pacing_threshold); % 1 x 1000
DAD_presence = (all_outputs_DAD_sum'>1); % 1 for EAD occurrence, 0 for no EAD

array_1 = P_DAD_array(DAD_presence>0.5); % EAD
array_1_mean = mean(array_1);
array_0 = P_DAD_array(DAD_presence<0.5); % no EAD
array_0_mean = mean(array_0);

% Tjur (2009)
R2logistic = abs(array_1_mean - array_0_mean);

% AF
all_parameters_AF = logistic_parameters_AF;
allpars_LOGISTIC_AF = all_parameters_AF;
X_LOG_AF = log(allpars_LOGISTIC_AF) ;
[N_trials_AF N_pars_AF] = size(all_parameters_AF);
temp_logistic_pacing_threshold_AF = logistic_pacing_threshold_AF';
find_twos_AF = find(temp_logistic_pacing_threshold_AF == 2); %DAD
find_ones_AF = find(temp_logistic_pacing_threshold_AF == 1); %no DAD
logistic_pacing_threshold_AF(logistic_pacing_threshold_AF==2) = 3;
logistic_pacing_threshold_AF(logistic_pacing_threshold_AF==1) = 2;
logistic_pacing_threshold_AF(logistic_pacing_threshold_AF==3) = 1;

for ii=1:N_pars_AF % z-score
    X_LOGISTIC_AF(:,ii)=(X_LOG_AF(:,ii)-mean(X_LOG_AF(:,ii)))/std(X_LOG_AF(:,ii));
end

[B_LOGISTIC_AF,dev,stats] = mnrfit (X_LOGISTIC_AF,logistic_pacing_threshold_AF);

%% Plots
minimum = min(B_LOGISTIC);
minimum = round(minimum, -2);
if minimum < 10
    minimum = min(B_LOGISTIC);
    minimum = round(minimum, -1);
end
minimum = -350;

%% Rotate
parameter_names_logisitc = [{'b0'},{'G_N_a'},{'G_C_l_C_a'},{'G_C_a_L'},{'G_t_o_f'},{'G_K_u_r'},{'G_K_r'},{'G_K_s'},{'G_K_1'},{'V_N_C_X'},{'V_N_K_A'}, {'G_K_2_P'},{'G_N_a_B'},{'G_C_a_B'}, ...
{'G_C_l_B'},{'G_K_A_c_h'},{'G_N_a_L'},{'G_K_p'},{'G_S_K'},{'V_P_M_C_A'},{'V_S_E_R_C_A'},{'V_R_y_R'},{'V_R_y_R_L_e_a_k'}];

set(gcf,'color','w') % with b0
hold on;
barWidth=0.8;

n_SR_AF = [B_LOGISTIC B_LOGISTIC_AF];
figure(1)
bh = barh([n_SR_AF(:,1) n_SR_AF(:,2)], 'BarWidth', barWidth); hold on;
bh(1).FaceColor = 'k';
bh(2).FaceColor = 'r';
set(gcf,'units','inch', 'position', [5.354166666666666,1.020833333333333,5.25,8.020833333333332])
N_pars = 22;
set(gca,'YTick',1:N_pars)
legend('nSR', 'AF')
set(gca,'YTickLabel',parameter_names_logisitc, 'ytick', 1:(length(parameter_names_logisitc)))
legend boxoff

for i = 1:N_pars
	plot ([min(B_LOGISTIC(i+1),0)-0.02, minimum], [i+1 i+1], ':','Color',[0.5 0.5 0.5],'linewidth',1.5, 'HandleVisibility','off'); 
end
set(gca, 'XTickLabelRotation', 90)
set(gca,'Ydir','reverse')
title ('Logistic Regression Analysis-DADs')
xlabel('Regression Coefficients')

box off
set(gca,'YAxisLocation', 'left', 'XAxisLocation', 'bottom')
xlim([-70 60])

%% AF
B0_AF = B_LOGISTIC_AF(1);
P_ead_array_AF = zeros(1,N_trials_AF);
for iii=1:N_trials_AF
    P_ead_array_AF(iii) = 1/(1+exp(-(B0_AF+sum(B_LOGISTIC_AF(2:end).*X_LOGISTIC_AF(iii,:)'))));
end

all_outputs_DAD_sum_AF = (logistic_pacing_threshold_AF); % 1 x 1000
DAD_presence_AF = (all_outputs_DAD_sum_AF'>1); % 1 for EAD occurrence, 0 for no EAD

array_1_AF = P_ead_array_AF(DAD_presence_AF>0.5); % EAD
array_1_mean_AF = mean(array_1_AF);
array_0_AF = P_ead_array_AF(DAD_presence_AF<0.5); % no EAD
array_0_mean_AF = mean(array_0_AF);

% Tjur (2009)
R2logistic_AF = abs(array_1_mean_AF - array_0_mean_AF);

number_of_DADs_AF_above_BCL = find(1 == logistic_pacing_threshold_AF);
number_of_DADs_AF_above_BCL = numel(number_of_DADs_AF_above_BCL);

number_of_DADs_nSR_above_BCL = find(1 == logistic_pacing_threshold);
number_of_DADs_nSR_above_BCL = numel(number_of_DADs_nSR_above_BCL);

str_1 = 'R^2 nSR = ' + string(R2logistic);
str_2 = 'R^2 AF = ' + string(R2logistic_AF);

str_3 = 'Cells with DADs nSR = ' + string(number_of_DADs_nSR_above_BCL);
str_4 = 'Cells with DADs AF = ' + string(number_of_DADs_AF_above_BCL);

figure(1)
text(10,20,str_1)
text(10,19.5,str_2)
text(10,19,str_3)
text(10,18.5,str_4)
