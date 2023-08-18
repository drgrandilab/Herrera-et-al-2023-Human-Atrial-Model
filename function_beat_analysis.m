function outputs = function_beat_analysis_2017(time,Vm,Ca,Na,dVm,period,AP_index)
% Beat analysis first or last AP (index = 1 or 2)

%% Check
if AP_index == 1,
    t1 = 0; t2 = t1+period-1.5; t3 = t1+2*period-1.5;
else
    t1 = time(end)-2*period; t2 = t1+period-1.5; t3 = t1+2*period-1.5;
end

% Alternans
t1_roi = find(time>t1); t1_index = t1_roi(1); Vm1 = Vm(t1_index);
t2_roi = find(time>t2); t2_index = t2_roi(1); Vm2 = Vm(t2_index);
t3_roi = find(time>t3); t3_index = t3_roi(1); Vm3 = Vm(t3_index);
Vmin1 = min(Vm(t1_index:t2_index));
Vmin2 = min(Vm(t2_index:t3_index));
%[abs(Vm1-Vm2) abs(Vm2-Vm3) abs(Vm1-Vm3) abs(Vm2-Vmin1) abs(Vm2-Vmin2)]
delta = [abs(Vm1-Vm2) abs(Vm2-Vm3) abs(Vm1-Vm3) abs(Vm2-Vmin1) abs(Vm2-Vmin2)];

% Triggered APs
t4_roi = find(time>t2+15); t4_index = t4_roi(1);
taps = find(dVm(t4_index:t3_index)>25);
%figure, plot(time(t4_index:t3_index),dVm(t4_index:t3_index))

if max(delta) > 5 %| length(taps) > 0,
    outputs = zeros(1,15);
else

    %% ROI
    if AP_index == 1,
        t_in = 0; t_fin = t_in+period-5;
    else
        t_in = time(end)-period; t_fin = t_in+period-5;
    end

    t_in_roi = find(time>t_in); t_in_index = t_in_roi(1);%-1;
    t_fin_roi = find(time>t_fin); t_fin_index = t_fin_roi(1);

    %% Em
    [dVm_max, index] = max(dVm(t_in_index:t_fin_index)); index = index-1; % max slope

    [Vm_max, index_max] = max(Vm(t_in_index:t_fin_index)); index_max = index_max-1; % peak
    Vm_min = min(Vm(t_in_index:t_fin_index)); % Em resting
    AP_amp = Vm_max-Vm_min;

    APfraction = 0.90;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD90 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.7;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD70 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.5;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    start_APD = time(t_in_index+index_max+APD_roi(1)); index_start = find (time == start_APD);
    vm_start = Vm(index_start);
    end_APD = time(t_in_index+index); index_end = find (time == end_APD);
    vm_end = Vm(index_end);
    APD50 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.3;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD30 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    % Ca
    [Ca_max, index_ca] = max(Ca(t_in_index:t_fin_index)); index_ca = index_ca-1; % peak CaT
    ca_max_index = find(Ca == Ca_max); time_max_ca = time(ca_max_index);
%     figure(2)
%     plot( time_max_ca,Ca_max, 'o')
    Ca_min = min(Ca(t_in_index:t_fin_index)); % diast Ca
    CaT_amp = Ca_max-Ca_min;
    CaT_rise = time(t_in_index+index_ca)-time(t_in_index);

    Cafraction = 0.5;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_50 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);

    Cafraction = 0.632;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_63 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);

    %% Na
    [Na_min, index_na] = min(Na(t_in_index:t_fin_index)); index_na = index_na-1; % diast Na

    %% Figures-Nate
%     start_APD = time(t_in_index+index_max+APD_roi(1)); end_APD = time(t_in_index+index); index_start_time = find(start_APD == time); index_end_time = find(start_APD == time);
%     vm_apd90_start = Vm(index_start_time); vm_apd90_end = Vm(index_end_time);
%     Vm_min_time = find(Vm == Vm_min); time_vm_min = time(Vm_min_time);
%     Ca_max_index = find(Ca == Ca_max); Time_Ca_peak = time(Ca_max_index);
%     Ca_min_index = find(Ca == Ca_min); Time_Ca_min = time(Ca_min_index);



%     figure(1)
%     subplot(2,1,1)
%     plot (start_APD, vm_start, 'bo')
%     hold on,
%     plot (end_APD, vm_end, 'bo')
%     plot (time_vm_min, Vm_min, 'rx')
%     subplot(3,1,2),
%     plot (Time_Ca_peak, Ca_max ,'bo')
%     plot (Time_Ca_min, Ca_min ,'bo')

    %% Collect all outputs
    outputs = [dVm_max Vm_max -Vm_min AP_amp APD90 APD70 APD50 APD30...
        Ca_max Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min];
%     outputs = [dVm_max Vm_max -Vm_min AP_amp APD90 APD70 APD50 APD30 Na_min];

end

% %% Plot
% figure, set(gcf,'color','w')
% subplot(2,2,1)
% hold on, plot(time,dVm,'k')
% plot(time(t_in_index:t_fin_index),dVm(t_in_index:t_fin_index),'r')
% plot(time(t_in_index+index),dVm_max,'*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('dEm (mV/ms)')
% xlim([time(t_in_index) time(t_in_index)+period])
%
% subplot(2,2,3)
% hold on,plot(time,Vm,'k')
% plot(time(t_in_index:t_fin_index),Vm(t_in_index:t_fin_index),'r')
% plot(time(t_in_index+index_max),Vm_max,'*')
% plot([time(t_in_index) time(t_in_index)+period],[Vm_min Vm_min],':')
% plot([time(t_in_index) time(t_in_index)+period],[Vm_max-APfraction*AP_amp Vm_max-APfraction*AP_amp],'g:')
% plot(time(t_in_index+index_max+APD_roi(1)),Vm(t_in_index+index_max+APD_roi(1)),'g*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('Em (mV)'),xlabel('Time (ms)')
% xlim([time(t_in_index) time(t_in_index)+period])
%
% subplot(2,2,2)
% hold on, plot(time,Ca,'k')
% plot(time(t_in_index:t_fin_index),Ca(t_in_index:t_fin_index),'r')
% plot(time(t_in_index+index_ca),Ca_max,'*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('[Ca] (mM)')
% xlim([time(t_in_index) time(t_in_index)+period])
%
% subplot(2,2,4)
% hold on, plot(time,Na,'k')
% plot(time(t_in_index:t_fin_index),Na(t_in_index:t_fin_index),'r')
% plot(time(t_in_index+index_na),Na_min,'*')
% set(gca,'box','off','tickdir','out','fontsize',12)
% ylabel('[Na] (mM)'),xlabel('Time (ms)')
% xlim([time(t_in_index) time(t_in_index)+period])
