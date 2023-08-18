%% Calculate timecourse for currents and other intermediates
function currents = calcCurrents(BLC_plot_array_thyme,BLC_plot_array_y,p)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents: [I_Na I_Catot];
currents=[];
for i=1:size(BLC_plot_array_thyme)
    if ceil(i/1000)==i/1000
%         disp(['t = ',num2str(ceil(BLC_plot_array_thyme(i)))]);
    end
    currents=[currents;Herrera_2023_model_SK_expression(BLC_plot_array_thyme(i),BLC_plot_array_y(i,:),p,'currents')];
end
% end calcCurrents