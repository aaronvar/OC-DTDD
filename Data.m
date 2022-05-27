n_runs=4000;

powers=21:-10:1;
powerdata_test_low=zeros(4,length(powers));
for change_variable = 1:length(powers)
    ind_var=powers(change_variable);
    SumRateData1=zeros(n_runs,1);
    SumRateData2=zeros(n_runs,1);
    SumRateData3=zeros(n_runs,1);
    SumRateData4=zeros(n_runs,1);
    for run_number=1:n_runs
        run("Starter.m")
        SumRateData1(run_number)=SumRate1;
        SumRateData2(run_number)=SumRate2;
        SumRateData3(run_number)=SumRate3;
        SumRateData4(run_number)=SumRate4;
    end
    
    powerdata_test_low(1,change_variable)=mean(SumRateData1);
    powerdata_test_low(2,change_variable)=mean(SumRateData2);
    powerdata_test_low(3,change_variable)=mean(SumRateData3);
    powerdata_test_low(4,change_variable)=mean(SumRateData4);
end

% plot(powers, allpowerdata_new,'--*')
% legend('Centralised','Conventional','Random','Duplex','Location','northwest')
% xlabel('Transmit Power (dBm)')
% ylabel('Sum Rate (bits/symbol)')
% title('Sum Rate vs Transmit Power of each scheme')
% 
% plot(sidelengths, sidelengthdata_test,'--*')
% legend('Centralised','Conventional','Random','Duplex')
% xlabel('Maximum Sidelength - rho (m)')
% ylabel('Sum Rate')
% title('Sum Rate vs Maximum Sidelength n=100')
% 
% plot(SIsuppressions, SIsuppressiondata_test,'--*')
% legend('Centralised','Conventional','Random','Duplex')
% xlabel('SI Suppression (dB)')
% ylabel('Sum Rate')
% title('Sum Rate vs SI Suppression n=100')