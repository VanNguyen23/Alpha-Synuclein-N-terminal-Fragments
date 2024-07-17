clear, clc
close all
%h = 250;
%% P1
clc
Energy = table2array(readtable('Etotal.dat'));
%Energy2 = table2array(readtable('E_P3L38AR2.dat'));
%Energy3 = table2array(readtable('E_P3L38AR3.dat'));
index1 = find(Energy(:,1) == 100000000);
index2 = find(Energy(:,1) == 150000000);
index3 = find(Energy(:,1) == 200000000);
index4 = find(Energy(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count21 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.186;
numbeads = 26;
numchain = 48;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count21:index(m)
        Energy(n,2) = Energy(n,2)*3.3/sqrt(temps(m)*12);
    end
    count21 = index(m)+1;
end
for l = count21:height(Energy)
    Energy(l,2) = (Energy(l,2)-Energy(count21,2))*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U = (Energy(:,3) - 0.5*Energy(:,4)*3*noptotal)*12.47;
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

Energy2 = table2array(readtable('Etotal2.dat'));
index1 = find(Energy2(:,1) == 100000000);
index2 = find(Energy2(:,1) == 150000000);
index3 = find(Energy2(:,1) == 200000000);
index4 = find(Energy2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count22 = 1;
for m = 1:length(index)
    for n = count22:index(m)
        Energy2(n,2) = Energy2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count22 = index(m)+1;
end
for l = count22:height(Energy2)
    Energy2(l,2) = (Energy2(l,2)-Energy2(count22,2))*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2 = (Energy2(:,3) - 0.5*Energy2(:,4)*3*noptotal)*12.47;
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

Energy3 = table2array(readtable('Etotal3.dat'));
index1 = find(Energy3(:,1) == 100000000);
index2 = find(Energy3(:,1) == 150000000);
index3 = find(Energy3(:,1) == 200000000);
index4 = find(Energy3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count23 = 1;
for m = 1:length(index)
    for n = count23:index(m)
        Energy3(n,2) = Energy3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count23 = index(m)+1;
end
for l = count23:height(Energy3)
    Energy3(l,2) = (Energy3(l,2)-Energy3(count23,2))*3.3/sqrt(setemp*12);
end
%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3 = (Energy3(:,3) - 0.5*Energy3(:,4)*3*noptotal)*12.47;
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

Energy4 = table2array(readtable('Etotal4.dat'));
index1 = find(Energy4(:,1) == 100000000);
index2 = find(Energy4(:,1) == 150000000);
index3 = find(Energy4(:,1) == 200000000);
index4 = find(Energy4(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count24 = 1;
for m = 1:length(index)
    for n = count24:index(m)
        Energy4(n,2) = Energy4(n,2)*3.3/sqrt(temps(m)*12);
    end
    count24 = index(m)+1;
end
for l = count24:height(Energy4)
    Energy4(l,2) = (Energy4(l,2)-Energy4(count24,2))*3.3/sqrt(setemp*12);
end
%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U4 = (Energy4(:,3) - 0.5*Energy4(:,4)*3*noptotal)*12.47;
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

h = plot(Energy(count21:end,2),U(count21:end),'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlim([0,35])
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Total Interaction Energy,kJ/mol','FontSize',20,'FontWeight','bold')

title('P3NExtend-WT','FontSize',20,'FontWeight','bold')
hold on
plot(Energy2(count22:end,2),U2(count22:end),'LineWidth',2)
plot(Energy3(count23:end,2),U3(count23:end),'LineWidth',2)
plot(Energy4(count24:end,2),U4(count24:end),'LineWidth',2)
legend('P1-WT','P1-L38A','P1-V40A','P1-L38M','FontSize',18)
hold off
% figure(3)
% % h2 = 20;
% % tEavg2 = arrayfun(@(i) mean(timeE(i:i+h2-1)),1:h:height(timeE)-h2+1)';
% % Uavg2 = arrayfun(@(i) mean(U(i:i+h2-1)),1:h:height(U)-h2+1)';
% plot(tEavg1,Uavg1,tEavg2,Uavg2,tEavg3,Uavg3)
% title('Potential Energy vs Time - P3L38A','FontSize',23,'FontWeight','bold')
% xlabel('Time, microsecond','FontSize',23,'FontWeight','bold')
% ylabel('U,kJ/mol','FontSize',23,'FontWeight','bold')

% %% L38A
% clc
% L38A1 = table2array(readtable('L38AR1.dat'));
% index1 = find(L38A1(:,1) == 100000000);
% index2 = find(L38A1(:,1) == 150000000);
% index3 = find(L38A1(:,1) == 200000000);
% index4 = find(L38A1(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count21 = 1;
% temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
% setemp = 0.195;
% numbeads = 119;
% numchain = 24;
% noptotal = numbeads*numchain;
% for m = 1:length(index)
%     for n = count21:index(m)
%         L38A1(n,2) = L38A1(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count21 = index(m)+1;
% end
% for l = count21:height(L38A1)
%     L38A1(l,2) = (L38A1(l,2)-L38A1(count21,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U_L38A1 = (L38A1(:,3) - 0.5*L38A1(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% L38A2 = table2array(readtable('L38AR2.dat'));
% index1 = find(L38A2(:,1) == 100000000);
% index2 = find(L38A2(:,1) == 150000000);
% index3 = find(L38A2(:,1) == 200000000);
% index4 = find(L38A2(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count22 = 1;
% for m = 1:length(index)
%     for n = count22:index(m)
%         L38A2(n,2) = L38A2(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count22 = index(m)+1;
% end
% for l = count22:height(L38A2)
%     L38A2(l,2) = (L38A2(l,2)-L38A2(count22,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U2_L38A = (L38A2(:,3) - 0.5*L38A2(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% L38A3 = table2array(readtable('L38AR3.dat'));
% index1 = find(L38A3(:,1) == 100000000);
% index2 = find(L38A3(:,1) == 150000000);
% index3 = find(L38A3(:,1) == 200000000);
% index4 = find(L38A3(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count23 = 1;
% for m = 1:length(index)
%     for n = count23:index(m)
%         L38A3(n,2) = L38A3(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count23 = index(m)+1;
% end
% for l = count23:height(L38A3)
%     L38A3(l,2) = (L38A3(l,2)-L38A3(count23,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U3_L38A = (L38A3(:,3) - 0.5*L38A3(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% figure(2)
% h = plot(L38A1(count21:end,2),U_L38A1(count21:end),'LineWidth',2);
% %ytickformat('%.0f')
% ax = ancestor(h, 'axes');
% ax.YAxis.Exponent = 0;
% set(ax,'FontSize',20)
% xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
% ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
% 
% title('P3NExtend-L38A','FontSize',20,'FontWeight','bold')
% hold on
% plot(L38A2(count22:end,2),U2_L38A(count22:end),'LineWidth',2)
% plot(L38A3(count23:end,2),U3_L38A(count23:end),'LineWidth',2)
% legend('Run 1','Run 2','Run 3','FontSize',18)
% hold off
% 
% %% L38M
% clc
% L38M1 = table2array(readtable('L38MR1.dat'));
% index1 = find(L38M1(:,1) == 100000000);
% index2 = find(L38M1(:,1) == 150000000);
% index3 = find(L38M1(:,1) == 200000000);
% index4 = find(L38M1(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
% setemp = 0.195;
% numbeads = 119;
% numchain = 24;
% noptotal = numbeads*numchain;
% for m = 1:length(index)
%     for n = count2:index(m)
%         L38M1(n,2) = L38M1(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(L38M1)
%     L38M1(l,2) = (L38M1(l,2)-L38M1(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U_L38M1 = (L38M1(:,3) - 0.5*L38M1(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% L38M2 = table2array(readtable('L38MR2.dat'));
% index1 = find(L38M2(:,1) == 100000000);
% index2 = find(L38M2(:,1) == 150000000);
% index3 = find(L38M2(:,1) == 200000000);
% index4 = find(L38M2(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         L38M2(n,2) = L38M2(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(L38M2)
%     L38M2(l,2) = (L38M2(l,2)-L38M2(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U2_L38M = (L38M2(:,3) - 0.5*L38M2(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% L38M3 = table2array(readtable('L38MR3.dat'));
% index1 = find(L38M3(:,1) == 100000000);
% index2 = find(L38M3(:,1) == 150000000);
% index3 = find(L38M3(:,1) == 200000000);
% index4 = find(L38M3(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         L38M3(n,2) = L38M3(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(L38M3)
%     L38M3(l,2) = (L38M3(l,2)-L38M3(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U3_L38M = (L38M3(:,3) - 0.5*L38M3(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% figure(3)
% h = plot(L38M1(count2:end,2),U_L38M1(count2:end),'LineWidth',2);
% %ytickformat('%.0f')
% ax = ancestor(h, 'axes');
% ax.YAxis.Exponent = 0;
% set(ax,'FontSize',20)
% xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
% ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
% 
% title('P3NExtend-L38M','FontSize',20,'FontWeight','bold')
% hold on
% plot(L38M2(count2:end,2),U2_L38M(count2:end),'LineWidth',2)
% plot(L38M3(count2:end,2),U3_L38M(count2:end),'LineWidth',2)
% legend('Run 1','Run 2','Run 3','FontSize',18)
% hold off
% 
% %% V40A
% clc
% V40A1 = table2array(readtable('V40AR1.dat'));
% index1 = find(V40A1(:,1) == 100000000);
% index2 = find(V40A1(:,1) == 150000000);
% index3 = find(V40A1(:,1) == 200000000);
% index4 = find(V40A1(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
% setemp = 0.195;
% numbeads = 119;
% numchain = 24;
% noptotal = numbeads*numchain;
% for m = 1:length(index)
%     for n = count2:index(m)
%         V40A1(n,2) = V40A1(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(V40A1)
%     V40A1(l,2) = (V40A1(l,2)-V40A1(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U_V40A1 = (V40A1(:,3) - 0.5*V40A1(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% V40A2 = table2array(readtable('V40AR2.dat'));
% index1 = find(V40A2(:,1) == 100000000);
% index2 = find(V40A2(:,1) == 150000000);
% index3 = find(V40A2(:,1) == 200000000);
% index4 = find(V40A2(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         V40A2(n,2) = V40A2(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(V40A2)
%     V40A2(l,2) = (V40A2(l,2)-V40A2(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U2_V40A = (V40A2(:,3) - 0.5*V40A2(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% 
% V40A3 = table2array(readtable('V40AR3.dat'));
% index1 = find(V40A3(:,1) == 100000000);
% index2 = find(V40A3(:,1) == 150000000);
% index3 = find(V40A3(:,1) == 200000000);
% index4 = find(V40A3(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         V40A3(n,2) = V40A3(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = index(m)+1;
% end
% for l = count2:height(V40A3)
%     V40A3(l,2) = (V40A3(l,2)-V40A3(count2,2))*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U3_V40A = (V40A3(:,3) - 0.5*V40A3(:,4)*3*noptotal)*12.47;
% %Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
% figure(4)
% h = plot(V40A1(count2:end,2),U_V40A1(count2:end),'LineWidth',2);
% %ytickformat('%.0f')
% ax = ancestor(h, 'axes');
% ax.YAxis.Exponent = 0;
% set(ax,'FontSize',20)
% xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
% ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
% xlim([0 70])
% title('P3NExtend-V40A','FontSize',20,'FontWeight','bold')
% hold on
% plot(V40A2(count2:end,2),U2_V40A(count2:end),'LineWidth',2)
% plot(V40A3(count2:end,2),U3_V40A(count2:end),'LineWidth',2)
% legend('Run 1','Run 2','Run 3','FontSize',18)
% hold off
