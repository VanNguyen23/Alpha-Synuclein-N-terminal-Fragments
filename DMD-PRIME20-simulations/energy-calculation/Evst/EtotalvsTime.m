clear, clc
close all
%h = 250;
%% WT
clc
Energy = table2array(readtable('Etotal.dat'));
%Energy2 = table2array(readtable('E_P3L38AR2.dat'));
%Energy3 = table2array(readtable('E_P3L38AR3.dat'));
index1 = find(Energy(:,1) == 100000000);
index2 = find(Energy(:,1) == 150000000);
index3 = find(Energy(:,1) == 200000000);
index4 = find(Energy(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 84;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        Energy(n,2) = Energy(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(Energy)
    Energy(l,2) = Energy(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U = (Energy(:,3) - 0.5*Energy(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';
U1noanneal = U-U(count2);
Energy2 = table2array(readtable('Etotal2.dat'));
index1 = find(Energy2(:,1) == 100000000);
index2 = find(Energy2(:,1) == 150000000);
index3 = find(Energy2(:,1) == 200000000);
index4 = find(Energy2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        Energy2(n,2) = Energy2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(Energy2)
    Energy2(l,2) = Energy2(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2 = (Energy2(:,3) - 0.5*Energy2(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

Energy3 = table2array(readtable('Etotal3.dat'));
index1 = find(Energy3(:,1) == 100000000);
index2 = find(Energy3(:,1) == 150000000);
index3 = find(Energy3(:,1) == 200000000);
index4 = find(Energy3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        Energy3(n,2) = Energy3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(Energy3)
    Energy3(l,2) = Energy3(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3 = (Energy3(:,3) - 0.5*Energy3(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

Energy4 = table2array(readtable('Etotal4.dat'));
index1 = find(Energy4(:,1) == 100000000);
index2 = find(Energy4(:,1) == 150000000);
index3 = find(Energy4(:,1) == 200000000);
index4 = find(Energy4(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        Energy4(n,2) = Energy4(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(Energy4)
    Energy4(l,2) = Energy4(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U4 = (Energy4(:,3) - 0.5*Energy4(:,4)*3*noptotal);

h = plot(Energy(count2:end,2)-Energy(count2,2),U1noanneal(count2:end),'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',40)
%xlim([0,40])
xlabel('Time, microsecond','FontSize',40,'FontWeight','bold')
ylabel('Total Potential Energy,kJ/mol','FontSize',40,'FontWeight','bold')
title('P3-WT','FontSize',50,'FontWeight','bold')
hold on
plot(Energy2(:,2),U2,'LineWidth',2)
plot(Energy3(:,2),U3,'LineWidth',2)
plot(Energy4(:,2),U4,'LineWidth',2)
%legend('Run 1','Run 2','Run 3','Run 4','FontSize',30,'LineWidth', 2.7)
[lh,hObj]=legend('Run 1','Run 2','Run 3','Run4');           % return the handles array
lh.FontSize = 38;
hL=findobj(hObj,'type','line');  % get the lines, not text
set(hL,'linewidth',7)            % set their width property
hold off

figure(6)
h = plot(Energy(:,2),Energy(:,4)./12,'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Reduced Temperature','FontSize',20,'FontWeight','bold')
title('P3-WT','FontSize',20,'FontWeight','bold')
ylim([0.15,0.25])
hold on
plot(Energy2(:,2),Energy2(:,4)./12,'LineWidth',2)
plot(Energy3(:,2),Energy3(:,4)./12,'LineWidth',2)
plot(Energy4(:,2),Energy4(:,4)./12,'LineWidth',2)

[~,hObj]=legend('Run 1','Run 2','Run 3','Run 4','FontSize',18);           % return the handles array
hL=findobj(hObj,'type','line');  % get the lines, not text
set(hL,'linewidth',2.7)            % set their width property
hold off

figure(7)
yyaxis left
plot(Energy(:,2),U)
yyaxis right
plot(Energy(:,2),Energy(:,4)./12)
ylim([0.15,0.25])
% figure(3)
% % h2 = 20;
% % tEavg2 = arrayfun(@(i) mean(timeE(i:i+h2-1)),1:h:height(timeE)-h2+1)';
% % Uavg2 = arrayfun(@(i) mean(U(i:i+h2-1)),1:h:height(U)-h2+1)';
% plot(tEavg1,Uavg1,tEavg2,Uavg2,tEavg3,Uavg3)
% title('Potential Energy vs Time - P3L38A','FontSize',23,'FontWeight','bold')
% xlabel('Time, microsecond','FontSize',23,'FontWeight','bold')
% ylabel('U,kJ/mol','FontSize',23,'FontWeight','bold')

%% L38A
clc
L38A1 = table2array(readtable('L38AR1.dat'));
index1 = find(L38A1(:,1) == 100000000);
index2 = find(L38A1(:,1) == 150000000);
index3 = find(L38A1(:,1) == 200000000);
index4 = find(L38A1(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 84;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        L38A1(n,2) = L38A1(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38A1)
    L38A1(l,2) = L38A1(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U_L38A1 = (L38A1(:,3) - 0.5*L38A1(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

L38A2 = table2array(readtable('L38AR2.dat'));
index1 = find(L38A2(:,1) == 100000000);
index2 = find(L38A2(:,1) == 150000000);
index3 = find(L38A2(:,1) == 200000000);
index4 = find(L38A2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38A2(n,2) = L38A2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38A2)
    L38A2(l,2) = L38A2(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2_L38A = (L38A2(:,3) - 0.5*L38A2(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

L38A3 = table2array(readtable('L38AR3.dat'));
index1 = find(L38A3(:,1) == 100000000);
index2 = find(L38A3(:,1) == 150000000);
index3 = find(L38A3(:,1) == 200000000);
index4 = find(L38A3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38A3(n,2) = L38A3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38A3)
    L38A3(l,2) = L38A3(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3_L38A = (L38A3(:,3) - 0.5*L38A3(:,4)*3*noptotal);

L38A4 = table2array(readtable('L38AR4.dat'));
index1 = find(L38A4(:,1) == 100000000);
index2 = find(L38A4(:,1) == 150000000);
index3 = find(L38A4(:,1) == 200000000);
index4 = find(L38A4(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38A4(n,2) = L38A4(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38A4)
    L38A4(l,2) = L38A4(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U4_L38A = (L38A4(:,3) - 0.5*L38A4(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

% L38A5 = table2array(readtable('L38AR5.dat'));
% index1 = find(L38A5(:,1) == 100000000);
% index2 = find(L38A5(:,1) == 150000000);
% index3 = find(L38A5(:,1) == 200000000);
% index4 = find(L38A5(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         L38A5(n,2) = L38A5(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = count2+index(m)+1;
% end
% for l = count2:height(L38A5)
%     L38A5(l,2) = L38A5(l,2)*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U5_L38A = (L38A5(:,3) - 0.5*L38A5(:,4)*3*noptotal);

figure(2)
h = plot(L38A1(:,2),U_L38A1,'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
title('P3-L38A','FontSize',20,'FontWeight','bold')
hold on
plot(L38A2(:,2),U2_L38A,'LineWidth',2)
plot(L38A3(:,2),U3_L38A,'LineWidth',2)
plot(L38A4(:,2),U4_L38A,'LineWidth',2)
%plot(L38A5(:,2),U5_L38A,'LineWidth',2)
legend('Run 1','Run 2','Run 3','Run 4','FontSize',18)
hold off

%% L38M
clc
L38M1 = table2array(readtable('L38MR1.dat'));
index1 = find(L38M1(:,1) == 100000000);
index2 = find(L38M1(:,1) == 150000000);
index3 = find(L38M1(:,1) == 200000000);
index4 = find(L38M1(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 84;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        L38M1(n,2) = L38M1(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38M1)
    L38M1(l,2) = L38M1(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U_L38M1 = (L38M1(:,3) - 0.5*L38M1(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

L38M2 = table2array(readtable('L38MR2.dat'));
index1 = find(L38M2(:,1) == 100000000);
index2 = find(L38M2(:,1) == 150000000);
index3 = find(L38M2(:,1) == 200000000);
index4 = find(L38M2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38M2(n,2) = L38M2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38M2)
    L38M2(l,2) = L38M2(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2_L38M = (L38M2(:,3) - 0.5*L38M2(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

L38M3 = table2array(readtable('L38MR3.dat'));
index1 = find(L38M3(:,1) == 100000000);
index2 = find(L38M3(:,1) == 150000000);
index3 = find(L38M3(:,1) == 200000000);
index4 = find(L38M3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38M3(n,2) = L38M3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38M3)
    L38M3(l,2) = L38M3(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3_L38M = (L38M3(:,3) - 0.5*L38M3(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

L38M4 = table2array(readtable('L38MR4.dat'));
index1 = find(L38M4(:,1) == 100000000);
index2 = find(L38M4(:,1) == 150000000);
index3 = find(L38M4(:,1) == 200000000);
index4 = find(L38M4(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        L38M4(n,2) = L38M4(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38M4)
    L38M4(l,2) = L38M4(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U4_L38M = (L38M4(:,3) - 0.5*L38M4(:,4)*3*noptotal);

figure(3)
h = plot(L38M1(:,2),U_L38M1,'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
title('P3-L38M','FontSize',20,'FontWeight','bold')
hold on
plot(L38M2(:,2),U2_L38M,'LineWidth',2)
plot(L38M3(:,2),U3_L38M,'LineWidth',2)
plot(L38M4(:,2),U4_L38M,'LineWidth',2)
legend('Run 1','Run 2','Run 3','Run4','FontSize',18)
hold off

%% V40A
clc
V40A1 = table2array(readtable('V40AR1.dat'));
index1 = find(V40A1(:,1) == 100000000);
index2 = find(V40A1(:,1) == 150000000);
index3 = find(V40A1(:,1) == 200000000);
index4 = find(V40A1(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 84;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        V40A1(n,2) = V40A1(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(V40A1)
    V40A1(l,2) = V40A1(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U_V40A1 = (V40A1(:,3) - 0.5*V40A1(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

V40A2 = table2array(readtable('V40AR2.dat'));
index1 = find(V40A2(:,1) == 100000000);
index2 = find(V40A2(:,1) == 150000000);
index3 = find(V40A2(:,1) == 200000000);
index4 = find(V40A2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        V40A2(n,2) = V40A2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(V40A2)
    V40A2(l,2) = V40A2(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2_V40A = (V40A2(:,3) - 0.5*V40A2(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

V40A3 = table2array(readtable('V40AR3.dat'));
index1 = find(V40A3(:,1) == 100000000);
index2 = find(V40A3(:,1) == 150000000);
index3 = find(V40A3(:,1) == 200000000);
index4 = find(V40A3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        V40A3(n,2) = V40A3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(V40A3)
    V40A3(l,2) = V40A3(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3_V40A = (V40A3(:,3) - 0.5*V40A3(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

V40A4 = table2array(readtable('V40AR4.dat'));
index1 = find(V40A4(:,1) == 100000000);
index2 = find(V40A4(:,1) == 150000000);
index3 = find(V40A4(:,1) == 200000000);
index4 = find(V40A4(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        V40A4(n,2) = V40A4(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(V40A4)
    V40A4(l,2) = V40A4(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U4_V40A = (V40A4(:,3) - 0.5*V40A4(:,4)*3*noptotal);

% V40A5 = table2array(readtable('V40AR5.dat'));
% index1 = find(V40A5(:,1) == 100000000);
% index2 = find(V40A5(:,1) == 150000000);
% index3 = find(V40A5(:,1) == 200000000);
% index4 = find(V40A5(:,1) == 250000000);
% index = [index1', index2', index3', index4'];
% count2 = 1;
% for m = 1:length(index)
%     for n = count2:index(m)
%         V40A5(n,2) = V40A5(n,2)*3.3/sqrt(temps(m)*12);
%     end
%     count2 = count2+index(m)+1;
% end
% for l = count2:height(V40A5)
%     V40A5(l,2) = V40A5(l,2)*3.3/sqrt(setemp*12);
% end
% 
% %timeE1 = Energy1(1:cutoff1,1);
% %tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
% U5_V40A = (V40A5(:,3) - 0.5*V40A5(:,4)*3*noptotal);

figure(4)
h = plot(V40A1(:,2),U_V40A1,'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')
title('P3-V40A','FontSize',20,'FontWeight','bold')
hold on
plot(V40A2(:,2),U2_V40A,'LineWidth',2)
plot(V40A3(:,2),U3_V40A,'LineWidth',2)
plot(V40A4(:,2),U4_V40A,'LineWidth',2)
%plot(V40A5(:,2),U5_V40A,'LineWidth',2)
legend('Run 1','Run 2','Run 3','Run 4','FontSize',18,'Location','Southwest')
hold off

%% S42A
clc
S42A1 = table2array(readtable('S42AR1.dat'));
index1 = find(S42A1(:,1) == 100000000);
index2 = find(S42A1(:,1) == 150000000);
index3 = find(S42A1(:,1) == 200000000);
index4 = find(S42A1(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 84;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        S42A1(n,2) = S42A1(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = count2+index(m)+1;
end
for l = count2:height(S42A1)
    S42A1(l,2) = S42A1(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U_S42A1 = (S42A1(:,3) - 0.5*S42A1(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

S42A2 = table2array(readtable('S42AR2.dat'));
index1 = find(S42A2(:,1) == 100000000);
index2 = find(S42A2(:,1) == 150000000);
index3 = find(S42A2(:,1) == 200000000);
index4 = find(S42A2(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        S42A2(n,2) = S42A2(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = count2+index(m)+1;
end
for l = count2:height(S42A2)
    S42A2(l,2) = S42A2(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U2_S42A = (S42A2(:,3) - 0.5*S42A2(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

S42A3 = table2array(readtable('S42AR3.dat'));
index1 = find(S42A3(:,1) == 100000000);
index2 = find(S42A3(:,1) == 150000000);
index3 = find(S42A3(:,1) == 200000000);
index4 = find(S42A3(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
for m = 1:length(index)
    for n = count2:index(m)
        S42A3(n,2) = S42A3(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = count2+index(m)+1;
end
for l = count2:height(S42A3)
    S42A3(l,2) = S42A3(l,2)*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U3_S42A = (S42A3(:,3) - 0.5*S42A3(:,4)*3*noptotal);
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';

figure(5)
h = plot(S42A1(:,2),U_S42A1,'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')

title('P3-S42A','FontSize',20,'FontWeight','bold')
hold on
plot(S42A2(:,2),U2_S42A,'LineWidth',2)
plot(S42A3(:,2),U3_S42A,'LineWidth',2)
legend('Run 1','Run 2','Run 3','FontSize',18,'Location','Southwest')
hold off

%% Aggregation Start Time:
clc
WT1ts = find(Energy(:,2) <= 5.35, 1, 'last');
count0 = 0;
for i = 1:WT1ts
    if (Energy(i,1)==0)
        count0 = count0+1;
    end
end

count1 = 0;
file = 78;

for i = 1:length(Energy(:,1))
    if (Energy(i,1)==0)
        count1 = count1+1;
    end
    if (count1==file)
        tindex = i;
        ts = Energy(i,2);
    end
end

%% E test potential 




