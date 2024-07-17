clear, clc
close all
%h = 250;

%% L38M
clc
L38M1 = table2array(readtable('300billions.dat'));
index1 = find(L38M1(:,1) == 100000000);
index2 = find(L38M1(:,1) == 150000000);
index3 = find(L38M1(:,1) == 200000000);
index4 = find(L38M1(:,1) == 250000000);
index = [index1', index2', index3', index4'];
count2 = 1;
temps= [0.50 0.45 0.40 0.35 0.30 0.28 0.26 0.24 0.22];
setemp = 0.195;
numbeads = 119;
numchain = 24;
noptotal = numbeads*numchain;
for m = 1:length(index)
    for n = count2:index(m)
        L38M1(n,2) = L38M1(n,2)*3.3/sqrt(temps(m)*12);
    end
    count2 = index(m)+1;
end
for l = count2:height(L38M1)
    L38M1(l,2) = (L38M1(l,2)-L38M1(count2,2))*3.3/sqrt(setemp*12);
end

%timeE1 = Energy1(1:cutoff1,1);
%tEavg1 = arrayfun(@(i) mean(timeE1(i:i+h-1)),1:h:height(timeE1)-h+1)';
U_L38M1 = (L38M1(:,3) - 0.5*L38M1(:,4)*3*noptotal)*12.47;
%Uavg1 = arrayfun(@(i) mean(U1(i:i+h-1)),1:h:height(U1)-h+1)';


h = plot(L38M1(count2:end,2),U_L38M1(count2:end),'LineWidth',2);
%ytickformat('%.0f')
ax = ancestor(h, 'axes');
ax.YAxis.Exponent = 0;
set(ax,'FontSize',20)
xlabel('Time, microsecond','FontSize',20,'FontWeight','bold')
ylabel('Potential Energy,kJ/mol','FontSize',20,'FontWeight','bold')

title('P3NExtend-L38M','FontSize',20,'FontWeight','bold')


