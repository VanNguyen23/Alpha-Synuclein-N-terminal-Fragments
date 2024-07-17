clear, clc
bpn = table2array(readtable('lastbpntr.dat'));
chains = 24;
parallel = 0;
antiparallel = 0;
numbeads1 = 119;
numbeads2 = 83;
chnln1 = 31;
chnln2 = 21;
count(1:(119*12+83*12)) = 0;
noptotal = numbeads1*chains/2+numbeads2*chains/2;
for j = 1:length(bpn)
    for i = 1:(119*12+83*12)
        if bpn(j,2) == i
            count(i) = count(i)+1;
        end
    end
end

plot([1:(119*12+83*12)],count)

            
