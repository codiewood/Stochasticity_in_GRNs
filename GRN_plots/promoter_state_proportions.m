function [prop] = promoter_state_proportions(times,simdata)
v = [1:5];
v(2) = [];
t = length(times);
molA = sum(simdata(1,v));
molB = sum(simdata(1,v+5));
molC = sum(simdata(1,v+10));
molD = sum(simdata(1,v+15));
max_promoters = max([molA,molB,molC,molD]);
prop = zeros(max_promoters+1,4);
c = 0;
for i = [4,9,14,19]
    c = c+1;
    for j = [0:max_promoters]
        count = sum(simdata(:,i) == j);
        prop(j+1,c) = count/t;
    end
end
end