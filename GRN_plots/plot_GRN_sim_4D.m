GRN_sim_4D.m

figure;
plot(t, simdata(:,[2,7,12,17]), 'Linewidth',1.5)
xlabel('Time')
ylabel('Molecule number')
title('Number vs Time')
legend('Protein a','Protein b','Protein c','Protein d')