%simulate SSA and ODE model
inputs = {1,1,1,1};
[ssa_t, ssa_simdata, ssa_names] = SSA_4D_GRN(inputs);
[ode_t, ode_simdata, ode_names] = ODE_4D_GRN(inputs);

%plots figure
figure;
    plot(ssa_t, ssa_simdata(:,[2,7,12,17]), '-' , 'LineWidth', 1.5)
    hold on
    plot(ode_t, ode_simdata(:,[2,7,12,17]), ':', 'LineWidth', 3)
    xlabel('Time (secs)')
    ylabel('Molecule number')
    title('Molecule numbers vs Time')
    legend('SSA Protein a','SSA Protein b','SSA Protein c','SSA Protein d', ...
        'ODE Protein a', 'ODE Protein b', 'ODE Protein c', 'ODE Protein d')
