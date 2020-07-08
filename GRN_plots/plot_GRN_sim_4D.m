%simulate SSA and ODE model
inputs = {1,1,100,100};
[ssa_t, ssa_simdata, ssa_names] = SSA_4D_GRN(inputs);
[ode_t, ode_simdata, ode_names] = ODE_4D_GRN(inputs);

%plots figure
figure;
    colororder([0 0.447058823529412 0.741176470588235;0.850980392156863 0.325490196078431 0.0980392156862745;0.929411764705882 0.694117647058824 0.125490196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.301960784313725 0.745098039215686 0.933333333333333;1 0.411764705882353 0.16078431372549;1 1 0.0666666666666667;0 0 0]);
    plot(ssa_t, ssa_simdata(:,[2,7,12,17]), 'LineWidth', 1.5)
    hold on
    plot(ode_t, ode_simdata(:,[2,7,12,17]), ':', 'LineWidth', 3)
    xlabel('Time (secs)')
    ylabel('Molecule number')
    title('Molecule numbers vs Time')
    legend({'SSA Protein a','SSA Protein b','SSA Protein c','SSA Protein d', ...
        'ODE Protein a', 'ODE Protein b', 'ODE Protein c', 'ODE Protein d'}, ...
        'Location', 'eastoutside')
