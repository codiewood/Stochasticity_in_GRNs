%ode plot and expected protein numbers based on promoter states
figure;
    colororder([0 0.447058823529412 0.741176470588235;0.850980392156863 0.325490196078431 0.0980392156862745;0.929411764705882 0.694117647058824 0.125490196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.301960784313725 0.745098039215686 0.933333333333333;1 0.411764705882353 0.16078431372549;1 1 0.0666666666666667;0 0 0;0.301960784313725 0.745098039215686 0.933333333333333;1 0.411764705882353 0.16078431372549;1 1 0.0666666666666667;0 0 0]);
    title('ODE Model of Molecule numbers vs Time, with expected values based on promoter state')
    plot(ode_t, ode_simdata(:,[2,7,12,17]), '-d', 'LineWidth', 1.5, 'MarkerSize', 5)
    hold on
    xlabel('Time (secs)')
    xlim([0 1000])
    ylabel('Molecule number')
    for i = [1:4]
        g_tot = g1*ode_simdata(:,5*i-1) + g0*(ode_simdata(:,5*i-4) + ode_simdata(:,5*i-2) + ode_simdata(:,5*i));
        tf = g_tot/k;
        plot(ode_t,tf, 'LineStyle', 'None', 'MarkerSize',15,'Marker','.')
    end
    legend({'ODE Protein a', 'ODE Protein b', 'ODE Protein c', 'ODE Protein d', ...
    'Expected protein a', 'Expected protein b', 'Expected protein c', 'Expected protein d'}, ...
    'Location', 'eastoutside')