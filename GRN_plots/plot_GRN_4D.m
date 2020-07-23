%   inputs = {molA,molB,molC,molD,mola,molb,molc,mold,g0,g1,k,hr,fr,ha,fa}
%       molA, molB, molC, molD; Number of molecules of unbound promoters A - D
%       mola, molb, molc, mold; Number of molecules of proteins a - d
%       g0; Rate parameter of unactivated protein production
%       g1; Rate parameter of activated protein production
%       k; Protein degradation rate
%       hr; Rate parameter of repressor binding
%       fr; Rate parameter of repressor unbinding
%       ha; Rate parameter of activator binding
%       fa; Rate parameter of activator unbinding
%% 
% runs simulations and sets default parameter values
inputs = {1,1,1,1,0,0,0,0,5,14,1,1e-4,1e-2,2,1e-1};

numargs = length(inputs);
args = {1,1,1,1,0,0,0,0,5,14,1,1e-4,1e-2,2,1e-1};
args(1:numargs) = inputs;
[molA,molB,molC,molD,mola,molb,molc,mold,g0,g1,k,hr,fr,ha,fa] = args{:};

Mobj = model_4D_GRN(inputs);
[ssa_t, ssa_simdata, ssa_names] = SSA_GRN_sim(Mobj);
[ode_t, ode_simdata, ode_names] = ODE_GRN_sim(Mobj);
%% 
%plots figures
figure;
    colororder([0 0.447058823529412 0.741176470588235;0.850980392156863 0.325490196078431 0.0980392156862745;0.929411764705882 0.694117647058824 0.125490196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.301960784313725 0.745098039215686 0.933333333333333;1 0.411764705882353 0.16078431372549;1 1 0.0666666666666667;0 0 0;0.301960784313725 0.745098039215686 0.933333333333333;1 0.411764705882353 0.16078431372549;1 1 0.0666666666666667;0 0 0]);
    title('Molecule numbers vs Time, with expected values based on promoter state')
    plot(ssa_t, ssa_simdata(:,[2,7,12,17]), 'LineWidth', 1.5)
    hold on
    plot(ode_t, ode_simdata(:,[2,7,12,17]), '-d', 'LineWidth', 1.5, 'MarkerSize', 5)
    xlabel('Time (secs)')
    xlim([0 1000])
    ylabel('Molecule number')
    for i = [1:4]
        g_tot = g1*ssa_simdata(:,[5*i-1]) + g0*(ssa_simdata(:,[5*i-4]) + ssa_simdata(:,[5*i-2]) + ssa_simdata(:,[5*i]));
        tf = g_tot/k;
        plot(ssa_t,tf, 'LineStyle', 'None', 'MarkerSize',15,'Marker','.')
    end
    legend({'SSA Protein a','SSA Protein b','SSA Protein c','SSA Protein d', ...
    'ODE Protein a', 'ODE Protein b', 'ODE Protein c', 'ODE Protein d', ...
    'Expected protein a', 'Expected protein b', 'Expected protein c', 'Expected protein d'}, ...
    'Location', 'eastoutside')
    
figure;
    t = tiledlayout(2,2);
    nexttile;
    a = histogram(ssa_simdata(:,2),15);
    title('Protein a');
    a.FaceColor = [0 0.447058823529412 0.741176470588235];
    nexttile;
    b = histogram(ssa_simdata(:,7),15);
    title('Protein b');
    b.FaceColor = [0.850980392156863 0.325490196078431 0.0980392156862745];
    nexttile;
    c = histogram(ssa_simdata(:,12),15);
    title('Protein c');
    c.FaceColor = [0.929411764705882 0.694117647058824 0.125490196078431];
    nexttile;
    d = histogram(ssa_simdata(:,17),15);
    title('Protein d')
    d.FaceColor = [0.494117647058824 0.184313725490196 0.556862745098039];
    xlabel(t,'Molecule number')
    ylabel(t,'Frequency')
%%
%calculate proportion of time promoters spend in each state

t = length(ssa_t);
v = zeros(1,4);
j = 1;
for i = [4,9,14,19];
    s = sum(ssa_simdata(:,[i]));
    p = s/t;
    v(j) = p;
    j = j+1;
end
v
