%% 
% runs simulation and sets default parameter values
inputs = {1,1,1,1,14,14,14,14,5,14,1,1e-4,1e-2,2,1e-1};

numargs = length(inputs);
args = {1,1,1,1,0,0,0,0,5,14,1,1e-4,1e-2,2,1e-1};
args(1:numargs) = inputs;
[molA,molB,molC,molD,mola,molb,molc,mold,g0,g1,k,hr,fr,ha,fa] = args{:};

Mobj = model_4D_GRN(inputs);
[ode_t, ode_simdata, ode_names] = ODE_GRN_sim(Mobj);

%%
figure;
    for i = [1:4]
        scatter(ode_t, ode_simdata(:,[5*i-3]), '*')
        hold on
    end
    xlim([0 1000])