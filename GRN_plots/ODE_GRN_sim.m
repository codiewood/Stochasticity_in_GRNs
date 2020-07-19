function [ode_t, ode_simdata, ode_names] = ODE_GRN_sim(Mobj)
%A function to simulate an ODE model of a GRN based on an input model
%object.

%simulate model deterministically
configset = getconfigset(Mobj);
configset.CompileOptions.DimensionalAnalysis = true;
configset.SolverType = 'ode15s';
configset.StopTime = 1000000;
solver = configset.SolverOptions;
set(configset.SolverOptions, 'AbsoluteTolerance', 1.0e-8);

%do sim
[ode_t, ode_simdata, ode_names] = sbiosimulate(Mobj);
end