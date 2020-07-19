function [ssa_t, ssa_simdata, ssa_names] = SSA_GRN_sim(Mobj)
%A function to simulate a stochastic model of a GRN based on an input model
%object.
    
%simulate model stochastically
configset = getconfigset(Mobj);
configset.CompileOptions.DimensionalAnalysis = true;
configset.SolverType = 'ssa';
configset.StopTime = 1000000;
solver = configset.SolverOptions;
solver.LogDecimation = 100;

%do sim
[ssa_t, ssa_simdata, ssa_names] = sbiosimulate(Mobj);
end

