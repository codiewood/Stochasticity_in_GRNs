%create model
Mobj = sbiomodel('cell');
%add compartment
compObj = addcompartment(Mobj,'comp');
compObj.CapacityUnits = 'liter';
%Reactions
Robj1 = addreaction(Mobj,'A0 -> A0 + a');
Robj2 = addreaction(Mobj,'A1 -> A1 + a');
Robj3 = addreaction(Mobj,'a -> phi');
Robj4 = addreaction(Mobj,'A0 + 2 a <-> A1');
%initialise species amount
Mobj.Species(1).InitialAmount = 1;
Mobj.Species(2).InitialAmount = 0;
Mobj.Species(3).InitialAmount = 0;
Mobj.Species(4).InitialAmount = 0;
%Units of species amounts
for i = 1:4
Mobj.Species(i).InitialAmountUnits = 'molecule';
end
%Mass action kinetics
Kobj1 = addkineticlaw(Robj1,'MassAction');
Kobj2 = addkineticlaw(Robj2,'MassAction');
Kobj3 = addkineticlaw(Robj3,'MassAction');
Kobj4 = addkineticlaw(Robj4,'MassAction');
%Rate parameters
g0 = 14;
g1 = 5;
k=1;
hr = 1e-4;
fr = 1e-2;
    %Rate parameter for Reaction 1
    Pobj1 = addparameter(Kobj1,'g0');
    Pobj1.Value = g0;
    Pobj1.ValueUnits = '1/second';
    Kobj1.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 2
    Pobj2 = addparameter(Kobj2,'g1');
    Pobj2.Value = g1;
    Pobj2.ValueUnits = '1/second';
    Kobj2.ParameterVariableNames = 'g1';
    %Rate parameter for Reaction 3
    Pobj3 = addparameter(Kobj3,'k');
    Pobj3.Value = k;
    Pobj3.ValueUnits = '1/second';
    Kobj3.ParameterVariableNames = 'k';   
    %Rate parameter for Reaction 4
    Pobj4 = addparameter(Kobj4,'hr');
    Pobj4.Value = hr;
    Pobj4.ValueUnits = '1/(molecule*molecule*second)';
    Pobj4r = addparameter(Kobj4,'fr');
    Pobj4r.Value = fr;
    Pobj4r.ValueUnits = '1/second';  
    Kobj4.ParameterVariableNames = {'hr','fr'};  

    
%simulate model stochastically
ssa_configset = getconfigset(Mobj);
ssa_configset.CompileOptions.DimensionalAnalysis = true;
ssa_configset.SolverType = 'ssa';
ssa_configset.StopTime = 1000000;
ssa_solver = ssa_configset.SolverOptions;
ssa_solver.LogDecimation = 10;
%ssa_configset.CompileOptions.DimensionalAnalysis = false;


%do sim
[t, simdata, names] = sbiosimulate(Mobj);
figure;
    histogram(simdata(:,2),15)
figure;
    plot(t, simdata(:,[2]), 'LineWidth', 1.5)
    xlim([0 1000])
    yyaxis right
    plot(t, simdata(:,[1]), 'LineWidth', 1.5)
    ylim([-0.5 2])
    legend('Protein a', 'Unbound promoter')
