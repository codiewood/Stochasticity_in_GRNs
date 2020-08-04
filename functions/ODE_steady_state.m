function [final_states] = ODE_steady_state(simdata)
dims = size(simdata);
time = dims(1);
final_states = zeros(1,4);
i = 1;
for j = [2,7,12,17]
    val = simdata(time,j);
    final_states(i) = val;
    i = i+1;
end

