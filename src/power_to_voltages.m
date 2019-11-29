function [Vl, Vg, Pg] = power_to_voltages(Y, Il_bar, Yl_bar, Pl_bar, Pg_bar, x0)
%UNTITLED Solve optimization problem for power to voltage translation
prob = optimproblem();
Vg = optimvar('Vg', 2, 'LowerBound', 0);
Vl = optimvar('Vl', 2, 'LowerBound', 0);
Pg = optimvar('Pg', 2, 'LowerBound', 0);
prob.Objective = sum((Pg - Pg_bar).^2);
prob.Constraints.cons1 = fg(Vg, Vl, Pg, Y) == 0;
prob.Constraints.cons2 = fl(Vg, Vl, Il_bar, Yl_bar, Pl_bar, Y) == 0;
sol = solve(prob, x0);
Vl = sol.Vl;
Vg = sol.Vg;
Pg = sol.Pg;
end

