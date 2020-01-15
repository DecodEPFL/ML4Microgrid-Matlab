function [c, ceq] = fl_constraint(Vg, Vl, I_bar, Y_bar, P_bar, Y)
%fl_constraint Power flow equality constraint
c = [];
ceq = fl(Vg, Vl, I_bar, Y_bar, P_bar, Y);
end

