function Vl = load_voltages_from_generator_voltages(Vg, I_bar, Y_bar, P_bar, Y)
%load_voltages_from_generator_voltages Get load voltages from generator
%voltages by inverting DC power flow equations
Vl_var = sym('Vl_var', [2, 1]);
Vl_symb = vpasolve(fl(Vg, Vl_var, I_bar, Y_bar, P_bar, Y) == 0, Vl_var);
[~, i_max] = max(Vl_symb.Vl_var1 + Vl_symb.Vl_var2);
Vl = [Vl_symb.Vl_var1(i_max), Vl_symb.Vl_var2(i_max)]';
end

