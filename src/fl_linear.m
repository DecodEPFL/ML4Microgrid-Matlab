function I = fl(Vg, Vl, Il_bar, Yl_bar, Pl_bar, Y)
% Current constraint on the grid without the non-linear contant-power term
n_gens = length(Vg);
I = Y(n_gens+1:end, 1:n_gens) * Vg + Y(n_gens+1:end, n_gens+1:end) * Vl + Il_bar + Yl_bar .* Vl;
end

