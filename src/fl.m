function I = fl(Vg, Vl, Il_bar, Yl_bar, Pl_bar, Y)
% Current constraint on the grid
n_gens = length(Vg);
I = Y(n_gens+1:end, 1:n_gens) * Vg + Y(n_gens+1:end, n_gens+1:end) * Vl + Il_bar + Yl_bar .* Vl; % + Vl.^-1 .* Pl_bar;
end

