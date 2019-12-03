function P = fg(Vg, Vl, Pg, Y)
% Power constraint on generators
n_gens = length(Vg);
P = Vg .* (Y(1:n_gens, 1:n_gens) * Vg) + Vg .* (Y(1:n_gens, n_gens+1:end) * Vl) - Pg;
end

