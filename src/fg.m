function P = fg(Vg, Vl, Pg, Y)
% Power constraint on generators
P = Vg .* (Y(1:2, 1:2) * Vg) + Vg .* (Y(1:2, 3:4) * Vl) - Pg;
end

