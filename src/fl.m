function I = fl(Vg, Vl, Il_bar, Yl_bar, Pl_bar, Y)
% Current constraint on the grid
I = Y(1:2, 3:4) * Vg + Y(3:4, 3:4) * Vl + Il_bar + Yl_bar .* Vl; % + Vl.^-1 .* Pl_bar;
end

