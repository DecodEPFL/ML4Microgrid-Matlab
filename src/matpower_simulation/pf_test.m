clear all
close all
clc
rng(11);

define_constants;
mpc = loadcase('case5_modified');

mpopt_doe = mpoption;
mpopt_doe = mpoption(mpopt_doe, 'fmincon.alg', 6, 'opf.ac.solver', 'FMINCON', 'opf.start', 3);

r = runpf(mpc, mpopt_doe);