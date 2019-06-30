close all;
clear all;

er = 11.44;
mu = 4 * pi * 10^(-7);
epsilon = 8.85 * 10^(-12);


S = 20e-6;
Wline = [36e-6 44e-6 52e-6];

k = S./(Wline);

Rs = 3.2;
kb = 8.62 .* 10^(-5);
Tc = 14.7;
del = 1.76 .* kb .* Tc;

Lks = (6.582 .* 10^(-16)) .* Rs ./ (pi .* del);

gc = 5.61 .* 10^4;
gg = 2.79 .* 10^4;


Lkc = gc .* Lks;
Lkg = gg .* Lks;

L_resonator = 400e-6;

Lkl = (Lkc + Lkg);

Lgl = (mu/4) .* (ellipke(sqrt(1 - k.^2)))./ellipke(k);
Cl = 4 * epsilon .* (1 + er)/2 .* ellipke(k) ./ (ellipke(sqrt(1 - k.^2)));

Z = sqrt((Lgl+Lkl)./Cl);


t = 300e-9;
A = S * t;
L = 400e-6;

Resistivity = Rs * t