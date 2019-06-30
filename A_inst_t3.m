clear all;
close all;

mu = 4 * pi * 10^(-7);
epsilon = 8.85 * 10^(-12);

gc = 6.71 .* 10^4;
gg = 2.79 .* 10^4;

S = 20e-6;
W = 8e-6;
k = (S)/(S + 2*W);



Lks = (6.582 .* 10^(-16) * 0.2788)/(pi .* 1.76 .* 8.62 .* 10^(-5) * 1.28);

er = 11.44;

Cl = 4 * epsilon .* (1 + er)/2 .* ellipke(k) ./ (ellipke(sqrt(1 - k.^2)));
Lg = (mu/4) .* (ellipke(sqrt(1 - k.^2)))./ellipke(k);

Lkc = Lks .* gc;
Lkg = Lks .* gg;

l = [160e-6 260e-6];

L_resonator = 5e-3;


Lalu = [1e-3 1e-3];
Lnbtin = L_resonator - Lalu;

Lk1 = (Lkc + Lkg) .* Lnbtin;
Lk2 = (Lkc + Lkg) .* Lalu;
Lk = Lk1 + Lk2;



C = Cl .* (Lnbtin + Lalu);
Lgl = Lg .* (Lnbtin + Lalu);

F0 = 1./(4 .* (sqrt((Lgl + Lk) .* C)))
d = 300e-9;

gc = (1/(4 .* S .* (1 - k.^2).* (ellipke(k)).^2 )) .* (pi + log(4 .* pi .* S ./ d) - k .* log((1 + k)./(1 - k)));

gg = (1/(4 .* S .* (1 - k.^2).* (ellipke(k)).^2 )) .* (pi + log(4 .* pi .* (S+2.*W) ./ d) - (1/k) .* log((1 + k)./(1 - k)));

f0 = 5.844005e9;
f1 = 5.844155e9;

S21min = 10^(-20.664260/20);
S21dB = 10 .* log10((1 + abs(S21min).^2)/2)

omega0 = 2 * pi * f0;
domega = 2 * pi * abs(f1 - f0);
Q = omega0/(2*domega);

Qc = Q/(1 - S21min)


% omega0 = 2 * pi * 5.842612e9;
% domega = 2 * pi * (5.842678e9 - 5.842612e9);
% 
%s21 = 10^(-14/20);

% Z = sqrt((Lgl + Lk)/C)
% alpha = Lk/(Lk + Lgl)

%Q = (omega0/(2 * domega)) * sqrt((1/abs(s21)^2 - 1))








