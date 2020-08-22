% analytical derivation
clearvars

syms rho positive real;
syms h positive real;
syms hT positive real;
syms rhop rhoT real;
syms hp real;
syms u real;
syms c real positive;
syms c1 real positive;


Gamma = [rhop 0 rhoT; ...
    rhop * u rho rhoT*u; ...
    rhop * (h + 0.5*u^2) + rho * hp - 1 rho*u rhoT*(h+0.5*u^2)+rho*hT];

A = [rhop * u rho rhoT*u;...
    rhop * u^2+1 2*rho*u rhoT*u^2;...
    (rhop*(h+0.5*u^2)+rho*hp)*u rho*(h + 1.5*u^2) u*(rhoT*(h+0.5*u^2)+rho*hT)];

% assume(c1 == rhoT + hT*rho*rhop - hp*rho*rhoT);
% assume(c == ((hT*rho)/(rhoT + hT*rho*rhop - hp*rho*rhoT))^(1/2));

% f=(hT*rho*(rhoT + hT*rho*rhop - hp*rho*rhoT))^(1/2)/(rhoT + hT*rho*rhop - hp*rho*rhoT);
% f=(hT*rho*c1)^(1/2)/c1;
% simplify(f,'All',true)

% assume(c == ((hT*rho)/(rhoT + hT*rho*rhop - hp*rho*rhoT))^(1/2));

[V,D] = eig(Gamma^(-1) * A);
assume(c1 == rhoT + hT*rho*rhop - hp*rho*rhoT);

% assume(c == ((hT*rho)/(rhoT + hT*rho*rhop - hp*rho*rhoT))^(1/2));
V = collect(simplify(V));
D = collect(simplify(D));

assume(c == (c1*hT*rho)^(1/2)/c1);
V = collect(simplify(V))
D = collect(simplify(D))

% e = eig(Gamma^(-1) * A);
% simplify(e)