% Parameters
W = 1600*9.81;       % Weight (N)
S = 40;              % Wing area (m^2)
rho0 = 1.225;        % Sea-level density (kg/m^3)
g = 9.81;
Rng = 7000e3;        % Range (m)
cost_budget = 184;   % Dollars
gal_cost = 7;        % $/gallon
Egal = 1.20106e8;    % J/gal (Jet A)
eta = 0.85;          % Propulsive efficiency
C_D0 = 0.02;         % Parasite drag coeff
e = 0.9;             % Span efficiency (elliptical ~1)

% Search over AR values
AR_vals = 10:1:100;
for AR = AR_vals
    min_cost = inf;
    best = [];
    for h = linspace(0,1000,101)       % altitude sweep
        % simple ISA troposphere (lapse 6.5K/km)
        T = 288.15 - 0.0065*h;
        p = 101325*(T/288.15)^(g/(0.0065*287));
        rho = p/(287*T);
        a = sqrt(1.4*287*T);
        % optimum lift for min drag
        k = 1/(pi*e*AR);
        C_L_opt = sqrt(C_D0/k);
        V_opt = sqrt(2*W/(rho*S*C_L_opt));
        % enforce Mach limit
        if V_opt > 0.3*a
            V = 0.3*a;
        else
            V = V_opt;
        end
        % Compute CL at chosen V
        C_L = W/(0.5*rho*V^2*S);
        % Compute drag
        C_D = C_D0 + k*C_L^2;
        D = 0.5*rho*V^2*S*C_D;
        % Cost = fuel energy / (Egal) * $7
        cost = D * Rng * gal_cost / (eta * Egal);
        if cost < min_cost
            min_cost = cost;
            best = [AR, h, V, C_L, D, cost];
        end
    end
    if best(6) <= cost_budget
        fprintf('AR=%.1f: cost=%.1f, alt=%.0fm, V=%.1fm/s\n', best(1), best(6), best(2), best(3));
        break;
    end
end
% (Output shows optimal AR and conditions)

