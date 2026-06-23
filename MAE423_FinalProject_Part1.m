%  MAE 423 FINAL PROJECT - PART I: ROCKET DESIGN
%       ASSUMPTIONS:
%      - The rocket flies through atmospheric air at sea level
%      - At the nozzle exit, the exhaust is ideally expanded to match
%        ambient conditions:
%      - The flow from chamber to exit is isentropic.
%      - Thin-wall pressure vessel: delta = D*P_o/(2*sigma_s)

clear; clc; close all;

%  DESIGN INPUT PARAMETERS  
M          = 10;        % Payload mass [kg]                    (5, 10, 20)
h_max_ft   = 20000;     % Target maximum altitude [ft]    (10k, 20k, 30k)
a_star_max = 10;        % Normalized max acceleration         (5, 10, 20)
SM         = 2;         % Static margin [calibers]
rho_s      = 1820;      % Shell density (carbon-epoxy) [kg/m^3]
rho_p      = 1772;      % Propellant density [kg/m^3]
sigma_s    = 90e6;      % Shell working stress [Pa]
N_fins     = 3;         % Number of fins

% Physical constants
g0  = 9.81;             % Gravitational acceleration [m/s^2]
P_a = 101325;           % Atmospheric pressure (sea level) [Pa]

% Atmospheric properties 
gamma        = 1.4;     % Specific heat ratio of air
T_e          = 298;     % Exit static temperature [K]  (= ambient)
gas_constant = 287;     % Specific gas constant for air [J/(kg*K)]
P_e     = P_a;          % Exit pressure (ideal expansion to atmosphere) [Pa]
c_p     = gamma * gas_constant / (gamma - 1);   % c_p of air [J/(kg*K)]

% Convert altitude to SI
h_max = h_max_ft * 0.3048;  % Target apogee [m]

%TRAJECTORY  
R_opt = a_star_max + 1;
fR    = log(R_opt)/2 * ( log(R_opt) - 2 ) + (R_opt - 1)/R_opt;
W_eq  = sqrt( h_max * g0 / fR );
t_b   = (R_opt - 1) * W_eq / ( g0 * R_opt );

% Burnout state (for reporting / verification)
v_b      = W_eq*log(R_opt) - g0*t_b;
h_b      = W_eq*t_b*(1 - log(R_opt)/(R_opt-1)) - 0.5*g0*t_b^2;
h_apogee = h_b + v_b^2/(2*g0);

%CHAMBER PRESSURE  (isentropic, via Mach number)
a_e = sqrt(gamma * gas_constant * T_e);
M_e = W_eq / a_e;
T_o = T_e * (1 + (gamma-1)/2 * M_e^2);
P_o = P_e * (1 + (gamma-1)/2 * M_e^2)^(gamma/(gamma-1));

%GEOMETRIC OPTIMIZATION  
D_grid = linspace(0.04, 0.50, 300);
nD     = length(D_grid);

% Storage
res_D   = NaN(nD,1);  res_L   = NaN(nD,1);
res_Ms  = NaN(nD,1);  res_Mp  = NaN(nD,1);  res_Mo  = NaN(nD,1);
res_Lp  = NaN(nD,1);  res_lam = NaN(nD,1);
res_XCG = NaN(nD,1);  res_XCP = NaN(nD,1);
res_del = NaN(nD,1);

for i = 1:nD
    D     = D_grid(i);
    delta = D * P_o / (2 * sigma_s);     % Thin-wall hoop-stress thickness

    stabResid = @(L) stability_residual(L, D, delta, R_opt, M, ...
                                         rho_s, rho_p, N_fins, SM);

    if stabResid(0.5*D) * stabResid(100*D) > 0
        continue;
    end
    L = fzero(stabResid, [0.5*D, 100*D]);
    if isnan(L) || L <= 0,  continue;  end

    [X_CG, X_CP, M_s, M_p, L_p] = compute_geometry( ...
        L, D, delta, R_opt, M, rho_s, rho_p, N_fins);

    % Feasibility checks
    if L_p > (L + D),       continue;  end
    if M_s < 0 || M_p < 0,  continue;  end

    M_o    = M + M_s + M_p;
    lambda = M / M_o;

    res_D(i)   = D;     res_L(i)   = L;
    res_Ms(i)  = M_s;   res_Mp(i)  = M_p;   res_Mo(i)  = M_o;
    res_Lp(i)  = L_p;   res_lam(i) = lambda;
    res_XCG(i) = X_CG;  res_XCP(i) = X_CP;
    res_del(i) = delta;
end

% optimum
valid = ~isnan(res_lam);
if ~any(valid)
    error('No feasible rocket design found. Check input parameters.');
end
[lambda_opt, idxOpt] = max(res_lam);

D_opt     = res_D(idxOpt);
L_opt     = res_L(idxOpt);
M_s_opt   = res_Ms(idxOpt);
M_p_opt   = res_Mp(idxOpt);
M_o_opt   = res_Mo(idxOpt);
L_p_opt   = res_Lp(idxOpt);
X_CG_opt  = res_XCG(idxOpt);
X_CP_opt  = res_XCP(idxOpt);
delta_opt = res_del(idxOpt);
epsilon   = M_s_opt / (M_s_opt + M_p_opt);

%   THRUST  
m_dot   = M_p_opt / t_b;
F       = m_dot * W_eq;
M_b     = M_o_opt / R_opt;
TW_0    = F / (M_o_opt * g0);
TW_b    = F / (M_b     * g0);
F_lbf   = F / 4.4482;

%   OUTPUTS
fprintf('  MAE 423 Final Project - Part I : Rocket Design Results\n');

fprintf('\n--- Inputs ---\n');
fprintf('  Payload mass        M          = %7.2f kg\n', M);
fprintf('  Target altitude     h_max      = %7.0f ft  (%6.1f m)\n', h_max_ft, h_max);
fprintf('  Norm. max accel     a*_max     = %7.1f\n', a_star_max);
fprintf('  Static margin       SM         = %7.1f calibers\n', SM);
fprintf('  Shell density       rho_s      = %7.0f kg/m^3\n', rho_s);
fprintf('  Propellant density  rho_p      = %7.0f kg/m^3\n', rho_p);
fprintf('  Shell stress        sigma_s    = %7.0f MPa\n', sigma_s/1e6);
fprintf('  Number of fins      N          = %7d\n', N_fins);

fprintf('\n--- Atmosphere & Exit Conditions (ideal expansion) ---\n');
fprintf('  gamma                          = %7.2f\n', gamma);
fprintf('  gas constant                   = %7.0f J/(kg*K)\n', gas_constant);
fprintf('  Exit temperature    T_e        = %7.1f K\n', T_e);
fprintf('  Exit pressure       p_e        = %7.3f kPa\n', P_e/1e3);

fprintf('\n--- Trajectory ---\n');
fprintf('  Optimal mass ratio  R_opt      = %7.4f  (= a*_max + 1)\n', R_opt);
fprintf('  Exhaust velocity    W_eq       = %7.1f m/s\n', W_eq);
fprintf('  Equivalent Isp                 = %7.1f s\n', W_eq/g0);
fprintf('  Burn time           t_b        = %7.2f s\n', t_b);
fprintf('  Burnout velocity    v_b        = %7.1f m/s\n', v_b);
fprintf('  Burnout altitude    h_b        = %7.1f m\n', h_b);
fprintf('  Predicted apogee    h_max(pred)= %7.1f m   (%.0f ft)\n', ...
        h_apogee, h_apogee/0.3048);

fprintf('\n--- Propulsion ---\n');
fprintf('  Speed of sound (exit) a_e      = %7.1f m/s\n', a_e);
fprintf('  Exit Mach number    M_e        = %7.3f\n', M_e);
fprintf('  Stagnation temp     T_o        = %7.1f K\n', T_o);
fprintf('  Chamber pressure    p_o        = %7.4f MPa\n', P_o/1e6);
fprintf('  Pressure ratio      p_o/p_a    = %7.3f\n', P_o/P_a);

fprintf('\n--- Geometry ---\n');
fprintf('  Body diameter       D          = %7.4f m   (%5.1f mm)\n', D_opt, D_opt*1000);
fprintf('  Body length         L          = %7.4f m\n', L_opt);
fprintf('  Total length        L_tot=2D+L = %7.4f m\n', 2*D_opt + L_opt);
fprintf('  Aspect ratio        L/D        = %7.2f\n', L_opt/D_opt);
fprintf('  Shell thickness     delta      = %8.6f m   (%5.3f mm)\n', delta_opt, delta_opt*1000);
fprintf('  delta / D                      = %8.6f\n', delta_opt/D_opt);
fprintf('  Propellant length   L_p        = %7.4f m\n', L_p_opt);
fprintf('  Center of pressure  X_CP       = %7.4f m\n', X_CP_opt);
fprintf('  Center of gravity   X_CG       = %7.4f m\n', X_CG_opt);
fprintf('  Stability margin    X_CP-X_CG  = %7.4f m  (%.2f cal.)\n', ...
        X_CP_opt-X_CG_opt, (X_CP_opt-X_CG_opt)/D_opt);

fprintf('\n--- Masses ---\n');
fprintf('  Structural mass     M_s        = %7.3f kg\n', M_s_opt);
fprintf('  Propellant mass     M_p        = %7.3f kg\n', M_p_opt);
fprintf('  Initial mass        M_o        = %7.3f kg\n', M_o_opt);
fprintf('  Payload ratio       lambda     = %7.4f\n', lambda_opt);
fprintf('  Structural ratio    epsilon    = %7.4f\n', epsilon);

fprintf('\n--- Thrust ---\n');
fprintf('  Mass flow rate      m_dot      = %7.4f kg/s\n', m_dot);
fprintf('  Thrust              F          = %7.1f N   (%6.1f lbf)\n', F, F_lbf);
fprintf('  Lift-off T/W        F/(M_o*g0) = %7.3f\n', TW_0);
fprintf('  Burn-out T/W        F/(M_b*g0) = %7.3f   (= a*_max + 1)\n', TW_b);

%   LOCAL FUNCTIONS

function res = stability_residual(L, D, delta, R, M, rho_s, rho_p, N_fins, SM)
%   STABILITY_RESIDUAL  Zero when the static-margin requirement is satisfied:
%   res = (X_CP - X_CG) - SM*D
    [X_CG, X_CP, ~, ~, ~] = compute_geometry( ...
        L, D, delta, R, M, rho_s, rho_p, N_fins);
    res = (X_CP - X_CG) - SM*D;
end


function [X_CG, X_CP, M_s, M_p, L_p] = compute_geometry( ...
            L, D, delta, R, M, rho_s, rho_p, N_fins)

%   MASSES
% Nose cone (thin conical shell)
%   r = D/2,  l_slant = sqrt((D/2)^2 + D^2) = D*sqrt(5)/2
A_nose = pi * D^2 * sqrt(5) / 4;
M_N    = rho_s * delta * A_nose;

% Body shell (thin cylinder, length L)
M_B  = rho_s * delta * pi * D * L;

% Fin-section shell (thin cylinder, length D)
M_FS = rho_s * delta * pi * D * D;

% Fins: N triangular thin plates, area each = (1/2)*c_r*s = D^2/2
A_fin   = 0.5 * D * D;
M_fins  = N_fins * rho_s * delta * A_fin;

% Structural mass
M_s = M_N + M_B + M_FS + M_fins;

% Propellant mass from R = M_o / M_b,  M_b = M + M_s
    M_p = (R - 1) * (M + M_s);

% Propellant column length (cylinder of diameter D)
    L_p = M_p / (rho_p * pi * D^2 / 4);

%   CENTER-OF-GRAVITY 
X_payload = D;                          % payload at base of nose
X_N_cg    = (2/3) * D;                  % cone shell centroid
X_B_cg    = D + L/2;                    % body cylinder centroid
X_FS_cg   = D + L + D/2;                % fin section cylinder centroid
X_fins_cg = D + L + (2/3)*D;            % right-triangle fin centroid
X_p_cg    = (2*D + L) - L_p/2;          % propellant centroid (filled
                                            % from the bottom upward)

M_total = M + M_s + M_p;
X_CG = ( M*X_payload + M_N*X_N_cg + M_B*X_B_cg + M_FS*X_FS_cg ...
    + M_fins*X_fins_cg + M_p*X_p_cg ) / M_total;

% CENTER-OF-PRESSURE 
% Nose cone:
%   (CN_alpha)_N = 2,  X_N_cp = (2/3)*L_nose
CNa_N  = 2;
X_N_cp = (2/3) * D;

% Fins 
%   c_r = D,  c_t = 0,  s = D,  m = D  (sweep length)
%     l = sqrt( (m + c_t/2 - c_r/2)^2 + s^2 )
c_r = D;
c_t = 0;
s   = D;
m   = D;
l_mid = sqrt( (m + c_t/2 - c_r/2)^2 + s^2 );

% Fin normal-force coefficient slope
%   (CN_alpha)_f = 4*N*(s/d)^2 / [1 + sqrt(1 + (2*l/(c_r+c_t))^2)]
CNa_f = 4*N_fins*(s/D)^2 / ( 1 + sqrt(1 + (2*l_mid/(c_r+c_t))^2) );

% Fin-body interference factor:
%   K_fb = 1 + R_body / (s + R_body),    R_body = D/2
R_body = D/2;
K_fb = 1 + R_body / (s + R_body);
CNa_fb = K_fb * CNa_f;

% Fin CP location:  X_f_cp = X_root_LE + delta_x_f
%   delta_x_f = (m/3)*(c_r + 2*c_t)/(c_r+c_t)
%               + (1/6)*( c_r + c_t - c_r*c_t/(c_r+c_t) )
X_root_LE = D + L;     % top of fin section (where root LE meets body)
    if (c_r + c_t) > 0
        dx_f = (m/3)*(c_r + 2*c_t)/(c_r + c_t) ...
             + (1/6)*( c_r + c_t - c_r*c_t/(c_r + c_t) );
    else
        dx_f = m/3 + c_r/6;     % limit for c_t -> 0
    end
    X_f_cp = X_root_LE + dx_f;

    % Combined center of pressure (moment balance about nose tip):
    %   X_CP = sum(CNa_i * X_i) / sum(CNa_i)
    X_CP = ( CNa_N*X_N_cp + CNa_fb*X_f_cp ) / ( CNa_N + CNa_fb );
end