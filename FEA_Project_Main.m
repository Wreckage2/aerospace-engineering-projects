%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D FE Linear Triangles - Complete Final Script    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

%%%%%%%%%%%%%
% INPUT FILES
%%%%%%%%%%%%%
globalnode = load('NodeFile.txt');
connect    = load('ConnectFile.txt');

M = size(globalnode,1);   % number of nodes
E = size(connect,1);      % number of elements

x = globalnode(:,2);
y = globalnode(:,3);

%%%%%%%%%%%%%
% MATERIALS
%%%%%%%%%%%%%
b  = 0.1;       % thickness
YM = 20;        % Young's modulus
PR = 1/4;       % Poisson's ratio

Dp = YM/(1-PR^2) * [1 PR 0; PR 1 0; 0 0 (1-PR)];

%%%%%%%%%%%%%
% BCs
%%%%%%%%%%%%%
% From assignment:
%   bottom-center fixed in both u & v
%   bottom-left and bottom-right have zero v only

uzerodisp = [20];         % zero u
vzerodisp = [19 20 21];   % zero v

%%%%%%%%%%%%%
% LOAD CASE
%%%%%%%%%%%%%

% LOADCASE = 1 → uniform traction (Part I)
% LOADCASE = 2 → point load (Part II)

LOADCASE = 2;   % <-- set this for whichever part you are running

% Initialize global force vector
F = zeros(2*M,1);

if LOADCASE == 1
    % ------------------- Part I: Uniform traction ----------------------
    p0 = 2;
    f0 = p0/(3*2+2);

    F(2)  = -f0;
    F(4)  = -2*f0;
    F(6)  = -2*f0;
    F(8)  = -2*f0;
    F(10) = -f0;

elseif LOADCASE == 2
    % ------------------- Part II: Point load ---------------------------
    p0 = 0.4;
    F(6) = -p0;   % vertical DOF of node #3
end

%%%%%%%%%%%%%
% GLOBAL STIFFNESS
%%%%%%%%%%%%%
Kg = zeros(2*M, 2*M);

for p = 1:E
    pconnect = connect(p,2:4);

    xp = x(pconnect);
    yp = y(pconnect);

    % B matrix and area
    [Bp,Area] = CalcBmatTri(xp,yp);

    % Element stiffness
    Kp = Bp' * Dp * Bp * Area * b;

    % Assembly into Kg
    for ip = 1:3
        ig = pconnect(ip);

        for jp = 1:3
            jg = pconnect(jp);

            Kg(2*ig-1, 2*jg-1) = Kg(2*ig-1, 2*jg-1) + Kp(2*ip-1, 2*jp-1);
            Kg(2*ig,   2*jg-1) = Kg(2*ig,   2*jg-1) + Kp(2*ip,   2*jp-1);
            Kg(2*ig-1, 2*jg)   = Kg(2*ig-1, 2*jg)   + Kp(2*ip-1, 2*jp);
            Kg(2*ig,   2*jg)   = Kg(2*ig,   2*jg)   + Kp(2*ip,   2*jp);
        end
    end
end

%%%%%%%%%%%%%
% APPLY BCs
%%%%%%%%%%%%%
for i = 1:length(uzerodisp)
    a = uzerodisp(i);
    Kg(2*a-1,:) = 0;     Kg(:,2*a-1) = 0;
    Kg(2*a-1,2*a-1) = 1;
    F(2*a-1) = 0;
end

for i = 1:length(vzerodisp)
    a = vzerodisp(i);
    Kg(2*a,:) = 0;       Kg(:,2*a) = 0;
    Kg(2*a,2*a) = 1;
    F(2*a) = 0;
end

%%%%%%%%%%%%%
% SOLVE SYSTEM
%%%%%%%%%%%%%
Ug = Kg \ F;

u = Ug(1:2:end);
v = Ug(2:2:end);

%%%%%%%%%%%%%
% PLOT DISPLACEMENTS
%%%%%%%%%%%%%
IPs = globalnode(:,2:3);
DispPs = [IPs(:,1)+u , IPs(:,2)+v];

figure(1)
scatter(IPs(:,1), IPs(:,2), 'r', 'filled'); hold on
scatter(DispPs(:,1), DispPs(:,2), 'k'); hold off
axis equal
title('Original vs. Displaced Node Locations');
xlabel('x'); ylabel('y');

%%%%%%%%%%%%%
% STRESS COMPUTATION
%%%%%%%%%%%%%
sigY = zeros(E,1);
xp_elem = zeros(3,E);
yp_elem = zeros(3,E);

for p = 1:E
    pconnect = connect(p,2:4);

    xp = x(pconnect);
    yp = y(pconnect);

    ue = [Ug(2*pconnect(1)-1); Ug(2*pconnect(1));
          Ug(2*pconnect(2)-1); Ug(2*pconnect(2));
          Ug(2*pconnect(3)-1); Ug(2*pconnect(3))];

    xp_elem(:,p) = xp + ue(1:2:5);
    yp_elem(:,p) = yp + ue(2:2:6);

    [Bp,Area] = CalcBmatTri(xp,yp);

    eps = Bp*ue;
    sig = Dp*eps;

    sigY(p) = sig(2);

end

%%%%%%%%%%%%%
% STRESS CONTOUR PLOT
%%%%%%%%%%%%%
figure(2)
fill(xp_elem, yp_elem, sigY)
colorbar; axis equal
title('\sigma_y Stress Distribution')
xlabel('x'); ylabel('y');

