%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%   Finite Element Code for Two-node Bar Elements                        %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clc
clear
format short
tic

fprintf('\nFE_Bar_2026_Project.m\n')

% inpflg = 0;  % Automatically generate nodes and elements along x-axis
inpflg = 1;  % Read nodes and elements from input files

% iplot = 0;  % No plotting
% iplot = 1;  % Plot geometry and matrices
iplot = 2;  % Plot geometry, matrices and deformed shape

if inpflg == 0
    L = 1;
    M = 128;
    m = M+1;
    dx = L/M;
    [Nodes,N] = Generate_Nodes(m,dx);
    [Elems,E] = Generate_Elements(m,M);
else
    [Nodes,N] = Read_Nodes;
    [Elems,E] = Read_Elements;
end

Attached = zeros(N,1);
for e=1:E
    Attached(Elems(e,4)) = Attached(Elems(e,4))+1;
    Attached(Elems(e,5)) = Attached(Elems(e,5))+1;
end
minAttached = min(Attached);
if minAttached < 1
    fprintf('\nWarning: The following Nodes are not attached to Elements\n')
    for i=1:N
        if Attached(i) == 0
            fprintf('\n %d ',i)
        end
    end
    fprintf('\n\n')
end

fprintf('\nNodes\n')
fprintf('  ID          X              Y\n')
for i=1:N
    fprintf('   %d   %12.4e   %12.4e\n',Nodes(i,1:3))
end

fprintf('\n\nElements\n')
fprintf('  ID     MatID          Area        Node 1  Node 2\n')
for i=1:E
    fprintf('   %d        %d      %12.4e       %d       %d\n',Elems(i,1:5))
end

tol = 1e-6;
xmax = max( abs(Nodes(:,2) ) );
ymax = max( abs(Nodes(:,3) ) );
% if ymax < tol*xmax
%     iplot = 1;
% else
%     iplot = 0;
% end

% Read material info
Mats = load('Materials.txt');
[M,l] = size(Mats);

fprintf('\n\nMaterials\n')
fprintf(' MatID         Y              Rho\n')
for i=1:M
    fprintf('   %d     %12.4e   %12.4e\n',Mats(i,1:3))
end

if inpflg == 0
    % Apply load on right end node in x-direction
    P_hat = 1;
else
    % Read Essential and Natural BCs
    EBC = load('EBC.txt');
    [P,l] = size(EBC);
    NBC = load('NBC.txt');
    [Q,l] = size(NBC);
end

% Determine total number of degrees-of-freedom
udof = 2;        % Displacement degrees-of-freedom per node
ndeg = N*udof;   % Total degrees-of-freedom in problem
eldof = 4;       % Total degrees-of-freedom per element

% Initialize global arrays
kgmax = E*eldof*eldof;
Ig    = zeros(kgmax,1);
Jg    = zeros(kgmax,1);
Kg    = zeros(kgmax,1);
kg = 0;
% U = zeros(ndeg,1);       % Generalized displacement vector
F = zeros(ndeg,1);       % Generalized force vector

TotalMass = 0;
Mass  = zeros(E,1);

% Loop over Bar Elements
for e = 1:E
    
    % Establish element connectivity and coordinates
    Nnums = Elems(e,4:5);
    xy = Nodes(Nnums(:),2:3);
    
    % Extract element cross-sectional area
    A = Elems(e,3);
    
    % Extract element Young's modulus and mass density
    Y   = Mats(Elems(e,2),2);
    Rho = Mats(Elems(e,2),3);

    % Construct element generalized stiffness matrix
    [Kuu,L] = El_Stiff(xy,A,Y);
    % Kuu
    
    % Determine element total mass
    Mass(e) = Rho*A*L;

    % Preassemble generalized stiffness matrix coefficients
    [kl,Il,Jl,Kl] = El_Assembly(N,Nnums,Kuu);
    kg0 = kg;
    kg  = kg + kl;
    Ig(kg0+1:kg) = Il;
    Jg(kg0+1:kg) = Jl;
    Kg(kg0+1:kg) = Kl;
    
    % total mass
    TotalMass = TotalMass + Mass(e);

end

fprintf('\n\nTotal Mass = %12.4e\n\n',TotalMass)

% Form global sparse generalized stiffness matrix
KG = sparse(Ig(1:kg),Jg(1:kg),Kg(1:kg),ndeg,ndeg);
if iplot > 0
    figure(101)
    title('Before applying constraints')
    spy(KG)
end

if inpflg == 0
    % Enforce u=0 conditions on left boundary node
    KG(1,:) = 0;
    KG(:,1) = 0;
    KG(1,1) = 1;
    % Enforce v=0 conditions on all nodes
    KG(2:2:ndeg,2:2:ndeg) = eye(N);
else
    for i=1:P
        nebc = EBC(i,1);
        debc = EBC(i,2);
        p = udof*(nebc-1)+debc;
        KG(p,:) = 0;
        KG(:,p) = 0;
        KG(p,p) = 1;
        F(p) = EBC(i,3);
    end
end
if iplot > 0
    figure(102)
    title('After applying constraints')
    spy(KG)
end

% Construct global force vector on right boundary node
if inpflg == 0
    F(ndeg-1) = P_hat;
    F(ndeg/2) = P_hat/2;
else
    for i=1:Q
        nnbc = NBC(i,2);
        dnbc = NBC(i,3);
        q = udof*(nnbc-1)+dnbc;
        F(q) = F(q) + NBC(i,4);
    end
end
% svdKG = svds(KG,10)
% KG
% F

% Solve system to determine displacements
U = KG\F;

fprintf('\n\nSOLUTION\n\nNodal Displacements\n\n')
fprintf('  ID         Ux             Uy\n')
for i=1:N
    fprintf('   %d   %12.4e   %12.4e\n',Nodes(i,1),U(2*i-1),U(2*i))
end

Uabs = zeros(N,1);
for i=1:N
    i0 = (i-1)*udof;
    Uabs(i) = sqrt( U(i0+1)^2 + U(i0+2)^2 );
end
Umax = max(Uabs);
fprintf('\nMaximum Nodal Displacement = %12.4e   %12.4e\n',Umax)

% if inpflg == 0
xmax = max(Nodes(:,2));
xmin = min(Nodes(:,2));
ymax = max(Nodes(:,3));
ymin = min(Nodes(:,3));
xysize = sqrt( (xmax-xmin)^2 + (ymax-ymin)^2 );
ufct = 0.05*xysize/Umax;
ux   = zeros(N,1);
uy   = zeros(N,1);
umag = zeros(N,1);
for k=1:N
    kloc = (k-1)*udof;
    ux(k) = U(kloc+1);
    uy(k) = U(kloc+2);
    umag(k) = sqrt( ux(k)^2 + uy(k)^2 );
end

if iplot > 0

    xp0 = zeros(2,E);
    yp0 = zeros(2,E);
    xp  = zeros(2,E);
    yp  = zeros(2,E);
    for e = 1:E
        % Establish element connectivity and coordinates
        Nnums = Elems(e,4:5);
        xy = Nodes(Nnums(:),2:3);
        for j=1:2
            node = Nnums(j);
            j1 = udof*(node-1)+1;
            j2 = j1+1;
            xp0(j,e) = xy(j,1);
            yp0(j,e) = xy(j,2);
            xp(j,e) = xy(j,1)+ufct*U(j1);
            yp(j,e) = xy(j,2)+ufct*U(j2);
        end
    end

    figure(1)
    for e = 1:E
        plot(xp0(:,e),yp0(:,e),'b','LineWidth',2)
        hold on
    end
    hold off
    axis equal

    if iplot > 1
        figure(2)
        for e = 1:E
            plot(xp0(:,e),yp0(:,e),'b--','LineWidth',2)
            hold on
            plot(xp(:,e),yp(:,e),'r','LineWidth',2)
            hold on
        end
        hold off
        axis equal
    end

    if iplot > 2
        figure(3)
        plot(Nodes(:,2),ux,'LineWidth',2)
        title('x-displacements')
        xlabel('x')
        ylabel('U_x')

        figure(4)
        plot(Nodes(:,2),uy,'LineWidth',2)
        title('y-displacements')
        xlabel('x')
        ylabel('U_y')
    end
end

% Recover nodal and element axial strains and stresses
Dsp = zeros(E,4);
Eps = zeros(E,1);
Sig = zeros(E,1);
% if inpflg == 0
    epsxx = zeros(N,1);
    sigxx = zeros(N,1);
    nnel  = zeros(N,1);
% end

fprintf('\n\n\nElement Displacements\n\n')
fprintf('  ID    Node 1      Ux             Uy         Node 2      Ux             Uy\n')
for e = 1:E
    
    % Establish element connectivity and coordinates
    Nnums = Elems(e,4:5);
    xy = Nodes(Nnums(:),2:3);
    
    % Extract element cross-sectional area
    A = Elems(e,3);
    
    % Extract element elastic Young's modulus and Poisson's ratio
    Y = Mats(Elems(e,2),2);
    
    % Extract element nodal displacements
    inode1 = Nnums(1);
    inode2 = Nnums(2);
    Dsp(e,1) = U(udof*(inode1-1)+1);
    Dsp(e,2) = U(udof*inode1);
    Dsp(e,3) = U(udof*(inode2-1)+1);
    Dsp(e,4) = U(udof*inode2);
    fprintf('   %d      %d   %12.4e   %12.4e       %d   %12.4e   %12.4e\n',...
        e,Nnums(1),Dsp(e,1:2),Nnums(2),Dsp(e,3:4))
    
    % Find element strains and stresses
    ue = Dsp(e,:)';
    [eps,sig] = El_Stress(e,xy,ue,A,Y);
%     ue
%     eps
%     sig
    
    % Store element strains
    Eps(e) = eps;
    
    % Store element stresses
    Sig(e) = sig;
    
    % if inpflg == 0
        for ij=1:2
            k = Nnums(ij);
            epsxx(k) = epsxx(k) + eps;
            sigxx(k) = sigxx(k) + sig;
            nnel(k) = nnel(k) + 1;
        end
    % end
    
end

if inpflg == 0
    % plot strains and force-stresses
    for i=1:N
        epsxx(i) = epsxx(i)/nnel(i);
        sigxx(i) = sigxx(i)/nnel(i);
    end
    
    if iplot > 0
        figure(11)
        plot(Nodes(:,2),epsxx,'LineWidth',2)
        title('strains-xx')
        figure(21)
        plot(Nodes(:,2),sigxx,'LineWidth',2)
        title('stresses-xx')
    end
    
    % nnel

    epsxxmax = max(epsxx)
    epsxxmin = min(epsxx)
    epsxx
    sigxxmax = max(sigxx)
    sigxxmin = min(sigxx)
    sigxx
    
end

fprintf('\n\nElement Strains\n\n')
fprintf('  ID          Eps\n')
for e=1:E
    fprintf('   %d    %12.4e\n',e,Eps(e,1))
end

fprintf('\n\nElement Stresses\n\n')
fprintf('  ID        Sig\n')
for e=1:E
    fprintf('   %d    %12.4e\n',e,Sig(e,1))
end
fprintf('\n')

toc
