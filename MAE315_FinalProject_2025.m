function MAE315_FinalProject_2025()
% ============================================================
% MAE315 Final Project (CST Plane Stress) - Screenshot Version
% Put in same folder as NodeFile.txt and ConnectFile.txt
% Run:
% >> MAE315_FinalProject_2025
% ============================================================

clc; close all;

%% -------------------- INPUTS --------------------
nodeFile = "NodeFile.txt";
connFile = "ConnectFile.txt";

% Plate / material (given)
b  = 0.1;      % thickness
w  = 2.0;      % width
H  = 3.0;      % height
E  = 20.0;     % Young's modulus
nu = 0.25;     % Poisson ratio

% Loads
p0     = 2.0;  % Part I uniform top traction magnitude (downward)
Fpoint = 0.4;  % Part II point force magnitude (downward)

% Plot-only toggle to match the handout figure's colorbar scale
MATCH_HANDOUT_SCALE = true;

tol = 1e-9;

%% -------------------- READ MESH --------------------
nodeData = readmatrix(nodeFile);      % [nodeID, x, y]
connData = readmatrix(connFile);      % [elemID, n1, n2, n3]

nodeIDs = nodeData(:,1);
nnode   = max(nodeIDs);

X = zeros(nnode,1);
Y = zeros(nnode,1);
X(nodeIDs) = nodeData(:,2);
Y(nodeIDs) = nodeData(:,3);

elemIDs   = connData(:,1);
elemNodes = connData(:,2:4);
nelem     = size(elemNodes,1);

ndof = 2*nnode;

%% -------------------- MATERIAL MATRIX (PLANE STRESS) --------------------
% Engineering shear convention (gamma_xy) => D(3,3) has /2
D = (E/(1-nu^2)) * [ 1,  nu, 0;
                    nu,  1,  0;
                     0,  0, (1-nu)/2 ];

% ============================================================
% PART I (a): CODE USED TO COMPUTE B MATRIX + PRINT B FOR ELEM #1
% Screenshot this block + the printed output in the command window
% ============================================================

%%ASSEMBLE GLOBAL STIFFNESS
K = sparse(ndof,ndof);

for e = 1:nelem
    n = elemNodes(e,:);
    xp = X(n);
    yp = Y(n);

    % B matrix computation 
    [Bp, Area] = CalcBmatTri(xp,yp);

    % Print B for element #1 
    if elemIDs(e) == 1
        fprintf("====================================================\n");
        fprintf("PART I (a): B matrix for element #1 (Area = %.6f)\n", Area);
        disp(Bp);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ke = (Bp.' * D * Bp) * Area * b;

    dofs = [2*n(1)-1; 2*n(1);
            2*n(2)-1; 2*n(2);
            2*n(3)-1; 2*n(3)];

    K(dofs,dofs) = K(dofs,dofs) + ke;
end
% BOUNDARY CONDITIONS 

bottomNodes = find(abs(Y - 0) < tol);

centerBottomCandidates = bottomNodes(abs(X(bottomNodes) - w/2) < tol);

centerBottom = centerBottomCandidates(1);

fixedDOF = [];
fixedVal = [];

fixedDOF = [fixedDOF; 2*bottomNodes];                % v=0 at all bottom nodes
fixedVal = [fixedVal; zeros(numel(bottomNodes),1)];

fixedDOF = [fixedDOF; 2*centerBottom - 1];           % u=0 at bottom center
fixedVal = [fixedVal; 0];

[fixedDOF, ia] = unique(fixedDOF, "stable");
fixedVal = fixedVal(ia);

% ============================================================
% PART I (b): CODE USED TO GENERATE FORCING TERM F FOR KU=F (TRACTION)
% Screenshot this block
% ============================================================

F1 = zeros(ndof,1);

boundaryE = boundaryEdges(elemNodes);
isTopEdge = abs(Y(boundaryE(:,1)) - H) < tol & abs(Y(boundaryE(:,2)) - H) < tol;
topEdges  = boundaryE(isTopEdge,:);

for k = 1:size(topEdges,1)
    i = topEdges(k,1);
    j = topEdges(k,2);

    L  = hypot(X(j)-X(i), Y(j)-Y(i));
    fy = -p0 * b * L/2;      % consistent nodal load per end node

    F1(2*i) = F1(2*i) + fy;
    F1(2*j) = F1(2*j) + fy;
end

% ============================================================
% PART I: SOLVE + PLOT + DISPLACEMENT COMPARISON
% (c) plot original vs deformed
% (d) compare FEM max displacement to delta = F H/(A E)
% ============================================================

U1 = applyDirichletSolve(K,F1,fixedDOF,fixedVal);

% ---------- PART I (c): Original vs Deformed mesh plot ----------


figure("Name","Original vs Deformed Mesh");
plotMeshOverlay(elemNodes,X,Y,U1,1.0);
title(" Original (black) and Deformed (red) mesh");
xlabel("x"); ylabel("y"); axis equal;


% ---------------------------------------------------------------

% ---------- PART I (d): compare displacement to delta = F H/(A E) ----------


v1 = U1(2:2:end);
vmin_FEM = min(v1);

Ftotal = sum(F1(2:2:end));  
Across  = w*b;
delta_MoM = (Ftotal*H)/(Across*E);
pctDiff = (vmin_FEM - delta_MoM)/delta_MoM * 100;

fprintf("====================================================\n");
fprintf("PART I (d):\n");
fprintf("  Total applied force Ftotal = %.6f\n", Ftotal);
fprintf("  FEM min v (max downward)    = %.6f\n", vmin_FEM);
fprintf("  MoM delta = F H/(A E)       = %.6f\n", delta_MoM);
fprintf("  Percent difference          = %.6f %%\n", pctDiff);
% -------------------------------------------------------------------------

% ============================================================
% PART II: POINT FORCE AT TOP CENTER NODE (F=0.4)
% ============================================================

% ---------- (setup) find top-center node and apply point load ----------


topNodes = find(abs(Y - H) < tol);
topCenterCandidates = topNodes(abs(X(topNodes) - w/2) < tol);
topCenter = topCenterCandidates(1);

F2 = zeros(ndof,1);
F2(2*topCenter) = -Fpoint;
U2 = applyDirichletSolve(K,F2,fixedDOF,fixedVal);

% ============================================================
% PART II (a): Plot original node locations + displaced node locations
% Screenshot this block and the figure
% ============================================================


figure("Name","Part II: Nodes before/after deformation");
plotNodeOverlay(X,Y,U2,1.0);
title("Original nodes (black o) and displaced nodes (red x)");
xlabel("x"); ylabel("y"); axis equal;


% ============================================================
% PART II (b)(c): Compute and output ALL stress components
% (b) elements #3 and #4
% (c) elements #27 and #28
% Screenshot print blocks
% ============================================================

sigmaAll = zeros(nelem,3); % [sigx sigy tauxy]
for e = 1:nelem
    n = elemNodes(e,:);
    xp = X(n);
    yp = Y(n);

    [Bp, ~] = CalcBmatTri(xp,yp);

    dofs = [2*n(1)-1; 2*n(1);
            2*n(2)-1; 2*n(2);
            2*n(3)-1; 2*n(3)];

    ue = U2(dofs);

    eps = Bp * ue;
    sig = D * eps;

    sigmaAll(e,:) = sig.';
end

fprintf("====================================================\n");
fprintf("PART II (b): stresses in elements #3 and #4\n");
printStressForElements([3 4], elemIDs, sigmaAll);

fprintf("----------------------------------------------------\n");
fprintf("PART II (c): stresses in elements #27 and #28\n");
printStressForElements([27 28], elemIDs, sigmaAll);

% ============================================================1
% PART II (d): Compare sigma_y in those elements to nominal sigma_y0 = F/A
% Screenshot this block
% ============================================================

sigma_y0 = -Fpoint/(w*b); % compressive negative
fprintf("----------------------------------------------------\n");
fprintf("PART II (d): nominal sigma_y0 = F/A = %.6f\n", sigma_y0);

targets = [3 4 27 28];
for eID = targets
    idx = find(elemIDs==eID,1,"first");
    sigy = sigmaAll(idx,2);
    pct  = (sigy - sigma_y0)/sigma_y0 * 100;
    fprintf("  Element %2d: sigma_y = %+10.6f, percent diff = %+10.6f %%\n", eID, sigy, pct);
end

% ============================================================
% PART II (e): Compute stresses in ALL elements + "fill" plot colored by sigma_y
% Screenshot this block and the plot
% ============================================================

sigmaY = sigmaAll(:,2);

sigmaY_plot = sigmaY;
plotTitle = "\sigma_y per element (Part II)";
if MATCH_HANDOUT_SCALE
    sigmaY_plot = sigmaY * b; % plot-only scaling to match the sample figure scale
    plotTitle = "\sigma_y*b per element ";
end

figure("Name","Part II: fill plot");
patch('Faces',elemNodes, ...
      'Vertices',[X Y], ...
      'FaceColor','flat', ...
      'CData',sigmaY_plot, ...
      'EdgeColor','k');
axis equal; colorbar;
xlabel("x"); ylabel("y");
title(plotTitle);

% Optional CSV output
out = [elemIDs, sigmaAll]; % [elemID sigx sigy tauxy]
writematrix(out, "PartII_ElementStresses.csv");

fprintf("====================================================\n");
fprintf("Wrote: PartII_ElementStresses.csv\n");
fprintf("Done.\n");

end

%% ===================== LOCAL FUNCTIONS =====================

function U = applyDirichletSolve(K,F,fixedDOF,fixedVal)
ndof = size(K,1);
U = zeros(ndof,1);
U(fixedDOF) = fixedVal;

free = true(ndof,1);
free(fixedDOF) = false;

Kff = K(free,free);
Ff  = F(free) - K(free,~free)*U(~free);

U(free) = Kff \ Ff;
end

function edgesB = boundaryEdges(elemNodes)
nelem = size(elemNodes,1);
allE = zeros(3*nelem,2);

row = 1;
for e = 1:nelem
    n = elemNodes(e,:);
    E = [n(1) n(2);
         n(2) n(3);
         n(3) n(1)];
    allE(row:row+2,:) = sort(E,2);
    row = row + 3;
end

[uniqE,~,ic] = unique(allE,"rows");
counts = accumarray(ic,1);
edgesB = uniqE(counts==1,:);
end

function plotMeshOverlay(elemNodes,X,Y,U,scale)
u = U(1:2:end);
v = U(2:2:end);

Xd = X + scale*u;
Yd = Y + scale*v;

triplot(elemNodes, X,  Y,  "k-"); hold on;
triplot(elemNodes, Xd, Yd, "r-"); hold off;
legend("Original","Deformed","Location","best");
end

function plotNodeOverlay(X,Y,U,scale)
u = U(1:2:end);
v = U(2:2:end);
Xd = X + scale*u;
Yd = Y + scale*v;

plot(X, Y, "ko", "MarkerSize",6, "LineWidth",1); hold on;
plot(Xd,Yd,"rx", "MarkerSize",7, "LineWidth",1.2); hold off;
legend("Original nodes","Displaced nodes","Location","best");
end

function printStressForElements(elemList, elemIDs, sigmaAll)
for eID = elemList
    idx = find(elemIDs==eID,1,"first");
    if isempty(idx)
        error("Element ID %d not found.", eID);
    end
    sig = sigmaAll(idx,:);
    fprintf("  Element %2d: sigma_x = %+10.6f, sigma_y = %+10.6f, tau_xy = %+10.6f\n", ...
        eID, sig(1), sig(2), sig(3));
end
end

function [Bp,Area] = CalcBmatTri(xp,yp)
% Engineering shear convention: gamma_xy = du/dy + dv/dx
x = xp; y = yp;

dNadx = -(y(3)-y(2)) / ( (y(1)-y(2))*(x(3)-x(2))-(x(1)-x(2))*(y(3)-y(2)) );
dNady =  (x(3)-x(2)) / ( (y(1)-y(2))*(x(3)-x(2))-(x(1)-x(2))*(y(3)-y(2)) );

dNbdx = -(y(1)-y(3)) / ( (y(2)-y(3))*(x(1)-x(3))-(x(2)-x(3))*(y(1)-y(3)) );
dNbdy =  (x(1)-x(3)) / ( (y(2)-y(3))*(x(1)-x(3))-(x(2)-x(3))*(y(1)-y(3)) );

dNcdx = -(y(2)-y(1)) / ( (y(3)-y(1))*(x(2)-x(1))-(x(3)-x(1))*(y(2)-y(1)) );
dNcdy =  (x(2)-x(1)) / ( (y(3)-y(1))*(x(2)-x(1))-(x(3)-x(1))*(y(2)-y(1)) );

Bp = zeros(3,6);
Bp(1,[1 3 5]) = [dNadx dNbdx dNcdx];
Bp(2,[2 4 6]) = [dNady dNbdy dNcdy];
Bp(3,:)       = [dNady dNadx dNbdy dNbdx dNcdy dNcdx];

Area = 0.5 * abs( (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) );
end
