function [nsct,isct,xsct,ysct,zsct,levs] = xfem8isct(ex,ey,ez,levs)

% Intersection of hex8 based on level-set values.
%
% Input:   ex    : global x-coordinates of nodes of 3D element
%          ey    : global y-coordinates of nodes of 3D element
%          ez    : global z-coordinates of nodes of 3D element
%          levs  : level set values randomly generated
% 
% Output:  nsct  : number of intersected edges
%          isect : vector flags of intersected edges
%          xsect : x-coordinates of intersections in global and local coordinates
%          ysect : y-coordinates of intersections in global and local coordinates
%          zsect : z-coordinates of intersections in global and local coordinates

% Map of node connections across edges in a 3D element.
% Edgmap is a 12 x 2 matrix.
edgmap = [ 1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5];

% Values of nodes in master element (local coordinates).

xp = [-1 -1 -1 -1  1  1  1  1];
yp = [-1 -1  1  1 -1 -1  1  1];
zp = [1  -1  -1 1 1  -1  -1 1];

% Set initial value of intersected edges to zero.
% Create a zero 12 x 1 matrix to flag edges with intersected edges.
% Create a zero 12 x 2 matrix to record x, y, z coordinates in global and
% local coordinates.
nsct = 0;
isct = zeros(12,1);
xsct = zeros(12,2);
ysct = zeros(12,2);
zsct = zeros(12,2);

for i = 1:12
    % ic1 and ic2 return the values of the first and second columns of
    % edgmap, respectively, for a determined row. ic1 and ic2 represent two
    % nodes connected by a cube edge.
    ic1 = edgmap(i,1);
    ic2 = edgmap(i,2);
    % The level set values at those nodes are multiplied to determined
    % intersection.
    if levs(ic1)*levs(ic2) < 0      % If so, then there is an intersection.
        nsct = nsct+1;              % Increase number of intersected edges by one.
        isct(i) = 1;                % Flag edge as having an intersection in isct matrix.
        sctr = -levs(ic1)/(levs(ic2)-levs(ic1));    % Dimensionless location (ratio) of intersection.
        % Insection in glb coordinates
        xsct(i,1) = ex(ic1)+sctr*(ex(ic2)-ex(ic1));
        ysct(i,1) = ey(ic1)+sctr*(ey(ic2)-ey(ic1));
        zsct(i,1) = ez(ic1)+sctr*(ez(ic2)-ez(ic1));
        % Insection in local coordinates
        xsct(i,2) = xp(ic1)+sctr*(xp(ic2)-xp(ic1));
        ysct(i,2) = yp(ic1)+sctr*(yp(ic2)-yp(ic1));
        zsct(i,2) = zp(ic1)+sctr*(zp(ic2)-zp(ic1));
    end
end

    display('Number of intersected edges: ');
    disp(nsct);
    display('Edges with intersections: ');
    disp(isct');
    display('x-coordinates of intersections in global coordinates: ');
    disp(xsct(:,1));
    display('x-coordinates of intersections in local coordinates: ');
    disp(xsct(:,2));
    display('y-coordinates of intersections in global coordinates: ');
    disp(ysct(:,1));
    display('y-coordinates of intersections in local coordinates: ');
    disp(ysct(:,2));
    display('z-coordinates of intersections in global coordinates: ');
    disp(zsct(:,1));
    display('z-coordinates of intersections in local coordinates: ');
    disp(zsct(:,2));

end