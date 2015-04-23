function [xtet,ytet,ztet,ctet,ptet,pnd,plist,Tetp1,Tetp2,tet1G,tet1L,tet2G,tet2L,size_phase1,size_phase2,volume1,volume2,total_volume1,total_volume2] = xfem8tet(isct,xsct,ysct,zsct,ex,ey,ez,levs)

% Intersection of hex8 based on level-set values.
%
% Input:  nsct   :  number of insected edges
%         isect  :  vector flags intersected edges
%         xsect  :  x-coordinates of insections in global and local coordinates
%         ysect  :  y-coordinates of insections in global and local coordinates
%         zsect  :  z-coordinates of insections in global and local coordinates
%         ex     :  global x-coordinates of nodes
%         ey     :  global y-coordinates of nodes
%         ez     :  global z-coordinates of nodes
%         levs   :  levs
% 
% Output: xtet   :  x-coordinates of nodes of triangulated element
%         ytet   :  y-coordinates of nodes of triangulated element
%         ztet   :  z-coordinates of nodes of triangulated element
%         ctet   :  connectivity of tets in triangulated element 
%         ptet   :  phase of each tetrahedron
%         pnd    :  nodal levelset values (actual values), # of nodes
%         plist  :  phase of each node (in terms of -1, 1 and 0)
%         Tetp1  :  nodes of tetrahedrons for main phase 1
%         Tetp2  :  nodes of tetrahedrons for main phase 2
%         tet1G  :  global coordinates of tetrahedron for phase 1
%         tet1L  :  local coordinates of tetrahedron for phase 1
%         tet2G  :  global coordinates of tetrahedron for phase 1
%         tet2L  :  local coordinates of tetrahedron for phase 2
%         size_phase1 : number of tetrahedrons in phase 1
%         size_phase2 : number of tetrahedrons in phase 2
%         volume1     : volume of a tetrahedron in phase 1
%         volume2     : volume of a tetrahedron in phase 2
%         total_volume1 : total volume of all tetrahedrons in phase 1
%         total_volume2 : total volume of all tetrahedrons in phase 2

%% Identify coordinates and phases to triangulate

% Nodes in local coordinates
xp = [-1 -1 -1 -1  1  1  1  1];
yp = [-1 -1  1  1 -1 -1  1  1];
zp = [1  -1  -1 1 1  -1  -1 1];

% Find intersected edges. Order of edges depends on configuration of edgmap
% isct flagged intersected edges with a 1
itr = find(isct>0);              
% Sort nodes by main phase. Determine which nodes have positive or negative
% level-set values.
ip1 = find(levs<0);
ip2 = find(levs>0);

% Id-numbers of node at edge intersections, i.e. create pseudonodes
ips = 9:8+length(itr);

% Create phase vector for triangulated element. Vector contains level-set
% values at the 8 original nodes plus zero values for the new pseudonodes
pnd = [levs zeros(1,length(itr))];

% Create nodes for main phase 1
% Coordinates where nodes are negative plus coordinates of intersections in
% global coordinates
exp1 = [ex(ip1) xsct(itr,1)'];  
eyp1 = [ey(ip1) ysct(itr,1)'];
ezp1 = [ez(ip1) zsct(itr,1)'];
% Nodes and pseudonodes numbers of main phase 1
ipx1 = [ip1 ips];

% Create nodes for main phase 2
% Coordinates where nodes are positive plus coordinates of intersections in
% global coordinates
exp2 = [ex(ip2) xsct(itr,1)'];
eyp2 = [ey(ip2) ysct(itr,1)'];
ezp2 = [ez(ip2) zsct(itr,1)'];
% Nodes and pseudonodes numbers of main phase 2
ipx2 = [ip2 ips];

% Combine triangulation points of main phase 1 and 2 in local and global
% coordinates
xtet = [ex xsct(itr,1)';xp xsct(itr,2)'];
ytet = [ey ysct(itr,1)';yp ysct(itr,2)'];
ztet = [ez zsct(itr,1)';zp zsct(itr,2)'];

% Same as xtet, ytet, ztet, but only with global coordinates. Useful for
% triangulation below.
exp = [ex xsct(itr,1)'];
eyp = [ey ysct(itr,1)']; 
ezp = [ez zsct(itr,1)']; 

%% Triangulate main phase 1 and 2 together
% If we triangulate main phase 1 and 2 separately,tetrahedron will
% superpose for some level-set combinations. We need to triangulate the
% entire element as a whole.

Tp = DelaunayTri(exp',eyp',ezp');
Tetp = Tp.Triangulation(:,:);

%% Separate triangulation into main phases 1 and 2

Tetp1 = zeros(1,4);
Tetp2 = zeros(1,4);
index1 = 1;
index2 = 1;
% If tetrahedron contains a negative node, it is phase 1. Phase 2,
% otherwise.
for i = 1:size(Tetp,1)
        if pnd(Tetp(i,1)) < 0||pnd(Tetp(i,2)) < 0||pnd(Tetp(i,3)) < 0||pnd(Tetp(i,4)) < 0
            Tetp1(index1,:) = Tetp(i,:);
            index1 = index1 + 1;
        else
            Tetp2(index2,:) = Tetp(i,:);
            index2 = index2 +1;
        end
end

%% Display coordinates of phase 1 tetrahedrons
% Display coordinates of phase 1 tetrahedrons in global coordinates
% Display volume of each phase 1 tetrahedron
total_volume1 = 0;
display('The global coordinates for the phase 1 tetrahedrons are: ')
for i = 1:size(Tetp1,1)
    tet1G = [xtet(1,Tetp1(i,1)),ytet(1,Tetp1(i,1)),ztet(1,Tetp1(i,1));xtet(1,Tetp1(i,2)),ytet(1,Tetp1(i,2)),ztet(1,Tetp1(i,2));xtet(1,Tetp1(i,3)),ytet(1,Tetp1(i,3)),ztet(1,Tetp1(i,3));xtet(1,Tetp1(i,4)),ytet(1,Tetp1(i,4)),ztet(1,Tetp1(i,4))];
    disp(tet1G)
    % Calculate volume. Algorithm: For 4 points in tetrahedron P,Q,R,S,
    % volume = abs(det([Q-P;R-Q;S-R]))/6
    a = tet1G(2,:)-tet1G(1,:);
    b = tet1G(3,:)-tet1G(2,:);
    c = tet1G(4,:)-tet1G(3,:);
    volume1 = abs(det([a;b;c]))/6;
    total_volume1 = total_volume1 + volume1;
    display('The volume of phase 1 tetrahedron is: ')
    disp(volume1)
end

% Display coordinates of phase 1 tetrahedrons in local coordinates
display('The local coordinates for the phase 1 tetrahedrons are: ')
for i = 1:size(Tetp1,1)
    tet1L = [xtet(2,Tetp1(i,1)),ytet(2,Tetp1(i,1)),ztet(2,Tetp1(i,1));xtet(2,Tetp1(i,2)),ytet(2,Tetp1(i,2)),ztet(2,Tetp1(i,2));xtet(2,Tetp1(i,3)),ytet(2,Tetp1(i,3)),ztet(2,Tetp1(i,3));xtet(2,Tetp1(i,4)),ytet(2,Tetp1(i,4)),ztet(2,Tetp1(i,4))];
    disp(tet1L)
end

% Display total volume of phase 1 tetrahedrons
    display('The total volume of phase 1 tetrahedrons is: ')
    disp(total_volume1)

%% Display coordinates of phase 2 tetrahedrons
% Display coordinates of phase 2 tetrahedrons in global coordinates
% Display volume of each phase 2 tetrahedron
total_volume2 = 0;
display('The global coordinates for the phase 2 tetrahedrons are: ')
for i = 1:size(Tetp2,1)
    tet2G = [xtet(1,Tetp2(i,1)),ytet(1,Tetp2(i,1)),ztet(1,Tetp2(i,1));xtet(1,Tetp2(i,2)),ytet(1,Tetp2(i,2)),ztet(1,Tetp2(i,2));xtet(1,Tetp2(i,3)),ytet(1,Tetp2(i,3)),ztet(1,Tetp2(i,3));xtet(1,Tetp2(i,4)),ytet(1,Tetp2(i,4)),ztet(1,Tetp2(i,4))];
    disp(tet2G)
    % Calculate volume.
    a = tet2G(2,:)-tet2G(1,:);
    b = tet2G(3,:)-tet2G(2,:);
    c = tet2G(4,:)-tet2G(3,:);
    volume2 = abs(det([a;b;c]))/6;
    total_volume2 = total_volume2 + volume2;
    display('The volume of phase 2 tetrahedron is: ')
    disp(volume2)
end

% Display coordinates of phase 2 tetrahedrons in local coordinates
display('The local coordinates for the phase 2 tetrahedrons are: ')
for i = 1:size(Tetp2,1)
    tet2L = [xtet(2,Tetp2(i,1)),ytet(2,Tetp2(i,1)),ztet(2,Tetp2(i,1));xtet(2,Tetp2(i,2)),ytet(2,Tetp2(i,2)),ztet(2,Tetp2(i,2));xtet(2,Tetp2(i,3)),ytet(2,Tetp2(i,3)),ztet(2,Tetp2(i,3));xtet(2,Tetp2(i,4)),ytet(2,Tetp2(i,4)),ztet(2,Tetp2(i,4))];
    disp(tet2L)
end

% Display total volume of phase 2 tetrahedrons
    display('The total volume of phase 2 tetrahedrons is: ')
    disp(total_volume2)
    
%% Display number of tetrahedrons in each phase.
size_phase1 = size(Tetp1(:,1));
size_phase1 = size_phase1(1);
size_phase2 = size(Tetp2(:,1));
size_phase2 = size_phase2(1);
display('The number of tetrahedrons in phase 1 is: ')
disp(size_phase1);
display('The number of tetrahedrons in phase 2 is: ')
disp(size_phase2);

%% Locate triangle interfaces - new algorithm
% Use function ismember to compare array vectors
% Display only result if three nodes repeat

display('The interfaces can be found on the triangles with nodes: ')
for i = 1:size(Tetp1,1)
    for j = 1:size(Tetp2,1)
        r = ismember(Tetp1(i,:),Tetp2(j,:));
        Tetp1_tri = Tetp1(i,:);
        Tetp1_tri = Tetp1_tri(r);
        if size(Tetp1_tri,2) == 3
            disp(Tetp1_tri)
        end
    end
end            
            
%% Combine triangulation of nodes in ctet
% Size of ptet depends on number of tetrahedrons
% Phase 1 tetrahedrons have value of -1, phase 2 value of 1 in ptet

% Original:
% ctet=[Tetp1;Tetp2];
% ptet=[-ones(size(Tetp1,1),1);ones(size(Tetp2,1),1)];

% New method:
% Provide two options for triangulation
% Option 1: based on volume
% If one phase is significantly larger than the other, switch tetramesh
if total_volume2 <= (0.2*(total_volume1 + total_volume2))
    ctet=[Tetp2;Tetp1];
    ptet=[-ones(size(Tetp2,1),1);ones(size(Tetp1,1),1)];
else
    ctet=[Tetp1;Tetp2];
    ptet=[-ones(size(Tetp1,1),1);ones(size(Tetp2,1),1)];
end

% Option 2: based on user's choice
% option = input('Choose triangulation option 1 or 2: ');
% if option == 1
%     ctet=[Tetp2;Tetp1];
%     ptet=[-ones(size(Tetp2,1),1);ones(size(Tetp1,1),1)];
% elseif
%     ctet=[Tetp1;Tetp2];
%     ptet=[-ones(size(Tetp1,1),1);ones(size(Tetp2,1),1)];
% end

%% Check connectivity between phases. PART of ORIGINAL CODE

% Obsolete? plist produces the same value as pnd
% Check connectivity for main phase 1 and creat sub-phase information

% Matrix plisp1 formed of (nodes + pseudonodes) x 1 elements with value of
% -2. ipx1 represents location of phase 1 nodes + pseudonodes. At these
% locations, values are replaced by -1
plisp1 = -2*ones(size(xtet,2),1);
plisp1(ipx1) = -1;

domp1 = 0;
% While there are phase 1 nodes
while ismember(-1,plisp1) > 0
    domp1 = domp1 + 1;
    % Find where the negative nodes are in the nodes + pseudonodes vector
    % Value are progressively changed by 0, one at a time
    ppp = find(plisp1 == -1);
    plisp1(ppp(1)) = 0;
    % Locate where values were changed to zero and create new variable ppp
    % pid represents the node where value has transformed into zero
    while ismember(0,plisp1) > 0
        ppp = find(plisp1 == 0);
        pid = ppp(1);  
        for it = 1:size(ctet,1)
            % If the node at pid belongs to Tetp1 and is negative
            if ismember(pid,ctet(it,:)) > 0 && ptet(it) < 0             
                plisp1(ctet(it,:)) = max(0,plisp1(ctet(it,:)));
            end
        end
        % At the location of phase 1 nodes + pseudonodes, the value will be
        % replaced to 1.
        plisp1(pid) = domp1;
    end
end

% Check connectivity for main phase 2 and creat sub-phase information
% Functions the same as previous routine, but for phase 2
plisp2 = -2*ones(size(xtet,2),1);
plisp2(ipx2) = -1;

domp2 = 0;
while ismember(-1,plisp2) > 0
    domp2 = domp2+1;
    ppp = find(plisp2 == -1);
    plisp2(ppp(1)) = 0;
    while ismember(0,plisp2) > 0
        ppp = find(plisp2 == 0);
        pid = ppp(1);  
        for it = 1:size(ctet,1)
            if ismember(pid,ctet(it,:)) > 0 && ptet(it) > 0
                plisp2(ctet(it,:)) = max(0,plisp2(ctet(it,:)));
            end
       end
       plisp2(pid) = domp2;
    end
end

% Build phase list including subphase information
% Replace all pseudonodes by 0
plisp1(ips) = 0;
plisp2(ips) = 0;

% Find location of the pseudonodes
idp1 = find(plisp1>0);
idp2 = find(plisp2>0);

% Create a zero (nodes + pseudonodes) x 1 matrix plist
plist = zeros(size(xtet,2),1);
% Create list that shows which original nodes belong to phase I and II
plist(idp1) = -plisp1(idp1);
plist(idp2) = plisp2(idp2);

end