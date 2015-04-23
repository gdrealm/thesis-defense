function[] = main()
% Main program
% Modify ex, ey, and ez to change shape of element
% Change levs to change level-set configuration

% Global x,y,z coordinates of element.
% This will form a 3D cubic element cube.
ex = [-0.5 -0.5 -0.5 -0.5 +0.5 +0.5 +0.5 +0.5];
ey = [-0.5 -0.5 +0.5 +0.5 -0.5 -0.5 +0.5 +0.5];
ez = [+0.5 -0.5 -0.5 +0.5 +0.5 -0.5 -0.5 +0.5];

%% Extract a particular case (set k manually).

% Initial test for a particular levelset function.
% Randn function will return a m x n matrix of random positive and negative
% numbers.
hexsect = cell(127,8);
k = round(1 + (127-1).*rand);
% levs = randn(1,8);
% levs = [-1 -1 -1 -1 -1 -1 -1 1];
levs = [+1 -1 +1 -1 -1 +1 -1 +1];
[nsct,isct,xsct,ysct,zsct,levs] = xfem8isct(ex,ey,ez,levs);
[xtet,ytet,ztet,ctet,ptet,pnd,plist,Tetp1,Tetp2,tet1G,tet1L,tet2G,tet2L,size_phase1,size_phase2,volume1,volume2,total_volume1,total_volume2] = xfem8tet(isct,xsct,ysct,zsct,ex,ey,ez,levs);
hexsect{k,2} = xtet;
hexsect{k,3} = ytet;
hexsect{k,4} = ztet;
hexsect{k,5} = ctet;
hexsect{k,6} = ptet;
hexsect{k,7} = pnd;
hexsect{k,8} = plist;
[max(plist) min(plist)];

figure(1)
tetramesh(Tetp1,[xtet(1,:)',ytet(1,:)',ztet(1,:)'],-ones(size(Tetp1,1),1));
xlabel('X')
ylabel('Y')
zlabel('Z')
grid
axis equal

figure(2)
tetramesh(Tetp2,[xtet(1,:)',ytet(1,:)',ztet(1,:)'],ones(size(Tetp2,1),1));
xlabel('X')
ylabel('Y')
zlabel('Z')
grid
axis equal

figure(3)
tetramesh(ctet,[xtet(1,:)',ytet(1,:)',ztet(1,:)'],ptet);
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal

figure(4)
plot3(xtet(1,:),ytet(1,:),ztet(1,:),'X')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal

end