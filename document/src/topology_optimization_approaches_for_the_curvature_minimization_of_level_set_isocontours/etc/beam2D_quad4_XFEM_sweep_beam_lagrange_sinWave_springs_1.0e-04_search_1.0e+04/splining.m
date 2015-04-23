function splining

close all
clear all

addpath topo

sweeptype=1;  % sweep type
polytype=3;   % type of polygon
plotflag=0;   % plot fitted splines (0/1)

switch sweeptype
    
    case 1   % sweep over number of segments for constant radius
        
        % define polygon
        ctr=[0 0];                    % center of circle
        if polytype<3
            np=4:2:50;                % number segments
        else
            np=1:1:80;
        end
        rad=ones(length(np),1);       % radius of circle
        
    case 2   % sweep over radius for constant number of segments
        
        % define polygon
        ctr=[0 0];                  % center of circle
        rad=0.25:0.05:0.75;         % radius of circle 
        np=50*ones(length(rad),1);  % number segments

    case 3   % sweep over radius for constant length of segments
        
        % define polygon
        ctr=[0 0];                         % center of circle
        rad=0.1:0.1:2;                     % radius of circle 
        np=2*floor(ceil(8*rad/rad(1))/2);  % number segments
end

numit=length(rad);  % number of analysis

strene=zeros(numit,1);
strmid=zeros(numit,1);

for it=1:numit
    
    % compute points on polygon
    switch polytype
        case 1
            [xy,reflen]=poly_halfcircle(ctr,rad(it),np(it));
        case 2
            [xy,reflen]=poly_sinwave_exact(rad(it),np(it));
            reflen=1e-3;
        case 3
            [xy,reflen]=poly_sinwave_femdoc(np(it));
            reflen=1e-3;
        case 4
            [xy]=poly_custom(6);
            reflen=0.021776049413800434;
        otherwise
            error('wrong polytype');
    end
            
    % compute strain energy
    [strene(it),strmid(it)]=beamfit(xy,reflen,plotflag*it);
    disp(strene(it));
end

switch sweeptype
    
    case 1  % sweep over number of segments for constant radius
        
        figure(plotflag*numit+1)
        plot(np,strene)
        xlabel('number of segments')
        title('total strain energy');
        
        figure(plotflag*numit+2)
        plot(np,strmid)
        xlabel('number of segments')
        title('strain energy of middle segments');
        
    case 2  % sweep over radius for constant number of segments
        
        figure(plotflag*numit+1)
        plot(rad,strene)
        xlabel('radius')
        title('total strain energy');
        
        figure(plotflag*numit+2)
        plot(rad,strmid)
        xlabel('radius')
        title('strain energy of middle segments');

    case 3  % sweep over radius for constant length of segments
        
        figure(plotflag*numit+1)
        plot(rad,strene)
        xlabel('radius')
        title('total strain energy');
        
        figure(plotflag*numit+2)
        plot(rad,strmid)
        xlabel('radius')
        title('strain energy of middle segments');
end
end

%==========================================================================

function [xy,reflen]=poly_halfcircle(xyc,rad,np)
%
% generate points on semi-circle

xy=zeros(np+1,2);    % coordinate of points on polygon

dphi=pi/np;          % angle between points

reflen=dphi*rad;

for i=1:np+1
    
    phi=-pi/2+(i-1)*dphi;
    
    xy(i,:)=xyc+rad*[sin(phi) cos(phi)];
end

end
%==========================================================================

function [xy,reflen]=poly_sinwave_exact(amp,np)
%
% generate points on semi-circle

xy=zeros(np+1,2);    % coordinate of points on polygon

dx=2/np;          % angle between points

reflen=dx;

for i=1:np+1
        
    xl=(i-1)*dx;
    xy(i,:)=[xl amp*sin(pi*xl)];
end

end
%==========================================================================

function [xy,reflen]=poly_custom(np)
%
% generate points on semi-circle

xy=zeros(np+1,2);    % coordinate of points on polygon

% Custom points half-circle r = 0.055
% xy( 1, 1 ) = +1.411666666667000e+00;
% xy( 1, 2 ) = +0.000000000000000e+00;
% 
% xy( 2, 1 ) = +1.466666666666666e+00;
% xy( 2, 2 ) = +5.499999999994164e-02;
% 
% xy( 3, 1 ) = +1.521666666666883e+00;
% xy( 3, 2 ) = +0.000000000000000e+00;

% Custom points half-circle r = 0.075
xy( 1, 1 ) = +1.391666666667000e+00;
xy( 1, 2 ) = +0.000000000000000e+00;

xy( 2, 1 ) = +1.400000000000000e+00;
xy( 2, 2 ) = +2.011844635237504e-02;

xy( 3, 1 ) = +1.446548220313729e+00;
xy( 3, 2 ) = +6.666666666666664e-02;

xy( 4, 1 ) = +1.466666666666666e+00;
xy( 4, 2 ) = +7.500000000000000e-02;

xy( 5, 1 ) = +1.486785113019947e+00;
xy( 5, 2 ) = +6.666666666666665e-02;

xy( 6, 1 ) = +1.533333333333333e+00;
xy( 6, 2 ) = +2.011844635384316e-02;

xy( 7, 1 ) = +1.541666666667000e+00;
xy( 7, 2 ) = +0.000000000000000e+00;

end
    
%==========================================================================

function [xy,reflen]=poly_sinwave_femdoc(np)

% load polygon date
srcfile=['./topo/elem_topo_' num2str(np)];
trgfile=[srcfile '.m'];
exefile=['elem_topo_' num2str(np)];

% display([ '... processing ' srcfile]);

copyfile(srcfile,trgfile,'f');

% polydata
warning off
feval(exefile);
warning on
delete(trgfile);

keySet=cell2mat(keySet);

% build polygon

numele=size(elem_topo,1);
numnod=size(node_coor,1);

topo=zeros(numele,2);
elen=zeros(numele,1);
for ie=1:numele
    na=find(keySet==elem_topo(ie,1));
    nb=find(keySet==elem_topo(ie,2));
    topo(ie,:)=[na nb];
    elen(ie)=norm(node_coor(na)-node_coor(nb));
end

reflen=min(elen);

iproc=ones(numele,1);
xy=zeros(numnod,2);

iele=-1;
for in=1:numnod
    [ix,ip]=find(topo==in);
    if length(ix)==1
        iele=ix;
        if ip==2
            topo(ix,:)=topo(ix,[2 1]);
        end
        break;
    end
end

if iele<1
    error('start of polygon not found');
end
kk=1;
xy(kk,:)=node_coor(topo(iele,1),1:2);

while sum(iproc)>0
    if iele==0
        error('end of polygon reached');
    end
    nb=topo(iele,2);
    [ix,ip]=find(topo==nb);
    if length(ix)>2; error('polygon does not fit'); end
    kk=kk+1;
    xy(kk,:)=node_coor(nb,1:2);
    iproc(iele)=0;
    
    if length(ix)==1; continue; end
    ixx=find(ix~=iele);
    iele=ix(ixx);
    ipp=ip(ixx);
    if ipp==2;
        topo(iele,:)=topo(iele,[2 1]);
    end
end

if kk~=numnod; error('something went wrong'); end

% figure(99)
% plot(xy(:,1),xy(:,2));
% maxsmo=5;
% newxy=xy;
% for is=1:maxsmo
%     totcuv=0;
%     for in=1:numnod
%         [ix,ip]=find(topo==in);
%         if length(ix)==1; continue; end
%         ny=0;
%         ll=0;
%         for ii=1:length(ix) 
%             iele=ix(ii);
%             elen=norm(xy(topo(iele,1),:)-xy(topo(iele,2),:));
%             ip2 =1;  % other node
%             ip1 =2;  % current node
%             if ip(ii)==1; ip2=2; ip1=1; end
%             in2=topo(iele,ip2);
%             in1=topo(iele,ip1);
%             if in1~=in; error('messed up nodes'); end
%             ny=ny+(xy(in2,2)-xy(in1,2)/elen);
%             ll=ll+elen;
%         end
%         newxy(in,2)=xy(in,2)+1e-5*ny/ll;
%         totcuv=totcuv+abs(ny/ll);
%     end
%     xy=newxy;
% %     clf
% %     plot(xy(:,1),xy(:,2),'-o');
%     fprintf('smoothing step %d: total curvatuve = %e\n',is,totcuv);
% end
 
% for in=1:numnod
% xy(in,2)=(0.25+(np-1)/99*0.5)*sin(pi*xy(in,1));
% xy(in,2)=xy(in,2)+1e-2*randn;
% end
% 
% pp=polyfit(xy(:,1),xy(:,2),7);

% figure(99)
% clf

% dxy=zeros(numnod,2);

% for in=1:numnod
% % dxy(in,1)=xy(in,1);    
% % dxy(in,2)=polyval(pp,xy(in,1))-(1+(0.25+(np-1)/99*0.5)*cos(pi*xy(in,1)));
% xy(in,2)=polyval(pp,xy(in,1));
% end

% plot(dxy(:,1),dxy(:,2),'b');

end

%==========================================================================

function [strene,strmid]=beamfit(xy,reflen,ifig)
%
% build fe problem and compute strain energy

% set element properties
h=0.1*reflen;  % dimension of square cross-section

ep(1)=1e11;    % Youngs modulus
ep(2)=h^2;     % cross-section area
ep(3)=h^4/12;  % 2nd area moment of inertia
ep(4)=1e-2;     % spring stiffness

% number of nodes
numnod=size(xy,1);

% build fe problem
[re,jac,edof,sol]=build_fepro(xy,ep);

% number of dofs
numdof=size(re,1);

% solve fe problem with displacements set to zero
pdof=[]; %1:numnod;
fdof=1:numdof;
fdof(pdof)=[];
sol(fdof)=jac(fdof,fdof)\re(fdof);

% plot deformed beam
if ifig>0
    plot_defbeam(xy,sol,edof,ifig);
end
    
% compute strain energies
[strene,strmid]=comp_strene(xy,sol,edof,ep);

end
    
%==========================================================================

function [re,jac,edof,sol]=build_fepro(xy,ep)

% number of nodes
numnod=size(xy,1);

% number of elements
numele=numnod-1;

% number of dofs in fe problem 
% displacements and rotations: 3 for all internal nodes; 2 for end nodes
% lagrange multiplier: 1 for all internal nodes

numdof=4*(numnod-2)+2*3+numnod-2;

% allocate arrays
jac=zeros(numdof);
re=zeros(numdof,1);
sol=zeros(numdof,1);
ndof=zeros(numnod,3);
ldof=zeros(numnod,1);
edof=zeros(numele,6);

% assign dofs to nodes
ndof(1:numnod,1)  =         1:numnod;      % ux displacement dofs
ndof(1:numnod,2)  =   numnod+1:2*numnod;   % uy displacement dofs
ndof(2:numnod,3)  = 2*numnod+1:3*numnod-1; % left rotation
ndof(1:numnod-1,4)= 3*numnod:4*numnod-2;   % right rotation

% assign lagr dofs to internal nodes
ldof(2:numnod-1)=4*numnod-1:5*numnod-4;

% assign dofs to beam elements
for ie=1:numele
    edof(ie,:)=[ndof(ie,[1 2 4]) ndof(ie+1,[1 2 3])];
end

% assemble elemental residual and jacobian
for ie=1:numele
    nids=[ie ie+1];                               % node ids
    edfs=edof(ie,:);                              % elemental dofs
    ke=beam2e(xy(nids,1),xy(nids,2),ep(1:3));     % elemental stiffness matrix
    jac(edfs,edfs)=jac(edfs,edfs)+ke;             % assemble into global stiffness matrix)
    elen=norm(xy(nids(1),:)-xy(nids(2),:));
    
    for in=1:2
        ndf=ndof(nids(in),1:2);
        jac(ndf,ndf)=jac(ndf,ndf)+elen*ep(4)*eye(2);
    end
end

% add spring stiffness
% for in=1:numnod
%     ndf=ndof(in,1:2);
%     jac(ndf,ndf)=jac(ndf,ndf)+ep(4)*eye(2);
% end

% assemble contributions from jumps in rotations
for in=2:numnod-1
    
    iel=in-1;            % left element
    ier=in;              % right elemnent
    
    nidsl=[iel iel+1];   % node ids of left element
    nidsr=[ier ier+1];   % node ids of left element
 
    % build normal vectors along elements from node A to node B
    vecl=[xy(nidsl(2),1)-xy(nidsl(1),1) xy(nidsl(2),2)-xy(nidsl(1),2)];
    vecr=[xy(nidsr(2),1)-xy(nidsr(1),1) xy(nidsr(2),2)-xy(nidsr(1),2)];

    vecl=vecl/norm(vecl);
    vecr=vecr/norm(vecr);
    
    mvec=0.5*(vecl+vecr);
    mvec=1/norm(mvec)*mvec;
    
    % compute angle between vectors
    rotvec1=cross([vecl 0],[mvec 0]);
    rotvec2=cross([mvec 0],[vecr 0]);
    sina1=rotvec1(3);
    sina2=rotvec2(3);
    
    % compute offsets in slope
    offset=sina1/sqrt(1-sina1^2)+sina2/sqrt(1-sina2^2);
    
    % add residual contributions
    jac(ndof(in,3),ldof(in))= 1;
    jac(ndof(in,4),ldof(in))=-1;
    jac(ldof(in),ndof(in,3))= 1;
    jac(ldof(in),ndof(in,4))=-1;
    re(ldof(in))            =offset;
end
    
end

%==========================================================================
function plot_defbeam(xy,sol,edof,ifig)

numnod=size(xy,1);

figure(ifig)
clf;
for ie=1:numnod-1
    ed=sol(edof(ie,:));                                  % elemental solution
    nids=[ie ie+1];                                      % node ids
    [excd,eycd]=beam2crd(xy(nids,1)',xy(nids,2)',ed',1); % compute points on deformed beam
    plot(excd',eycd');
    hold on;
    plot(xy(ie,1),xy(ie,2),'o'); hold on;                % end points
    plot(xy(ie+1,1),xy(ie+1,2),'o'); hold on;
end

axis equal

end

%==========================================================================

function [strene,strmid]=comp_strene(xy,sol,edof,ep)
% compute strain energy of entire beam

strene=0;
strmid=0;
lenmid=0;

numele=size(edof,1);

strpa=zeros(numele,3);

for ie=1:numele
    nids=[ie ie+1];                               % node ids
    edfs=edof(ie,:);                              % elemental dofs
    esol=sol(edfs);                               % elemental solution
    ke=beam2e(xy(nids,1),xy(nids,2),ep);          % elemental stiffness matrix
    se=0.5*esol'*ke*esol;                         % elemental strain energy
    elen=norm(xy(nids(1),:)-xy(nids(2),:));
    
    strpa(ie,:)=[elen se se/elen];
    
    strene=strene+se;
   
    disp( se );
    
    if ie == floor(numele/2) || ie == floor(numele/2)+1
        strmid=strmid+se;
        lenmid=lenmid+elen;
    end
        
end

strmid=strmid/lenmid;

% figure(99)
% clf
% plot(strpa(:,1),strpa(:,2),'o');

end

%==========================================================================

function [Ke,fe]=beam2e(ex,ey,ep,eq)
% Ke=beam2e(ex,ey,ep)
% [Ke,fe]=beam2e(ex,ey,ep,eq)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the stiffness matrix for a two dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%
%            ep = [E A I]       element properties
%                                  E: Young's modulus
%                                  A: Cross section area
%                                  I: Moment of inertia
%
%            eq = [qx qy]       distributed loads, local directions
% 
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%
%            fe : element load vector (6 x 1)
%--------------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;

  E=ep(1);  A=ep(2);  I=ep(3);
 
  qx=0; qy=0;  if nargin>3; qx=eq(1); qy=eq(2); end

  Kle=[E*A/L   0            0      -E*A/L      0          0 ;
         0   12*E*I/L^3   6*E*I/L^2  0   -12*E*I/L^3  6*E*I/L^2;
         0   6*E*I/L^2    4*E*I/L    0   -6*E*I/L^2   2*E*I/L;
       -E*A/L  0            0       E*A/L      0          0 ;
         0   -12*E*I/L^3 -6*E*I/L^2  0   12*E*I/L^3  -6*E*I/L^2;
         0   6*E*I/L^2    2*E*I/L    0   -6*E*I/L^2   4*E*I/L];
   
  fle=L*[qx/2 qy/2 qy*L/12 qx/2 qy/2 -qy*L/12]';

  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];

  Ke=G'*Kle*G;   fe=G'*fle; 
%--------------------------end--------------------------------
end
%==========================================================================

function [excd,eycd]=beam2crd(ex,ey,ed,mag)
%-------------------------------------------------------------
%    PURPOSE
%     Calculate the element continous displacements for a 
%     number of identical 2D Bernoulli beam elements. 
%  
%    INPUT:  ex,ey,
%            ed,
%            mag 
%
%    OUTPUT: excd,eycd 
%-------------------------------------------------------------

% LAST MODIFIED: P-E AUSTRELL 1993-10-15
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%   
  [nie,ned]=size(ed);
  
  for i=1:nie
  
      b=[ex(i,2)-ex(i,1) ey(i,2)-ey(i,1)];   L=sqrt(b*b');   n=b/L;
 
      G=[n(1) n(2)  0    0    0   0;
        -n(2) n(1)  0    0    0   0;
          0    0    1    0    0   0;
          0    0    0   n(1) n(2) 0;
          0    0    0  -n(2) n(1) 0; 
          0    0    0    0    0   1];
   
     d=ed(i,:)';   dl=G*d ;
 %
     xl=[0:L/20:L]';    one=ones(size(xl));
 %
     Cis=[-1 1;
           L 0]/L; ds=[dl(1);dl(4)];
 %          
     ul=([xl one]*Cis*ds)';      
 %          
     Cib=[ 12   6*L    -12  6*L;
          -6*L -4*L^2  6*L -2*L^2;
           0    L^3     0     0;
          L^3     0     0     0]/L^3; 
 %           
     db=[dl(2);dl(3);dl(5);dl(6)];
 %          
     vl=([xl.^3/6 xl.^2/2 xl one]*Cib*db)';
 %   
     cld=[ ul ;
           vl ]; A=[n(1) -n(2);
                    n(2) n(1)];  cd=A*cld;
 %
     xyc=A(:,1)*xl'+ [ex(i,1);ey(i,1)]*one';
 %                 
     excd(i,:)=xyc(1,:)+mag*cd(1,:);   
     eycd(i,:)=xyc(2,:)+mag*cd(2,:);     
  end
  
end