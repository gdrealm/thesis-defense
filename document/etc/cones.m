clc
clear all
close all

x = 0:0.05:3;
y = 0:0.05:2;

x = x';
y = y';

%% ------------------------------------------------------------------------

xc11 = 0.7500;
yc11 = 0.5000;

xc12 = 2.2500;
yc12 = 1.5000;

%% ------------------------------------------------------------------------

xc21 = 1.1250;
yc21 = 0.7500;

xc22 = 1.8750;
yc22 = 1.2500;

%% ------------------------------------------------------------------------

xc31 = 1.3125;
yc31 = 0.8750;

xc32 = 1.6875;
yc32 = 1.1250;

%% ------------------------------------------------------------------------

rc1 = 2.0 / 3 / 2;
rc2 = 1.0 / 3 / 2;

[X,Y] = meshgrid( x, y );

%% ------------------------------------------------------------------------

Z1 = rc1 - sqrt( ( X - xc11 ) .^ 2 + ( Y - yc11 ) .^ 2 );
Z2 = rc2 - sqrt( ( X - xc12 ) .^ 2 + ( Y - yc12 ) .^ 2 );
Z = max( Z1, Z2 );
C = zeros( size( Z, 1 ), size( Z, 2 ) );

% phi
figure 
hold on
h = surf( X, Y, Z );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, Z, [0 0], 'Linecolor', [0 0 0]);
h = surf( X, Y, C );
set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])

%% ------------------------------------------------------------------------

Z1 = rc1 - sqrt( ( X - xc21 ) .^ 2 + ( Y - yc21 ) .^ 2 );
Z2 = rc2 - sqrt( ( X - xc22 ) .^ 2 + ( Y - yc22 ) .^ 2 );
Z = max( Z1, Z2 );
C = zeros( size( Z, 1 ), size( Z, 2 ) );

% phi
figure 
hold on
h = surf( X, Y, Z );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, Z, [0 0], 'Linecolor', [0 0 0]);
h = surf( X, Y, C );
set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])

%% ------------------------------------------------------------------------

Z1 = rc1 - sqrt( ( X - xc31 ) .^ 2 + ( Y - yc31 ) .^ 2 );
Z2 = rc2 - sqrt( ( X - xc32 ) .^ 2 + ( Y - yc32 ) .^ 2 );
Z = max( Z1, Z2 );
C = zeros( size( Z, 1 ), size( Z, 2 ) );

% phi
figure 
hold on
h = surf( X, Y, Z );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, Z, [0 0], 'Linecolor', [0 0 0]);
h = surf( X, Y, C );
set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])
