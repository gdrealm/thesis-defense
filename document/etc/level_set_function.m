clc
clear all
close all

x = 0:0.05:3;
y = 0:0.05:2;

x = x';
y = y';

xc = 1.5;
yc = 1.0;

rc = 0.5;

e1 = 1.0;
e2 = 2.0;
e4 = 4.0;

[X,Y] = meshgrid( x, y );
Z = rc - sqrt( ( X - xc ) .^ 2 + ( Y - yc ) .^ 2 );
C = zeros( size( Z, 1 ), size( Z, 2 ) );

D = zeros( size( Z, 1 ), size( Z, 2 ) ) + 0.5;

psi1 = 0.5 * ( tanh( Z ./ ( 2 * e1 ) ) + 1 );
psi2 = 0.5 * ( tanh( Z ./ ( 2 * e2 ) ) + 1 );
psi4 = 0.5 * ( tanh( Z ./ ( 2 * e4 ) ) + 1 );

% phi
figure 
hold on
h = surf( X, Y, Z );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, Z, [0 0], 'Linecolor', [0 0 0]);
% h = surf( X, Y, C );
% set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
axis([0 3 0 2 -1.5 0.7])

% epsilon = 1.0
figure
hold on
h = surf( X, Y, psi1 );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, psi1, [0.5 0.5], 'Linecolor', [0 0 0]);
% h = surf( X, Y, D );
% set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])

% epsilon = 2.0
figure
hold on
h = surf( X, Y, psi2 );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, psi2, [0.5 0.5], 'Linecolor', [0 0 0]);
% h = surf( X, Y, D );
% set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])

% epsilon = 4.0
figure
hold on
h = surf( X, Y, psi4 );
set(h, 'EdgeColor', 'none');
h = contour(X, Y, psi4, [0.5 0.5], 'Linecolor', [0 0 0]);
% h = surf( X, Y, D );
% set(h, 'EdgeColor', 'none');
view(30,12)
grid
axis square
% axis([0 3 0 2 -1.5 0.7])

% h =  legend( '$\phi\left(\mathbf{x}\right)$', '$\phi=0$' );
% set( h, 'Interpreter', 'Latex', 'FontSize', 20 );
