clc
close all

pwd = ['.'];

curvature = [pwd, '/', 'GSOL_INTERFACE_CURVATURE'];
curvature = load(curvature);

r_min = 0.25;
r_max = 0.95;
steps = 100;

r = [r_min : (r_max - r_min ) / steps : r_max];
r = r';

exact_curvature = ones(size(r,1),1);
exact_curvature = exact_curvature * 2 * pi ./ r;

exact_zeros = zeros(size(r,1),1);

% Error curvature with perimeter, and curvature and real.
figure
% plot(r,abs(curvature-exact_curvature)./exact_curvature,'-d',r, abs(perimeter-exact_perimeter)./exact_perimeter,'-',r, exact_zeros,'-')
plot(r,abs(curvature-exact_curvature),'-d')
axis tight
xlabel('$r_{c}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\mathrm{Error}$', 'Interpreter', 'Latex', 'FontSize', 20);

h =  legend( '$\vert \int_{\Gamma_{\phi=0}} \kappa_{n}\,\mathrm{d}\Gamma - 2\pi \vert$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );

% Curvature and real
figure
plot(r,curvature,'-o',r,exact_curvature,'-')
axis tight
xlabel('$r_{c}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\kappa_{n}$', 'Interpreter', 'Latex', 'FontSize', 20);
h =  legend( '$\kappa_{n}$',  '$2\pi$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );

