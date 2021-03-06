clc
close all

pwd = ['.'];

curvature = [pwd, '/', 'GSOL_NORMAL_INTERFACE_CURVATURE_SQUARED'];
curvature = load(curvature);

r_min = 0.25;
r_max = 0.95;
steps = 100;

r = [r_min : (r_max - r_min ) / steps : r_max];
r = r';

exact_curvature = ones(size(r,1),1);
exact_curvature = exact_curvature * 2 * pi ./ r;

exact_zeros = zeros(size(r,1),1);

% Curvature and real
figure
plot(r,curvature./2,'-o',r,exact_curvature,'-')
axis tight
xlabel('$r_{c}$', 'Interpreter', 'Latex', 'FontSize', 20);
% ylabel('h', 'Interpreter', 'Latex', 'FontSize', 20);
h =  legend( '$\int_{\phi=0} \Vert \frac{\mathrm{d}\mathbf{n}_{u}}{\mathrm{d}s} \Vert^2 \,\mathrm{d}\Gamma$',  '$\frac{2\pi}{r_{c}}$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );


% Error curvature with perimeter, and curvature and real.
figure
% plot(r,abs(curvature-exact_curvature)./exact_curvature,'-d',r, abs(perimeter-exact_perimeter)./exact_perimeter,'-',r, exact_zeros,'-')
plot(r,abs(curvature./2-exact_curvature),'-d')
axis tight
xlabel('$r_{c}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\mathrm{Error}$', 'Interpreter', 'Latex', 'FontSize', 20);

h =  legend( '$\vert\int_{\phi=0} \Vert \frac{\mathrm{d}\mathbf{n}_{u}}{\mathrm{d}s} \Vert^2 \,\mathrm{d}\Gamma - \frac{2\pi}{r_{c}} \vert$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );

