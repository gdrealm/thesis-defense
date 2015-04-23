clc
close all

pwd = ['/data/villanueva/work/optimization/curvature/paper_study/beam2D_quad4_XFEM_sweep_nnrm_InterfaceMeanCurvature_circle'];

perimeter = [pwd, '/', 'GSOL_INTERFACE_PERIMETER'];
curvature = [pwd, '/', 'GSOL_NORMAL_INTERFACE_MEAN_CURVATURE'];

perimeter = load(perimeter);
curvature = load(curvature);

r_min = 0.25;
r_max = 0.95;
steps = 100;

r = [r_min : (r_max - r_min ) / steps : r_max];
r = r';

exact_perimeter = ones(size(r,1),1);
exact_perimeter = r * 2 * pi;

exact_curvature = ones(size(r,1),1);
exact_curvature = exact_curvature * 2 * pi;

exact_zeros = zeros(size(r,1),1);

% Error curvature with perimeter, and curvature and real.
figure
plot(r,(curvature-exact_curvature)./exact_curvature,'-d',r, (perimeter-exact_perimeter)./exact_perimeter,'-',r, exact_zeros,'-')
axis tight

% Absolute error curvature with perimeter, and curvature and real.
figure
plot(r,abs((curvature-exact_curvature)./exact_curvature),'-d',r, abs((perimeter-exact_perimeter)./exact_perimeter),'-',r, exact_zeros,'-')
axis tight

% Curvature and real
figure
plot(r,curvature,'-o',r,exact_curvature,'-')
axis tight

% Perimeter and real
% figure
% plot(r,perimeter,'-o',r,exact_perimeter,'-')
