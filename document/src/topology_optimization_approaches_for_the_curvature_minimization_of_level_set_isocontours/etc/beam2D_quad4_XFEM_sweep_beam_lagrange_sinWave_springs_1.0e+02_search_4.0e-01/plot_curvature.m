clc
close all

pwd1 = ['/data/villanueva/work/optimization/curvature/beam2D_structuralCompliance_quad4_XFEM_Beam_Lagrange_sin_springs_1.0e+02'];
pwd2 = ['/data/villanueva/work/optimization/curvature/beam2D_structuralCompliance_quad4_XFEM_Beam_Lagrange_sin_springs_1.0e+00'];
pwd3 = ['/data/villanueva/work/optimization/curvature/beam2D_structuralCompliance_quad4_XFEM_Beam_Lagrange_sin_springs_1.0e-02'];
pwd4 = ['/data/villanueva/work/optimization/curvature/beam2D_structuralCompliance_quad4_XFEM_Beam_Lagrange_sin_springs_1.0e-04'];

% perimeter1 = [pwd1, '/', 'GSOL_INTERFACE_PERIMETER'];
curvature1 = [pwd1, '/', 'GSOL_INTERFACE_CURVATURE_BEAM_LAGRANGE'];
curvature2 = [pwd2, '/', 'GSOL_INTERFACE_CURVATURE_BEAM_LAGRANGE'];
curvature3 = [pwd3, '/', 'GSOL_INTERFACE_CURVATURE_BEAM_LAGRANGE'];
curvature4 = [pwd4, '/', 'GSOL_INTERFACE_CURVATURE_BEAM_LAGRANGE'];

% perimeter1 = load(perimeter1);
curvature1 = load(curvature1);
curvature2 = load(curvature2);
curvature3 = load(curvature3);
curvature4 = load(curvature4);

r_min = 0.055;
r_max = 1.375;
steps = 100;

r = [r_min : (r_max - r_min ) / steps : r_max];
r = r';

% Error curvature with perimeter, and curvature and real.
figure
plot(r,curvature1,'-d',r,curvature2,'-s',r,curvature3,'-o',r,curvature4,'-+')
axis tight

% Curvature.
figure
plot(r,curvature1,'-o')
axis tight

figure
plot(r,curvature2,'-o')
axis tight

figure
plot(r,curvature3,'-o')
axis tight

figure
plot(r,curvature4,'-o')
axis tight