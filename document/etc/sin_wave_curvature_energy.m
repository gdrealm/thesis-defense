clc
clear all
close all

A_min = 0.25;
A_max = 0.95;
steps = 100;

A = [A_min : (A_max - A_min ) / steps : A_max];
A = A';

clc
close all

pwd = ['C:\Users\Hernan\Desktop\comps\src\topology_optimization_approaches_for_the_curvature_minimization_of_level_set_isocontours\etc\beam2D_quad4_XFEM_sweep_nnrm_InterfaceCurvatureSquared_sin'];
curvature = [pwd, '/', 'GSOL_NORMAL_INTERFACE_CURVATURE_SQUARED'];
curvature = load(curvature);

kappa = zeros( size( A, 1 ), 1 );

for i = 1 : size( A, 1 )
    a = A( i );
    numerator = ( a ^ 2 ) * ( pi ^ 4 ) * ( 5 + 3 * a * pi ) * ( 5 + ( ( 1 / ( 1 + a * pi ) ) ^ ( 1 / 2 ) ) * ( ( 1 + a * pi ) ^ ( 1 / 2 ) ) );
    denominator = 16 * ( 1 + a * pi ) ^ ( 3 / 2 );
    kappa( i ) = numerator / denominator;
end

hold all
plot( A, kappa );
plot( A, curvature );