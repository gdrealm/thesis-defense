rho = [0.0:0.01:1.0];

beta0 = 0.0;
beta1 = 1.0;
beta2 = 2.0;
beta4 = 4.0;
beta8 = 8.0;
betainf = 1000.0;

rho0 = 1.0 - exp( -beta0 * rho ) + rho * exp( -beta0 );
rho1 = 1.0 - exp( -beta1 * rho ) + rho * exp( -beta1 );
rho2 = 1.0 - exp( -beta2 * rho ) + rho * exp( -beta2 );
rho4 = 1.0 - exp( -beta4 * rho ) + rho * exp( -beta4 );
rho8 = 1.0 - exp( -beta8 * rho ) + rho * exp( -beta8 );
rhoinf = 1.0 - exp( -betainf * rho ) + rho * exp( -betainf );

figure
hold all
plot( rho, rho0, '-x' );
plot( rho, rho1, '-+' );
plot( rho, rho2, '-o' );
plot( rho, rho4, '-s' );
plot( rho, rho8, '-d' );
plot( rho, rhoinf, '-*' );

xlabel( '$\tilde{\rho}$', 'Interpreter', 'Latex', 'FontSize', 20 );
ylabel( '$\hat{\rho}$', 'Interpreter', 'Latex', 'FontSize', 20 );

h =  legend( '$\beta=0.0$', '$\beta=1.0$', '$\beta=2.0$', '$\beta=4.0$', '$\beta=8.0$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );