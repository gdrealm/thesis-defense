gamma = [0.0:0.01:1.0];

% alpha_pen0 = 0;
alpha_pen1 = 10^(+0);
alpha_pen2 = 10^(-1);
alpha_pen3 = 10^(-2);
alpha_pen4 = 10^(-4);

alpha_max = 10^4;
alpha_min = 0;

% alpha0 = alpha_max + ( alpha_min - alpha_max ) .* gamma .* ( 1 + alpha_pen0 ) ./ ( gamma + alpha_pen0 );
alpha1 = alpha_max + ( alpha_min - alpha_max ) .* gamma .* ( 1 + alpha_pen1 ) ./ ( gamma + alpha_pen1 );
alpha2 = alpha_max + ( alpha_min - alpha_max ) .* gamma .* ( 1 + alpha_pen2 ) ./ ( gamma + alpha_pen2 );
alpha3 = alpha_max + ( alpha_min - alpha_max ) .* gamma .* ( 1 + alpha_pen3 ) ./ ( gamma + alpha_pen3 );
alpha4 = alpha_max + ( alpha_min - alpha_max ) .* gamma .* ( 1 + alpha_pen4 ) ./ ( gamma + alpha_pen4 );

figure
hold all
% plot( gamma, alpha0, '-x' );
plot( gamma, alpha1, '-+' );
plot( gamma, alpha2, '-o' );
plot( gamma, alpha3, '-s' );
plot( gamma, alpha4, '-d' );

xlabel( '$\gamma$', 'Interpreter', 'Latex', 'FontSize', 20 );
ylabel( '$\alpha$', 'Interpreter', 'Latex', 'FontSize', 20 );

h =  legend( '$\alpha_{p}=10^{+0}$', '$\alpha_{p}=10^{-1}$', '$\alpha_{p}=10^{-2}$', '$\alpha_{p}=10^{-4}$' );
set( h, 'Interpreter', 'Latex', 'FontSize', 20 );