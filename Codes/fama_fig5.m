clc
close all
clear all

rgb = @(x,y,z) [x,y,z]/255;
colors = [ rgb(71, 147, 175); rgb(255, 196, 112); rgb(221, 87, 70) ];

%% Parameters
gamma_th_db = [2, 3];
gamma_th = db2pow( gamma_th_db );

sigma_s = 1; % Squared symbol power
% sigma_g = sqrt( db2pow( 20 ) ); % Squared channel power
sigma_n = 1; % Squared noise power
avg_snr = db2pow( 20 );
% input_snr = ( sigma_s^2 * sigma_g^2 ) / ( sigma_n^2 );

% Length factor
w_v = 1 : 10;
% Number of ports
num_ports = 100;
% Fading factor
m = 2;
sigma_g = sqrt( avg_snr / ( 2 * m ) ) * ( sigma_n / sigma_s );

% Number of users
num_users = [3, 4];
% Number of integral samples
num_points = 200;
% Maximum mc samples
max_mc_s = 500000;

% Outage prob
an_out_prob = zeros( length( num_users ), length( w_v ), length( gamma_th ) );
sim_out_prob = zeros( length( num_users ), length( w_v ), length( gamma_th ) );
gl_out_prob = zeros( length( num_users ), length( w_v ), length( gamma_th ) );
n_gl = 10;

step_c = 0;

for g_p = 1 : length( gamma_th )
    
    for w_p = 1 : length( w_v )

        w = w_v( w_p );
        % Correlation factor
        corr_factor = get_corr_factor( num_ports, w );

        for u_p = 1 : length( num_users )
            step_c = step_c + 1;
            % Analytical
            an_out_prob( u_p, w_p, g_p ) = exact_fama_op_nakagami( num_ports, num_users( u_p ), gamma_th( g_p ), corr_factor, sigma_g, sigma_s, sigma_n, m, num_points );
            % Monte carlo
            num_en_th = min( round( 1000 * ( 1 / an_out_prob( u_p, w_p, g_p ) ) ), max_mc_s );
            num_en = round( 1000 * ( 1 / an_out_prob( u_p, w_p, g_p ) ) );
            num_int = ceil( num_en / max_mc_s );
            for int = 1 : num_int

                sim_out_prob( u_p, w_p, g_p ) = sim_out_prob( u_p, w_p, g_p ) + mc_fama_op_nakagami( num_en_th, num_ports, num_users( u_p ), gamma_th( g_p ), corr_factor, sigma_g, sigma_n, m );
            end
            sim_out_prob( u_p, w_p, g_p ) = sim_out_prob( u_p, w_p, g_p ) / num_int;

            num_en = min( round( 1000 * ( 1 / an_out_prob( u_p, w_p ) ) ), max_mc_s );
            sim_out_prob( u_p, w_p, g_p ) = mc_fama_op_nakagami( num_en, num_ports, num_users( u_p ), gamma_th( g_p ), corr_factor, sigma_g, sigma_n, m );


            % Gauss Laguerre
            gl_out_prob( u_p, w_p, g_p ) = gs_fama_op_nakagami( num_ports, num_users( u_p ), gamma_th( g_p ), corr_factor, sigma_g, sigma_s, sigma_n, m, n_gl );

            fprintf( 'Step: (%d/%d), An. OP: %f, MC OP (%d): %f, GL OP: %f, ratio: %f\n', ...
                step_c, length( num_ports ) * length( m ) * length( num_users ), ...
                an_out_prob( u_p, w_p, g_p ), num_en, sim_out_prob( u_p, w_p, g_p ), ...
                gl_out_prob( u_p, w_p, g_p ), ...
                an_out_prob( u_p, w_p, g_p ) / sim_out_prob( u_p, w_p, g_p ) );
        end

    end
    figure( 1 )
    he(1, g_p) = semilogy( w_v, an_out_prob( 1, :, g_p ), 'Color', colors( 1, : ), 'Linewidth', 2 );
    hold on
    he(2, g_p) = semilogy( w_v, an_out_prob( 2, :, g_p ), 'Color', colors( 2, : ), 'Linewidth', 2 );

    hs(1, g_p) = semilogy( w_v, sim_out_prob( 1, :, g_p ), '+', 'Color', 'k', 'Linewidth', 2 );
    hs(2, g_p) = semilogy( w_v, sim_out_prob( 2, :, g_p ), '+', 'Color', 'k', 'Linewidth', 2 );

    hgl(1, g_p) = semilogy( w_v, gl_out_prob( 1, :, g_p ), 'square', 'Color', 'r', 'Linewidth', 1 );
    hgl(2, g_p) = semilogy( w_v, gl_out_prob( 2, :, g_p ), 'square', 'Color', 'r', 'Linewidth', 1 );

end

axx = gca;
axx.TickLabelInterpreter = 'latex';
axx.FontSize = 15;

xlabel( 'Antenna Size -- $W$', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 15 );
ylim([1e-4, 1]);
xlim([1, 10]);

grid on
legend( [he(1, 1), he(2, 2), hs(1, 1), hgl(1, 1)], {'$U = 3$', '$U = 4$', 'Simulation', 'Gauss-Laguerre'}, 'Interpreter', 'Latex', 'FontSize', 14 );

