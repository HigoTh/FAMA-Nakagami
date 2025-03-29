clc
close all
clear all

rgb = @(x,y,z) [x,y,z]/255;
colors = [ rgb(238, 153, 102); rgb(51, 115, 230); rgb(51, 115, 87) ];

%% Parameters
gamma_th_db = 0;
gamma_th = db2pow( gamma_th_db );

sigma_s = 1; % Squared symbol power
sigma_g = 1; % Squared channel power
avg_snr = db2pow( linspace( -10, 25, 30 ) );
avg_snr_sim = db2pow( linspace( -10, 25, 5 ) );
% Length factor
w = 1;
% Number of ports
num_ports_v = [100, 150];
   

% Fading factor
m = 2;
% Number of users
num_users = [3,4,5];
% Number of integral samples
num_points = 100;
% Maximum mc samples
max_mc_s = 100000;

% Outage prob
gl_out_prob = zeros( length( num_ports_v ) * length( num_users ), length( avg_snr_sim ) );
an_out_prob = zeros( length( num_ports_v ) * length( num_users ), length( avg_snr ) );
sim_out_prob = zeros( length( num_ports_v ) * length( num_users ), length( avg_snr_sim ) );

step_c = 0;

id = 1;
for n_p = 1 : length( num_ports_v )
    
    num_ports = num_ports_v( n_p );
    % Correlation factor
    corr_factor = get_corr_factor( num_ports, w ); % Correlation factor 
    
    for u_p = 1 : length( num_users )

        % SIR OP
        out_prob_sir = exact_fama_op_sir_nakagami( num_ports, num_users( u_p ), gamma_th, corr_factor, m, num_points );
%         out_prob_sir = 0;
        
        for g_p = 1 : length( avg_snr )

            sigma_n = sqrt( 2 * m / avg_snr( g_p ) ) * sigma_g * sigma_s;

            % Analytical
            an_out_prob( id, g_p ) = exact_fama_op_nakagami( num_ports, num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m, num_points );
            step_c = step_c + 1;
            fprintf( 'Step: %d\n', step_c );
            
        end

        for g_l = 1 : length( avg_snr_sim )
            
            sigma_n = sqrt( 2 * m / avg_snr_sim( g_l ) ) * sigma_g * sigma_s;
            
            % Gauss-Laguerre
            gl_out_prob( id, g_l ) = gs_fama_op_nakagami( num_ports, num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m, 20 );

            % Monte carlo
            num_en_th = min( round( 200 * ( 1 / gl_out_prob( id, g_l ) ) ), max_mc_s );
            num_en = round( 200 * ( 1 / gl_out_prob( id, g_l ) ) );
            num_int = ceil( num_en / max_mc_s );
            for int = 1 : num_int

                sim_out_prob( id, g_l ) = sim_out_prob( id, g_l ) + mc_fama_op_nakagami( num_en_th, num_ports, num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_n, m );
            end
            sim_out_prob( id, g_l ) = sim_out_prob( id, g_l ) / num_int;
            
            step_c = step_c + 1;
            fprintf( 'Step: %d\n', step_c );
            
        end
        
        h1(id, n_p) = semilogy( 10 * log10( avg_snr ), an_out_prob( id, : ), 'Color', colors( n_p, : ), 'Linewidth', 2 );
        hold on
        h2(id, n_p) = semilogy( 10 * log10( avg_snr_sim ), sim_out_prob( id, : ), '+', 'Color', 'k', 'Linewidth', 2 );
        h3(id, n_p) = semilogy( 10 * log10( avg_snr_sim ), gl_out_prob( id, : ), 's', 'Color', 'r', 'Linewidth', 1 );
        h4(id, n_p) = semilogy( [ -10, 25 ], [ out_prob_sir, out_prob_sir ], '--', 'Color', 'r', 'Linewidth', 1 );

        id = id + 1;
    end
end

xlim( [-10, 25] );
ylim( [1e-6, 1] );

axx = gca;
axx.TickLabelInterpreter = 'latex';
axx.FontSize = 15;

xlabel( 'Average SNR (dB)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 15 );
grid on

legend( [h4(1, 1), h1(2, 1), h1(4, 2), h3(1, 1), h2(1, 1)], {'(10)', '(11)', 'x', 'Simulation', 'Gauss-Laguerre'}, 'Interpreter', 'Latex', 'FontSize', 14, 'NumColumns', 2 );


