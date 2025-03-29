clc
close all
clear all

rgb = @(x,y,z) [x,y,z]/255;
colors = [ rgb(71, 147, 175); rgb(255, 196, 112); rgb(221, 87, 70) ];

%% Parameters
gamma_th_db = 1;
gamma_th = db2pow( gamma_th_db );

sigma_s = 1; % Squared symbol power
sigma_g = 1; % Squared channel power
avg_snr = db2pow( linspace( -10, 25, 20 ) );
% Length factor
w = 1;
% Number of ports
num_ports = 150;
% Correlation factor
corr_factor = get_corr_factor( num_ports, w ); % Correlation factor    

% Fading factor
m = [1,2,3];
% Number of users
num_users = 4;
% Number of integral samples
num_points = 150;
% Maximum mc samples
max_mc_s = 500000;

% Outage prob
gl_out_prob = zeros( length( m ), length( avg_snr ) );
an_out_prob = zeros( length( m ), length( avg_snr ) );
sim_out_prob = zeros( length( m ), length( avg_snr ) );

step_c = 0;

for m_p = 1 : length( m )
    
    % SIR OP
    out_prob_sir = exact_fama_op_sir_nakagami( num_ports, num_users, gamma_th, corr_factor, m( m_p ), num_points );
    
    for g_p = 1 : length( avg_snr )

        sigma_n = sqrt( 2 * m( m_p ) / avg_snr( g_p ) ) * sigma_g * sigma_s;
        
        step_c = step_c + 1;
        % Analytical
        an_out_prob( m_p, g_p ) = exact_fama_op_nakagami( num_ports, num_users, gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m( m_p ), num_points );
        
        % Monte carlo
        num_en_th = min( round( 1000 * ( 1 / an_out_prob( m_p, g_p ) ) ), max_mc_s );
        num_en = round( 1000 * ( 1 / an_out_prob( m_p, g_p ) ) );
        num_int = ceil( num_en / max_mc_s );
        for int = 1 : num_int

            sim_out_prob( m_p, g_p ) = sim_out_prob( m_p, g_p ) + mc_fama_op_nakagami( num_en_th, num_ports, num_users, gamma_th, corr_factor, sigma_g, sigma_n, m( m_p ) );
        end
        sim_out_prob( m_p, g_p ) = sim_out_prob( m_p, g_p ) / num_int;
        
        % Gauss-Laguerre
        gl_out_prob( m_p, g_p ) = gs_fama_op_nakagami( num_ports, num_users, gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m( m_p ), 10 );
        
        fprintf( 'Step: (%d/%d), An. OP: %f, MC OP: %f, SIR OP: %f, GL OP: %f, ratio: %f\n', ...
            step_c, length( sigma_g ) * length( m ), ...
            an_out_prob( m_p, g_p ), sim_out_prob( m_p, g_p ), out_prob_sir, gl_out_prob( m_p, g_p ), ...
            an_out_prob( m_p, g_p ) / sim_out_prob( m_p, g_p ) );

    end

    h(m_p, 1) = semilogy( 10 * log10( avg_snr ), an_out_prob( m_p, : ), 'Color', colors( 1, : ), 'Linewidth', 2 );
    hold on
    h(m_p, 2) = semilogy( 10 * log10( avg_snr ), sim_out_prob( m_p, : ), '+', 'Color', 'k', 'Linewidth', 2 );
    h(m_p, 3) = semilogy( 10 * log10( avg_snr ), gl_out_prob( m_p, : ), 's', 'Color', 'r', 'Linewidth', 1 );
    h(m_p, 4) = semilogy( [ -10, 25 ], [ out_prob_sir, out_prob_sir ], '--', 'Color', 'r', 'Linewidth', 2 );
    
    
end

xlim( [-10, 25] );
ylim( [1e-3, 1] );

axx = gca;
axx.TickLabelInterpreter = 'latex';
axx.FontSize = 15;

xlabel( 'Average SNR (dB)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 15 );
grid on

legend( [h(1, 4), h(1, 1), h(1, 2), h(1, 3)], {'(10)', '(11)', 'Simulation', 'Gauss-Laguerre'}, 'Interpreter', 'Latex', 'FontSize', 14 );