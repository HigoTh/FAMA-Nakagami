clc
close all
clear all

rgb = @(x,y,z) [x,y,z]/255;
colors = [ rgb(71, 147, 175); rgb(255, 196, 112); rgb(221, 87, 70) ];

%% Parameters
gamma_th_db = 1;
gamma_th = db2pow( gamma_th_db );

sigma_s = 1; % Squared symbol power
% sigma_g = sqrt( db2pow( 20 ) ); % Squared channel power
sigma_n = 1; % Squared noise power
avg_snr = db2pow( 20 );

% Length factor
w = 1;

% Number of ports
num_ports = round( logspace( 1, log10( 150 ), 10 ) );
% Fading factor
m = [1,2,3];
% Number of users
num_users = [3,5];
% Number of integral samples
num_points = 100;
% Maximum mc samples
max_mc_s = 1000000;

% Outage prob
an_out_prob = zeros( length( num_ports ), length( m ), length( num_users ) );
sim_out_prob = zeros( length( num_ports ), length( m ), length( num_users ) );
gl_out_prob = zeros( length( num_ports ), length( m ), length( num_users ) );
n_gl = 10;

step_c = 0;
for u_p = 1 : length( num_users )

    for m_p = 1 : length( m )
    
        sigma_g = sqrt( avg_snr / ( 2 * m( m_p ) ) ) * ( sigma_n / sigma_s );
        
        for n = 1 : length( num_ports )
            
            corr_factor = get_corr_factor( num_ports( n ), w );
            
            step_c = step_c + 1;
            % Analytical
            an_out_prob( n, m_p, u_p ) = exact_fama_op_nakagami( num_ports( n ), num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m( m_p ), num_points );
            % Monte carlo
            num_en = round( 100 * ( 1 / an_out_prob( n, m_p, u_p ) ) );
            num_int = ceil( num_en / max_mc_s );
            for int = 1 : num_int
                
                sim_out_prob( n, m_p, u_p ) = sim_out_prob( n, m_p, u_p ) + mc_fama_op_nakagami( max_mc_s, num_ports( n ), num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_n, m( m_p ) );
            end
            sim_out_prob( n, m_p, u_p ) = sim_out_prob( n, m_p, u_p ) / num_int;

            % Gauss Laguerre
            gl_out_prob( n, m_p, u_p ) = gs_fama_op_nakagami( num_ports( n ), num_users( u_p ), gamma_th, corr_factor, sigma_g, sigma_s, sigma_n, m( m_p ), n_gl );
            
            fprintf( 'Step: (%d/%d), An. OP: %f, MC OP (%d): %f, GL OP: %f, ratio: %f\n', ...
                step_c, length( num_ports ) * length( m ) * length( num_users ), ...
                an_out_prob( n, m_p, u_p ), 0, sim_out_prob( n, m_p, u_p ), ...
                gl_out_prob( n, m_p, u_p ), ...
                an_out_prob( n, m_p, u_p ) / sim_out_prob( n, m_p, u_p ) );
        end
        figure( 1 )
        he(u_p, m_p) = loglog( num_ports, an_out_prob( :, m_p, u_p ), 'Color', colors( m_p, : ), 'Linewidth', 2 );
        hold on
        hs(u_p, m_p) = loglog( num_ports, sim_out_prob( :, m_p, u_p ), '+', 'Color', 'k', 'Linewidth', 2 );
        hgl(u_p, m_p) = loglog( num_ports, gl_out_prob( :, m_p, u_p ), 'square', 'Color', 'r', 'Linewidth', 1 );
        axx = gca;
        axx.TickLabelInterpreter = 'latex';
        axx.FontSize = 15;

        xlabel( 'Number of Ports -- $N$', 'Interpreter', 'Latex', 'FontSize', 15 );
        ylabel( 'Outage Probability', 'Interpreter', 'Latex', 'FontSize', 15 );
        
    end
end
grid on
legend( [he(1, 1), he(1, 2), he(1, 3), hs(1, 1), hgl(1, 1)], {'$m=1$', '$m=2$', '$m=3$', 'Simulation', 'Gauss-Laguerre'}, 'Interpreter', 'Latex', 'FontSize', 14 );


