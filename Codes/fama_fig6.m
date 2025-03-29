clc
close all
clear all

rgb = @(x,y,z) [x,y,z]/255;
colors = [ rgb(71, 147, 175); rgb(255, 196, 112); rgb(221, 87, 70) ];
trac = ['-s';'-x';'-o'];
%% Parameters
gamma_th_db = 6;
gamma_th = db2pow( gamma_th_db );

sigma_s = 1; % Squared symbol power
% sigma_g = sqrt( db2pow( 20 ) ); % Squared channel power
sigma_n = 1; % Squared noise power
avg_snr = db2pow( 20 );
% input_snr = ( sigma_s^2 * sigma_g^2 ) / ( sigma_n^2 );

% Length factor
w_v = 2;
% Number of ports
num_ports = [ 200, 300, 400 ];
% Fading factor
m = [1,2,3];

% Number of users
num_users = 2 : 8;
% Number of integral samples
num_points = 200;
% Maximum mc samples
max_mc_s = 500000;

% Correlation factor
w = 2;


upp_bound = zeros( length( num_users ), length( num_ports ), length( m ) );
an_out_prob = zeros( length( gamma_th ), length( num_users ) );
mult_gain = zeros( length( num_users ), length( num_ports ), length( m ) );

step_c = 1;

for n = 1 : length( num_ports )
    
    corr_factor = get_corr_factor( num_ports( n ), w );
    
    for mm = 1 : length( m )
        
        for u = 1 : length( num_users )

            an_out_prob = exact_fama_op_sir_nakagami( num_ports( n ), num_users( u ), gamma_th, corr_factor, m( mm ), 200 );
            mult_gain( u, n, mm ) = num_users( u ) * ( 1 - an_out_prob );    
            
            fprintf( 'Step: (%d/%d), An. OP: %f, GL OP: %f, Mult. Gain: %f\n', ...
                step_c, length( num_ports ) * length( m ) * length( num_users ), ...
                an_out_prob,...
                0, ...
                mult_gain( u, n, mm ) );
            step_c = step_c + 1;
        end
        
        u_opt = opt_num_users( num_ports( n ), gamma_th, w, m( mm ) );
        h(n,mm) = plot( num_users, mult_gain( :, n, mm ), '-s', 'Color', colors( n, : ), 'Linewidth', 1.2 );
        hold on
        hs = plot( u_opt, mult_gain( u_opt - 2 + 1, n, mm ), 'p', 'MarkerSize',10, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6] );
        
    end
    
end
axx = gca;
axx.TickLabelInterpreter = 'latex';
axx.FontSize = 15;
grid on
xlabel( 'Number of Users -- $U$', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Multiplexing gain -- $m_{g}$', 'Interpreter', 'Latex', 'FontSize', 15 );
xlim([2,8]);
legend( [h(1,1), h(2,2), h(3,3), hs], {'$N=200$', '$N=300$', '$N=400$', '$m_{g}(U^{*})$'}, 'Interpreter', 'Latex', 'FontSize', 14 );

