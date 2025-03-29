function out_prob = exact_fama_op_nakagami( num_ports, num_users, gamma_th, delta, sigma_g, sigma_s, sigma_n, m_nkg, num_points )

    % Average SNR
    gamma_avg = ( sigma_s^2 * sigma_g^2 ) / ( sigma_n^2 );
    
    % r-domains
    r1 = linspace( 1e-3, 10 * num_users, num_points ); % r~ (Não pode ser nulo - Divisão por zero)
    r2 = linspace( 1e-3, 10 * num_users, num_points ); % r (Não pode ser nulo - Divisão por zero)

    % Equivalent number of users
    u_eq = (num_users - 1) * m_nkg;
    
    % Constants
    k1 = 1 / ( 2^( u_eq + m_nkg  ) * gamma( m_nkg ) * gamma( u_eq ) );

    fr2_aux = zeros( length( r2 ), 1 );
    for i = 1 : length( r2 )
        
        arg_a = sqrt( ( delta.* r2( i ) )./( 1-delta ) );
        % Find the maximum yk
        y_max = find_maxyk( arg_a, gamma_th, gamma_avg, delta, m_nkg, 1e-5 );
        % yk domain
        yk = linspace( 0, y_max, num_points );
        % Marcum-Q function
        arg_b = sqrt( ( gamma_th * yk ) + ( ( ( 2 * gamma_th * m_nkg ) / ( gamma_avg * ( 1 - delta ) ) ) ) );
        marc_t = marcumq( arg_a, arg_b, m_nkg );
        
        fr1_aux = zeros( length( r1 ), 1 ); 
        for j = 1 : length( r1 )
           
            % Auxiliary functions
            f_a = ( yk / ( ( delta.* r1( j ) )./( 1 - delta ) ) ).^( 0.5 * ( u_eq - 1 ) );
            f_b = exp( -0.5 * ( yk + ( ( delta.* r1( j ) )./( 1 - delta ) ) ) );
            f_c = besseli( u_eq - 1, sqrt( yk * ( delta.* r1( j ) ) / ( 1 - delta ) ) );
            f_total = f_a .* f_b .* f_c .* marc_t;
            
            % Integral Marcum term
            int_m_t = ( 1 - 0.5 * trapz( yk, f_total ) )^( num_ports );
            % Product with exponential term
            fr1_aux( j, 1 ) = int_m_t * exp( -r1( j ) / 2 ) .* r1( j )^( u_eq - 1 );
        end
        % Integral term
        Fr1 = trapz( r1, fr1_aux );
        fr2_aux( i, 1 ) = Fr1 * exp( -r2( i ) / 2 ) * r2( i )^( m_nkg - 1 );
    end
    out_prob = k1 * trapz( r2, fr2_aux );

end

function y_max = find_maxyk( arg_a, gamma_th, gamma_avg, delta, m_nkg, gl_th )

    % Test range
    yk = linspace( 0, 1000, 20 );
    arg_b = sqrt( ( gamma_th * yk ) + ( ( gamma_th / ( gamma_avg * ( 1 - delta ) ) ) ) );
    % Marcum function
    f = marcumq( arg_a, arg_b, m_nkg );
    % Threshold
    f_th = gl_th * max( f );
    % Find the index
    [ ~, i_th ] = min( abs( f - f_th ) );
    y_max = yk( i_th );

end