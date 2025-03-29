function mc_out_prob = mc_fama_op_nakagami( num_en, num_ports, num_users, gamma_th, corr_factor, sigma_g, sigma_n, m_nkg )

    mc_out_prob = 0;
    % Generate channels for each UE
    G_u = zeros( num_users, num_en, num_ports );
    % Interfenrece for each UE
    I_u = zeros( num_users, num_en, num_ports );
    % SINR for each UE
    SINR_u = zeros( num_users, num_en, num_ports );
    for u = 1 : num_users
       
        G_u( u, :, : ) = ch_gain_sp_nkg( sqrt( 1 / m_nkg ) * sigma_g, sqrt( corr_factor ), m_nkg, num_ports, num_en );
    end
    
    for u = 1 : num_users
        for u_p = 1 : num_users
            if( u_p ~= u )
                % Interference plus noise between the u_p x u
                I_u( u, :, : ) = I_u( u, :, : ) + abs( G_u( u_p, :, : ) ).^2;
            end
        end
        SINR_u( u, :, : ) = abs( G_u( u, :, : ) ).^2 ./ ( I_u( u, :, : ) + sigma_n^2 );
    end
    
    for u = 1 : num_users
        local_sinr = reshape( SINR_u( u, :, : ), num_en, num_ports );
        mc_out_prob = mc_out_prob + sum( max( local_sinr, [], 2 ) < gamma_th, 'all' ) / ( num_en * num_users );
    end


end

function h_ch = ch_gain_sp_nkg( sigma, corr_factor, m_nkg, num_ports, num_curves )

    % Alloc
    g_ch = zeros( num_curves, num_ports, m_nkg );
    h_ch = zeros( num_curves, num_ports );
    % Gaussian matrix

    % Correlation factor
    mu_k = corr_factor;
    comp_mu_k = sqrt( 1 - mu_k^2 );
    % Channels for 1 to num_ports
    X_0 = normrnd( 0, 1/sqrt( 2 ), [ num_curves, 1, m_nkg ] );
    Y_0 = normrnd( 0, 1/sqrt( 2 ), [ num_curves, 1, m_nkg ] );
    for k = 1 : num_ports

        X_k = normrnd( 0, 1/sqrt( 2 ), [ num_curves, 1, m_nkg ] );
        Y_k = normrnd( 0, 1/sqrt( 2 ), [ num_curves, 1, m_nkg ] );
        % Channel gain
        g_ch( :, k, : ) = sigma * ( ( comp_mu_k * X_k + mu_k * X_0 ) + 1i * ( comp_mu_k * Y_k + mu_k * Y_0 ) );
    end
    
    for m = 1 : m_nkg
        
        h_ch = h_ch + abs( g_ch( :, :, m ) ).^2;        
    end
    h_ch = sqrt( h_ch );
    
end