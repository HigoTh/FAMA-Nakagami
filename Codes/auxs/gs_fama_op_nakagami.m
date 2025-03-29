function gl_out_prob = gs_fama_op_nakagami( num_ports, num_users, gamma_th, delta, sigma_g, sigma_s, sigma_n, m_nkg, n_gl )

    % Gaussa laguerre
	[a_i, w_i] = GaussLaguerre( n_gl, 0 );
    a_j = a_i;
    a_l = a_i;
    w_j = w_i;
    w_l = w_i;
    
    % Constants
    d_f = delta / ( 1 - delta );
    gamma_avg = ( sigma_s^2 * sigma_g^2 ) / ( sigma_n^2 );
    u_eq = (num_users - 1) * m_nkg;
    k1 = 1 / ( gamma( m_nkg ) * gamma( u_eq ) );
    
    fr2_aux = zeros( n_gl, 1 ); 
    for l = 1 : n_gl

        % Marcum-Q Function
        arg_a = sqrt( 2 * d_f * a_l( l ) );
        arg_b = sqrt( gamma_th * a_i + ( ( 2 * gamma_th * m_nkg ) / ( gamma_avg * ( 1 - delta ) ) ) );
        marc_t = marcumq( arg_a, arg_b, m_nkg );
        
        fr1_aux = zeros( n_gl, 1 ); 
        for j = 1 : n_gl
           
            % Auxiliary functions
            f_a = ( a_i / ( 2 * d_f.* a_j( j ) ) ).^( 0.5 * ( u_eq - 1 ) );
            f_b = exp( 0.5 * ( a_i - ( 2 * d_f.* a_j( j ) ) ) );
            f_c = besseli( u_eq - 1, sqrt( 2 * d_f * a_i * a_j( j ) ) );
            f_total = f_a .* f_b .* f_c .* marc_t;
            
            % Sum Marcum term
            int_m_t = ( 1 - 0.5 * sum( f_total .* w_i ) ).^( num_ports );
            % Product with exponential term
            fr1_aux( j, 1 ) = int_m_t * a_j( j )^( u_eq - 1 ) * w_j( j );
        end
        % Integral term
        Fr1 = sum( fr1_aux );
        fr2_aux( l, 1 ) = Fr1 * a_l( l )^( m_nkg - 1 ) * w_l( l );
    end
    gl_out_prob = k1 * sum( fr2_aux );
    
end

function [x, w] = GaussLaguerre(n, alpha)
% This function determines the abscisas (x) and weights (w) for the
% Gauss-Laguerre quadrature of order n>1, on the interval [0, +infinity].
    % Unlike the function 'GaussLaguerre', this function is valid for
    % n>=34. This is due to the fact that the companion matrix (of the n'th
    % degree Laguerre polynomial) is now constructed as a symmetrical
    % matrix, guaranteeing that all the eigenvalues (roots) will be real.
    
    
% © Geert Van Damme
% geert@vandamme-iliano.be
% February 21, 2010    
% Building the companion matrix CM
    % CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial
    % under consideration. Moreover, CM will be constructed in such a way
    % that it is symmetrical.
i   = 1:n;
a   = (2*i-1) + alpha;
b   = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
CM  = diag(a) + diag(b,1) + diag(b,-1);
% Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
[V L]   = eig(CM);
[x ind] = sort(diag(L));
V       = V(:,ind)';
w       = gamma(alpha+1) .* V(:,1).^2;

end
