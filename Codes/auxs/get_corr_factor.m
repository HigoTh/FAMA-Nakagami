function corr_factor = get_corr_factor( num_ports, w )

    k = 1 : num_ports - 1;
    c1 = 2 / ( num_ports * ( num_ports - 1 ) );
    
    corr_factor = abs( c1 * sum( ( num_ports - k ) .* besselj( 0, ( 2 * pi * k * w ) / ( num_ports - 1 ) ) ) );

end