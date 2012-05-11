function accept = metropolis_hastings( old_ll, new_ll, proposal_bias )

p_accept = min( 1, exp( new_ll - old_ll ) / proposal_bias ) ;
accept   = rand() < p_accept ;

end