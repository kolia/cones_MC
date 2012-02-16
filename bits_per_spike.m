function ll = bits_per_spike( ll, cone_map )

ll = ll / ( sum( cone_map.N_spikes ) * log(2) ) ;

end