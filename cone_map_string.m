function base_str = cone_map_string(cone_map)

base_str  = sprintf('%s__%.0eiters', cone_map.datafolder,cone_map.N_iterations) ;

base_str(base_str=='+') = '' ;
base_str(base_str=='.') = '' ;
end