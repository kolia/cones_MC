function sta = normalize_sta( sta )

sta = sta-min(sta(:)) ;
sta = sta / max(sta(:)) ;

end