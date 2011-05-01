function default( s , name , default_value )
% Set value of argout to s.(name), unless s doesn't have a field called
% 'name', in which case set it to default_value.  Uses the
% matlab 'assignin' function to assign the output value to argout.

if isfield( s , name )
    assignin( 'caller' , name , s.(name) ) ;
else
    assignin( 'caller' , name , default_value ) ;
end

end