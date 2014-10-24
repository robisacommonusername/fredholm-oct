%Matlab doesn't have this function, octave does
function ret = is_function_handle(f)
    ret = isa(f,'function_handle');
end