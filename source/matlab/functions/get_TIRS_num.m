function [val] = get_TIRS_num(filename)
test = strfind(filename,'prefire_01');
val = 0; %zero is returned if the filename does not contain either string
if (test > 0) 
    val =1;
else
    test = strfind(filename,'prefire_02');
    if (test > 0)
        val=2;
    end
end

end
