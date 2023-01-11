function [ tdens] = gfvar2tdens( transformation ,power,gfvar)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

switch transformation
    case 'identity'
        tdens=gfvar;
    case 'square'
        tdens=(gfvar.^2)/4;
    case 'squareplus'
        tdens=real((gfvar.^2.1)/2.1);
    case 'cube'
        tdens=(gfvar.^3)/3;  
    case 'power'
        tdens=real((gfvar.^power))/power;
    case 'exp'
        tdens=exp(gfvar);
    otherwise
        disp('Not supported')
end

if ( any(not (isreal(tdens))))
    return
end

end

