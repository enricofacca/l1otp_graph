function [ gfvar ] = tdens2gfvar( transformation ,power, tdens )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

switch transformation
    case 'identity'
        gfvar=tdens;
    case 'square'
        gfvar=sqrt(4*tdens);
    case 'squareplus'
        gfvar=real((2.1*tdens).^(1/2.1));
    case 'cube'
        gfvar=(3*tdens).^(1/3);
    case 'power'
        gfvar=real((power*tdens).^(1/power));
    case 'exp'
        gfvar=log(tdens);
    otherwise
        disp('Not supported')
end

if ( any(not (isreal(gfvar))))
    return
end
end

