function [trans_prime,trans_second ] = gfvar2f(transformation ,power, gfvar)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

ntdens=size(gfvar,1);
switch transformation
    case 'identity'
        trans_prime  = ones(ntdens,1);
        trans_second = zeros(ntdens,1);
    case 'square'
        trans_prime  = gfvar/2;
        trans_second = ones(ntdens,1)/2;
    case 'squareplus'
        trans_prime  = real(gfvar.^(1.1));
        trans_second = real(1.1*gfvar.^(0.1));
    case 'power'
        trans_prime  = real(gfvar.^(power-1));
        trans_second = real((power-1)*gfvar.^(power-2));
    case 'cube'
        trans_prime  = gfvar.^2;
        trans_second = 2*gfvar;
    case 'exp'
        trans_prime  = exp(gfvar);
        trans_second = exp(gfvar);
    otherwise
        disp('Not supported')
end

if ( any(not (isreal(trans_prime))))
    return
end

if ( any(not (isreal(trans_second))))
    return
end

end

