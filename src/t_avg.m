function [ t_avg ] = t_avg(B)
% t_avg Approximate Tavg/tau when VFS ~ VDD (see Calculating Teasy, but instead uses full residual range)
    t_avg = B .* (1 + log(2)/2) + B.^2 .* log(2) ./ 2;
end