function [ t_easy ] = t_easy(B)
% t_easy Approximate Teasy/tau when VFS ~ VDD (see Calculating Teasy)
    t_easy = (1 + log(2)) .* (2.^(-B) - 1) + B .* (1 + (2.^(-B) - 0.5) .* log(2)) + B.^2 .* log(2) ./ 2;
end