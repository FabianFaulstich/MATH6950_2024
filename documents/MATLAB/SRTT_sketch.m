function [c] = SRTT_sketch(b,d)

    n = length(b);
    signs = 2*randi(2,n,1)-3; % diagonal entries of D (Rademacher)
    idx = randsample(n,d); % indices i_1,...,i_d defining R

    % Multiply S against b
    c = signs .* b; % multiply by D
    c = dct(c); % multiply by F
    c = c(idx); % multiply by R
    c = sqrt(n/d) * c; % scale

end