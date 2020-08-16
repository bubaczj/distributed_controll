function [A, B, X0, U0, range] = linAB_TS_X(data, n, rng)

range = rng * data.h20;

for i = 1 : n  
    h20(i) = range * (i - 1 - (n - 1) / 2) / (n - 1) + data.h20;  
    [A(:,:,i), B(:,:,i), X0(:,i), U0(:,i)] = linAB(h20(i), data);
end
end