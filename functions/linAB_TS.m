function [A, B, X0, U0, range] = linAB_TS(data, n, rng)

range = rng * data.F10;

for i = 1 : n  
    F10(i) = range * (i - 1 - (n - 1) / 2) / (n - 1) + data.F10;
    [A(:,:,i), B(:,:,i), X0(:,i), U0(:,i)] = linAB_U(F10(i), data);
end
end

