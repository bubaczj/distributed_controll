function [A, B, X0, U0, range, activ] = linAB_TS_U_q(data, n, rng, q)

if nargin < 4
    q = 2;
end

range = rng * data.F10;
d0 = 2 * range / sum(q.^(1:n - 2));
band = data.F10 *(1 - rng) - d0;
activ = zeros(2,n);
for i = 1 : n
    ad = 0;
    if i > 1 
        ad = d0*sum(q.^(0:i-2));
    end
    F10 = band + ad + (d0/2)*q^(i-1); 
    [A(:,:,i), B(:,:,i), X0(:,i), U0(:,i)] = linAB_U(F10, data);
    activ(1,i) = F10;
    activ(2,i) = 1.25 * d0 * q^(i - 1);
end
end

