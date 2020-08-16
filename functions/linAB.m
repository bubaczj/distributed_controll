function [A, B, X0, U0] = linAB(h20, data)
%X = [h1; h2]
%U = [F1; Fd]

h10 = h20 * (data.a2 / data.a1)^2;
X0 = [h10; h20];

U0 = [data.a1 * sqrt(h10) - data.Fd; data.Fd];

B = [1/data.A1, 1/data.A1; 0, 0];

A = [-data.a1 / (2 * data.A1 * sqrt(X0(1))), 0;
    data.a1 / (4 * data.C2 * X0(2) * sqrt(X0(1))), -data.a2 / (4 * data.C2 * X0(2)^(3/2))];

end

