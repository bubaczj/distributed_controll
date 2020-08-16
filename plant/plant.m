function [dX] = plant(t, X, U, U0, data)
%X = [h1; h2]
%U = [F1; Fd]

U = U + U0;                                         % przyrostowe sterowanie, unifikacja z uk³¹dem zlinearyzownaym

for i = 1 : 2
    X(i) = max(eps, X(i));                          % nie mo¿e byæ bardziej pusty, ni¿ pusty, dodatkowo dzielenie przez zero
end


for i = 1 : 2
    U(i) = max(0, U(i));                            % wp³yw, ani zak³ócenie nie s¹ ujemne
end

dX(1, 1) = (U(1) + U(2) - data.a1 * (sqrt(X(1)))) / data.A1;

dX(2, 1) = (data.a1 * sqrt(X(1)) / X(2) - data.a2 / sqrt(X(2))) / (2 * data.C2);

end

