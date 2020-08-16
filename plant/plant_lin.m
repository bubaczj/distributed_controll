function [dX] = plant_lin(t, X, U, X0, A, B)
%X = [h1; h2]
%U = [F1; Fd]

dX = A*(X - X0) + B*U;                  % model liniowy

if X(1) <= 0                            % zbiorniki nie mog¹ byæ bardziej puste, ni¿ puste
    dX(1, 1) = max(dX(1,1), 0);
end

if X(2) <= 0
    dX(2, 1) = max(dX(2,1), 0);
end

end

