function [dX] = plant_TS_X(t, X, U, X0, U0, A, B, n, U0_glob, range)
%fuzed by U(1)
% n > 1

dX = [0; 0];
for i = 1 : n
    if i == 1
       
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(X0(2, i), 2 * range / (n - 1), X(2), -inf);

    elseif i == n
        
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(X0(2, i), 2 * range / (n - 1), X(2), inf);

    else
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(X0(2, i), 2 * range / (n - 1), X(2));

    end
end

end

