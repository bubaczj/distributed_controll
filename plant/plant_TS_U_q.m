function [dX] = plant_TS_U_q(t, X, U, X0, U0, A, B, n, U0_glob, activ)
%fuzed by U(1)
% n > 1

dX = [0; 0];
act_sum = 0;
for i = 1 : n
    if i == 1
        act_sum = act_sum + activation(activ(1,i), activ(2,i), U(1) + U0_glob(1), -inf);
    elseif i == n
        act_sum = act_sum + activation(activ(1,i), activ(2,i), U(1) + U0_glob(1), inf);
    else
        act_sum = act_sum + activation(activ(1,i), activ(2,i), U(1) + U0_glob(1));
    end
end
for i = 1 : n
    if i == 1
       
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(activ(1,i), activ(2,i), U(1) + U0_glob(1), -inf)/act_sum;

    elseif i == n
        
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(activ(1,i), activ(2,i), U(1) + U0_glob(1), inf)/act_sum;

    else
        dX = dX + plant_lin(t, X, U + U0_glob - U0(:,i) , X0(:,i), A(:,:,i), B(:,:,i))...
            * activation(activ(1,i), activ(2,i), U(1) + U0_glob(1))/act_sum;

    end
end

end

