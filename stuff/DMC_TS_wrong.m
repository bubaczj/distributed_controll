function [U] = DMC_TS(x, des_h2, mpc, start)

persistent  U_last dU1_hist 

if nargin > 3
        U_last = 0;
        dU1_hist = zeros(mpc.D(1) - 1, 1);
else
    dU = 0;
    for i = 1 : mpc.n
        if i == 1

            dU = dU + (mpc.ke(i) * (des_h2 - x(2)) - mpc.ku(:,i)' * flipud(dU1_hist))...
                * activation(mpc.U0_TS(1, i), 2 * mpc.range / (mpc.n - 1), U_last + mpc.U0(1), -inf);
            
        elseif i == mpc.n
            dU = dU + (mpc.ke(i) * (des_h2 - x(2)) - mpc.ku(:,i)' * flipud(dU1_hist))...
                * activation(mpc.U0_TS(1, i), 2 * mpc.range / (mpc.n - 1), U_last + mpc.U0(1), inf);
        else
            dU = dU + (mpc.ke(i) * (des_h2 - x(2)) - mpc.ku(:,i)' * flipud(dU1_hist))...
                * activation(mpc.U0_TS(1, i), 2 * mpc.range / (mpc.n - 1), U_last + mpc.U0(1));
        end
        
        
    end
    U = U_last + dU;
    U_last = U;
    dU1_hist = [dU1_hist(2:end); dU];
    
    
    U = [U; 0];
end
