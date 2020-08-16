function [U] = DMC(x, des_h2, mpc, start)

persistent  U_last dU_hist

if nargin > 3
        U_last = 0;
        dU_hist = zeros(mpc.D - 1, 1);
else
    dU = (mpc.ke * (des_h2 - x(2)) - mpc.ku * flipud(dU_hist));
    U = U_last + dU;
    
    U_last = U;
    dU_hist = [dU_hist(2:end); dU];
    
    U = [U; 0];
end

end

