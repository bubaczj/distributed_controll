function [U] = DMC_TS_U_q(t, x, des_h2, mpc, activ, start)

persistent  U_last dU1_hist  t_hist

if nargin > 5
        U_last = zeros(mpc.D - 1, 1);
        dU1_hist = zeros(mpc.D - 1, 1);
        t_hist = -(mpc.D - 2) - 1 : 1 : -1;
else
    t_run = floor(t);
    total_act = 0;
    for i = 1 : mpc.n
        if i == 1
            total_act = total_act + activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), -inf);
        elseif i == mpc.n 
            total_act = total_act + activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), inf);
        else
            total_act = total_act + activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1));
        end
    end
    if t_run > max(t_hist)
        mx = max(t_hist);
        t_hist = [t_hist, mx + 1 : t_run];
        [~, tidx] = size(t_hist);
        
        dU1_hist = [dU1_hist; zeros(t_run - mx, 1)];
        U_last = [U_last; U_last(end) * ones(t_run - mx, 1)];
        
        dU = 0;
    for i = 1 : mpc.n
        if i == 1

            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(tidx - (mpc.D - 1) : tidx - 1)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), -inf)/total_act;
            
        elseif i == mpc.n
            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(tidx - (mpc.D - 1) : tidx - 1)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), inf)/total_act;
        else
            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(tidx - (mpc.D - 1) : tidx - 1)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1))/total_act;
        end
        
        
    end
               
        U = U_last(tidx - 1) + dU;

        U_last = [U_last; U];
        dU1_hist = [dU1_hist; dU];

        U = [U; 0];
        
    else
        t_idx = find(t_hist <= t, 1, 'last') - 1;
        
        dU = 0;
    for i = 1 : mpc.n
        if i == 1

            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(t_idx - (mpc.D - 1) +1 : t_idx)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), -inf)/total_act;
            
        elseif i == mpc.n
            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(t_idx - (mpc.D - 1) +1 : t_idx)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1), inf)/total_act;
        else
            dU = dU + (mpc.ke_U(i) * (des_h2 - x(2)) - mpc.ku_U(:,i)' * flipud(dU1_hist(t_idx - (mpc.D - 1) +1 : t_idx)))...
                * activation(activ(1, i), activ(2,i), U_last(end) + mpc.U0(1))/total_act;
        end
        
        
    end
        
        
        
        U = U_last(t_idx) + dU;
        
        U_last(t_idx + 1) = U;
        dU1_hist(t_idx + 1) = dU;
        U = [U; 0];
    end
end

end