function [U] = DMC(t, x, des_h2, mpc, start)

persistent  U_last dU1_hist  t_hist

if nargin > 4
        U_last = zeros(mpc.D(1) - 1, 1);
        dU1_hist = zeros(mpc.D(1) - 1, 1);
        t_hist = -(mpc.D(1) - 2) - 1 : 1 : -1;
else
    t_run = floor(t);
    if t_run > max(t_hist)
        mx = max(t_hist);
        t_hist = [t_hist, mx + 1 : t_run];
        [~, tidx] = size(t_hist);
        
        dU1_hist = [dU1_hist; zeros(t_run - mx, 1)];
        U_last = [U_last; U_last(end) * ones(t_run - mx, 1)];
        
        dU = mpc.ke * (des_h2 - x(2)) - mpc.ku * flipud(dU1_hist(tidx - (mpc.D(1) - 1) : tidx - 1));
        U = U_last(tidx - 1) + dU;

        U_last = [U_last; U];
        dU1_hist = [dU1_hist; dU];

        U = [U; 0];
    else
        t_idx = find(t_hist <= t, 1, 'last') - 1;
        dU = mpc.ke * (des_h2 - x(2)) - mpc.ku * flipud(dU1_hist(t_idx - (mpc.D(1) - 1) +1 : t_idx));
        U = U_last(t_idx) + dU;
        
        U_last(t_idx + 1) = U;
        dU1_hist(t_idx + 1) = dU;
        U = [U; 0];
    end
end

end

