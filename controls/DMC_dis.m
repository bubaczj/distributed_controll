function [U] = DMC_dis(x, des_h2, mpc, start)

dis_amplitude = 0.01;
persistent  U_last dU1_hist dis_hist dis_last

if nargin > 3
        U_last = 0;
        dU1_hist = zeros(mpc.D - 1, 1);
        dis_hist = zeros(mpc.D_dis - 1, 1);
        dis_last = 0;
else
    t_run = floor(t);
    if t_run > max(t_hist)
        t_hist = [t_hist, max(t_hist) + 1 : t_run];
        [~, tidx] = size(t_hist);
        
        dU1_hist = [dU1_hist; zeros(t_run - max(t_hist), 1)];
        U_last = [U_last; U_last(end) * ones(t_run - max(t_hist), 1)];
        
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
    dis_hist = [dis_hist(2:end); ddis];
    dis_last = dis;
    dis = dis_amplitude * randn(1);
    ddis = dis - dis_last;    
    U = [U; dis];
end

end
