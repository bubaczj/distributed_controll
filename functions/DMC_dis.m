function [U] = DMC_dis(x, des_h2, mpc, start)

dis_amplitude = 15;
persistent  U_last dU1_hist dis_hist dis_last

if nargin > 3
        U_last = 0;
        dU1_hist = zeros(mpc.D - 1, 1);
        dis_hist = zeros(mpc.D_dis - 1, 1);
        dis_last = 0;
else

    dU = mpc.ke * (des_h2 - x(2)) - mpc.ku * flipud(dU1_hist) - mpc.ku_dis * flipud(dis_hist);
    U = U_last + dU;
    
    dis = dis_amplitude * randn(1);
    ddis = dis - dis_last;
    
    U_last = U;
    dU1_hist = [dU1_hist(2:end); dU];
    dis_hist = [dis_hist(2:end); ddis];
    dis_last = dis;
    
    U = [U; dis];
end

end
