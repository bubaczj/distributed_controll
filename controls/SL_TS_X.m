function [U] = SL_TS_X(t, x, des_h2, mpc, start)

persistent  U_last dU1_hist t_hist dU0 options
%b ---- -y0
if nargin > 4
        U_last = zeros(mpc.D_X(1) - 1, 1);
        dU1_hist = zeros(mpc.D_X(1) - 1, 1);
        t_hist = -(mpc.D_X(1) - 2) - 1 : 1 : -1;
        dU0 = zeros(mpc.Nu_X(1), 1);
        options = optimoptions('fmincon','Display','off','StepTolerance',1e-2);
        
else
    t_run = floor(t);
    
    if t_run > max(t_hist)
        t_hist = [t_hist, max(t_hist) + 1 : t_run];
        [~, tidx] = size(t_hist);
        
        dU1_hist = [dU1_hist; zeros(t_run - max(t_hist), 1)];
        U_last = [U_last; U_last(end) * ones(t_run - max(t_hist), 1)];
        
        sl.Mp = zeros(mpc.N_X(1), mpc.D_X(1) - 1);
        sl.M = zeros(mpc.N_X(1), mpc.Nu_X(1));
        sl.w = zeros(mpc.n_X(1));
            for i = 1 : mpc.n_X(1)
                if i == 1
                    w(i) = activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2), -inf);
                elseif i == mpc.n_X           
                    w(i) =  activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2), inf);
                else
                    w(i) = activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2));
                end
            end 
            
        for i = 1 : mpc.n_X(1)
            if w(i) > 0
                sl.M = sl.M + w(i) .* mpc.M_X(:,:,i);
                sl.Mp = sl.Mp + w(i) .* mpc.Mp_X(:,:,i);
            end
        end
        
        sl.N = mpc.N_X(1);
        sl.lambda = mpc.lambda_X(1);
        
        y0 = x(2) * ones(mpc.N_X(1), 1) + sl.Mp * flipud(dU1_hist(tidx - (mpc.D_X(1) - 1) : tidx - 1));

        A = [eye(mpc.Nu_X(1)); tril(ones(mpc.Nu_X(1))); sl.M];
        A = [A; -A];

        b = [15 * ones(mpc.Nu_X(1), 1); (150  - (U_last(end) + 80)) * ones(mpc.Nu_X(1), 1); (50 - y0);...
            15 * ones(mpc.Nu_X(1), 1); (-0 + (U_last(end) + 80)) * ones(mpc.Nu_X(1), 1);  y0 .* 1];

        dU0 = fmincon(@(dU)dU_min(dU, x, flipud(dU1_hist(tidx - (mpc.D_X(1) - 1) : tidx - 1)), des_h2, sl),...
            dU0, A, b,[],[],[],[],[], options);
        dU = dU0(1);
     
    
        U = U_last(tidx - 1) + dU;
        U_last = [U_last; U];
        dU1_hist = [dU1_hist; dU];
        U = [U; 0];
        
    else
        t_idx = find(t_hist <= t, 1, 'last') - 1;
        
        sl.Mp = zeros(mpc.N_X(1), mpc.D_X(1) - 1);
        sl.M = zeros(mpc.N_X(1), mpc.Nu_X(1));
        w = zeros(mpc.n_X(1));
        
        for i = 1 : mpc.n_X(1)
            if i == 1
                w(i) = activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2), -inf);
            elseif i == mpc.n_X
                w(i) =  activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2), inf);
            else
                w(i) =  activation(mpc.X0_TS_X(2, i), 2 * mpc.range_X / (mpc.n_X - 1), x(2));
            end
        end
        
        for i = 1 : mpc.n_X(1)
            if w(i) > 0
                sl.M = sl.M + w(i) .* mpc.M_X(:,:,i);
                sl.Mp = sl.Mp + w(i) .* mpc.Mp_X(:,:,i);
            end
        end
        sl.N = mpc.N_X(1);
        sl.lambda = mpc.lambda_X(1);
        
        y0 = x(2) * ones(mpc.N_X(1), 1) + sl.Mp * flipud(dU1_hist(t_idx - (mpc.D_X(1) - 1) + 1 : t_idx));

        A = [eye(mpc.Nu_X(1)); tril(ones(mpc.Nu_X(1))); sl.M];
        A = [A; -A];

        b = [15 * ones(mpc.Nu_X(1), 1); (150  - (U_last(end) + 80)) * ones(mpc.Nu_X(1), 1); (50 - y0);...
            15 * ones(mpc.Nu_X(1), 1); (-0 + (U_last(end) + 80)) * ones(mpc.Nu_X(1), 1);  y0.*1];

        dU0 = fmincon(@(dU)dU_min(dU, x, flipud(dU1_hist(t_idx - (mpc.D_X(1) - 1) + 1 : t_idx)), des_h2, sl),...
            dU0, A, b,[],[],[],[],[], options);
        dU = dU0(1);
        
        
        U = U_last(t_idx) + dU;
        
        U_last(t_idx + 1) = U;
        dU1_hist(t_idx + 1) = dU;
        U = [U; 0];
    end
    
end

end