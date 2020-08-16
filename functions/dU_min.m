function [val] = dU_min(dU_f, x, dU_hist, des_h2, mpc)
    x0 = x(2) * ones(mpc.N, 1) + mpc.Mp * dU_hist;

    val = (des_h2 * ones(mpc.N, 1) - x0 - mpc.M * dU_f)' * (des_h2 * ones(mpc.N, 1) - x0 - mpc.M * dU_f) + mpc.lambda * dU_f' * dU_f;
end

