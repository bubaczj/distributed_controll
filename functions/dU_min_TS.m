function [val] = dU_min_TS(dU_f, x, dU_hist, des_h2, mpc, i)
    x0 = x(2) * ones(mpc.N(i), 1) + mpc.Mp(:, :, i) * dU_hist;

    val = (des_h2 * ones(mpc.N(i), 1) - x0 - mpc.M(:,:,i) * dU_f)' * (des_h2 * ones(mpc.N(i), 1) - x0 - mpc.M(:,:,i) * dU_f) + mpc.lambda(i) * dU_f' * dU_f;
end