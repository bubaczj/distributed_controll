function [des_h2] = desired_h2(t, h20)
    n = 2000;
    des_h2 = floor(t / n) * h20 / 5 + h20/2;
end

