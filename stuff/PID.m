function [U] = PID(x, des_h2, mpc, start)

persistent  uch last_x

if nargin > 3
        last_x = 0;
        uch = 0;
else
    uch = uch + (des_h2 - x(2));
    Kp = 5.3/3.9;   % Ku = 5.3
    Ki = 0.0012;
    Kd = 0.01;
    U = Kp * (des_h2 - x(2)) + Ki * uch + Kd * (x(2) - last_x);
    last_x = x(2);
    U = [U; 0];
end

end