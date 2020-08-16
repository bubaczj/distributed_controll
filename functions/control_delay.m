function [U_del, t_hist, U_hist] = control_delay(t, U, start)

persistent U_val U_time

if nargin > 2
    
    if start == 0
        U_time = [0];
        U_val = [0; 0];
    elseif start == 1
        U_hist = U_val;
        t_hist = U_time;
        U_del = [];
    end
    
else
    
    U_time = [U_time, t];
    U_val = [U_val, U];
    
    if t < 80 
        U_del = [0; U(2)];
    else
        t_index = find(U_time <= t-80, 1, 'last');
        U_del = [U_val(1, t_index); U(2)];
    end
    
end

end

