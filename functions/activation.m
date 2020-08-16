function [ x_out ] = activation( mid, range, x_in, edge )

% mid - œrodek przedzia³u
% range - szerokoœæ przedzia³u
% edge - warunek krañcowoœci edge = -inf lewy kraniec, edge = +inf prawy
% koniec


A = mid - range/2;
B = mid + range/2;

dr = 2/(B-A);
df = -dr;
if nargin == 3
    if x_in <= A || x_in >= B
        x_out = 0;
    else
        if x_in < (A + B)/2
            x_out = dr * (x_in - A);
        else
            x_out = 1 + df * (x_in - (A + B)/2);
        end
    end
end
if nargin > 3
    
   dr = 2/(B-A);
   if edge == -inf
       A = A + range/2;
       if x_in <= A
          x_out = 1; 
       else
           if x_in >= B
               x_out = 0;
           else
               x_out = 1 - dr * (x_in - A);
           end
       end
   end

   if edge == inf
       B = B -range / 2;
       if x_in <= A
          x_out = 0; 
       else
           if x_in >= B
               x_out = 1;
           else
               x_out = dr * (x_in - A);
           end
       end
   end
end

end

