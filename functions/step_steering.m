function des_h = step_steering(t)
des_h = 0;
for i = 1 : 100
   if t > 1000 * (i - 1)
      des_h = (-1)^(i + 1) * 15;
   else
       break
   end
end
%des_h = 33/2;
des_h = des_h + 33;
end