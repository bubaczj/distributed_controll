function F = F_steps(t)
%%%%% Wymuszenie wieloma skokami
F = 0;
for i = 1 : 100
   if t > 1500 * (i - 1)
      F = (-1)^i * 30 * log(i + 1);
   else
       break
   end
end
F = [F;0];
end