function [U] = P(X, des_h2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
U = 2 * [des_h2 - X(2); 0];
end

