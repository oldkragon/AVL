%% AVL Script
clear
close all
clc

N = 10;
b2 = 1;
theta = -2;

function y_panels = CalculatePanelLenghts(b2, N)
y_panels = b2.*(cos(pi/2 .* (0:N-1)/(N)));
y_panels = flip(y_panels)';
end

function twists = CalculateTwists(b2, y_pan, theta, fun)
twists = theta.*fun(y_pan)./fun(b2);
end

f = @(x) x;
g = @(x) x.^2;

y = CalculatePanelLenghts(b2, N)
twists_f = CalculateTwists(b2, y, theta, f)
twists_g = CalculateTwists(b2, y, theta, g)

plot(y, twists_g, y, twists_f)