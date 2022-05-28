%% purpose: show how bayesian fusion works
% from: https://www.mathworks.com/matlabcentral/answers/370386-how-to-plot-a-gaussian-1d-in-matlab
close all

gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

x = linspace(-5,25,100);

mu1 = 9.6;
sig1 = 1.2;
amp1 = 0.5;
vo1 = 0; 
y = gaus(x,mu1,sig1,amp1,vo1);
plot(x, y, 'b-', 'LineWidth',3)
hold on
grid on

mu2 = 4.2;
sig2 = 1.5;
amp2 = 0.4;
vo2 = 0; 
y = gaus(x,mu2,sig2,amp2,vo2);
plot(x, y, 'g-', 'LineWidth',3)

mu3 = mu1*sig1^2/(sig1^2+sig2^2) + mu2*sig2^2/(sig1^2+sig2^2);
sig3 = (sig1^-2 + sig2^-2)^-0.5*(5.6^2/(5.6^2-(mu1-mu2)^2))^0.5;
amp3 = 0.15;
vo3 = 0; 
y = gaus(x,mu3,sig3,amp3,vo3);
plot(x, y, 'r-', 'LineWidth',3)

% title(sprintf('Gaussian with \\mu=%.1f \\sigma=%.1f amp=%.1f vo=%.1f', ...
%     mu, sig, amp, vo))
title("Bayesian fusion of RSSI and ToF data. True value = 4m")
axis([0 14 0 0.6])
legend("RSSI (spurious)","ToF (accepted)","Combined (rejected)");


