%% purpose: show how kalman filtering works
% based on: https://www.kalmanfilter.net/kalman1d.html
close all
clear all
rng('default')

% declare
N = 20;
measurements = zeros(1,N);
measurements = measurements + 5 + 2*(rand(1,N)-0.5);
uncert = 1.5;
process_noise_var = 0.1;
kalman_gains = zeros(1,N);
x_hats = zeros(1,N+1);
error_est = zeros(1,N+1);

% seed
x_hats(1) = 10;
error_est(1) = 100;

%loop
for i=1:N
    kalman_gains(i) = error_est(i)/(error_est(i)+uncert);
    x_hats(i+1) = x_hats(i)+kalman_gains(i)*(measurements(i)-x_hats(i));
    error_est(i+1) = (1-kalman_gains(i))*error_est(i) + process_noise_var;
end

figure(1)
subplot(3,1,1)
plot(x_hats(2:end))
hold on
grid on
plot(measurements)
axis([1 N 2 8])
title("Kalman filter effect on link distance. True value = 5m")
ylabel("Meters (m)")
xlabel("Measurement num")
legend("Filter output", "Measurements")

subplot(3,1,2)
plot(kalman_gains)
axis([1 N 0 1])
title("Kalman filter gain")
ylabel("Kalman gain")
xlabel("Measurement num")

subplot(3,1,3)
plot(error_est(2:end))
axis([1 N 0 1.6])
title("Estimate uncertainty")
ylabel("Estimate uncertainty")
xlabel("Measurement num")









