%State-Space Signal Model for Kalman Filter.
% Inputs:
% N: Number of time steps.
% F: State transition matrix of size Nx x Nx.
% G: Input matrix of size Nx x Nu.
% C: Measurement matrix of size Nr x Nx.
% Q: Process noise covariance matrix
% R: Measurement noise covariance matrix
% x_Initial: Initial state vector of size Nx x 1.
% x_Estimate_Initial: Initial Estimate of state which is [0 ; 0]
% P_Estimate: the Error is estimated to be 1000
seed = (16+16+16+16);  %(Akshay Khanna) - a,a,a,a
rng(seed, 'twister');

i= 100;  % Number of iterations or time steps, value considered from the book 9.396
T = 0.1; %Time
max_Stepsize = 0:T:(i-1)*T; %Max Increment per second
F = [1 T ; 0 1]; %State Transition Matrix, value considered from the book 9.386
G = [T^2/2; T] ; %Control Transition Matrix, value considered from the book 9.387
sig2_Q = 40; %Variance of process noise
sig2_V = 40; %Variance of velocity, value considered from the book 9.395
sig2_R = 100; %Variance of measurement noise, value considered from the book 9.397
Q = sig2_Q; %Assigning the value to Q, u(k-1) ~ N(0,Q)
R = sig2_R; %Assigning the value to R 
C = [1 0]; %Measurement matrix of size, value considered from the book 9.390

x_trueVal = zeros(2, i);  %Initialize matrix 2xN to store true state vector values of the position and velocity.
x_noiseVal = zeros(1, i); %Initialize matrix 1xN to store position measurements with noise.
x_Initial = [1000; -50]; %Initialize [position(0); velocity(0)], values 9.393 from book

x_Estimate = zeros(2, i);
k_gainValue = zeros(2, i);
MSE = zeros(2, i);

x_Estimate_Initial = [0 ; 0]; % Initialize the estimate of x at time k-1
P_Estimate = 1000; % the Error is estimated to be 1000 i.e. 10*sig2_R

%-Calling the function and returning two values, True val, noisy val, kalman gain values, Mean squared error, new estimates valus of x-%


[x_trueVal, x_noiseVal, k_gainValue, MSE, x_Estimate] = kalman_signal_model_input(i, C, Q, R, F, G, x_Initial, x_trueVal, x_noiseVal, x_Estimate_Initial, P_Estimate);


%----Plotting the graph similar to the one in the textbook Figure 9.38----%

figure 
set(0,'DefaultLineLineWidth', 1.5, 'DefaultAxesFontSize',12, 'DefaultTextFontSize',12)
plot(max_Stepsize, x_trueVal(1, :), max_Stepsize, x_noiseVal, ':k', max_Stepsize, x_Estimate(1,:), '--')         
legend('True State Value', 'Noisy Value', 'Estimate')
xlabel('t(s)')
ylabel('Position(m)')
axis([0, 10, 400, 1000]) %Scaling similar to the fugure in the textbook 9.38                     
grid
%----Plotting the graph similar to the one in the textbook Figure 9.38----%
figure
plot(max_Stepsize, x_trueVal(2,:), max_Stepsize, x_Estimate(2,:), '--') 
legend('True','Estimate')
xlabel('t (s)')
ylabel('Velocity (m/s)')
grid
%----Plotting the graph similar to the one in the textbook Figure 9.39----%
figure  %Kalman Gain
plot(max_Stepsize, k_gainValue(1,:)) 
title('Kalman Gain - Position')
xlabel('t (s)')
ylabel('Gain for position')
axis([0,10,0,1.5])
grid
%----Plotting the graph similar to the one in the textbook Figure 9.39----%
figure
plot(max_Stepsize, k_gainValue(2,:)) 
xlabel('t (s)')
ylabel('Gain for velocity')
title('Kalman Gain - Velocity')
axis([0,10,0,1.5])
grid
%----Plotting the graph similar to the one in the textbook Figure 9.40----%
figure                                         
plot(max_Stepsize, x_trueVal(2, :))                              
legend('True')
xlabel('t(s)')
ylabel('Velocity(m/s)')
axis([0, 10, -100, 200]) %Scaling similar to the fugure in the textbook 9.38                           
grid
%----Plotting the graph similar to the one in the textbook Figure 9.40----%
figure                                         
plot(max_Stepsize, MSE(1, :))                              
legend('MSE for position')
xlabel('t(s)')
ylabel('Position(m)')
axis([0, 10, 0, 1000]) %Scaling similar to the fugure in the textbook 9.38                           
grid
%----Plotting the graph similar to the one in the textbook Figure 9.40----%
figure                                         
plot(max_Stepsize, MSE(2, :))                              
legend('MSE for velocity')
xlabel('t(s)')
ylabel('Velocity(m/s)')
axis([0, 10, 0, 1000]) %Scaling similar to the fugure in the textbook 9.38                           
grid

function [x_trueVal, x_noiseVal, k_gainValue, MSE, x_Estimate] = kalman_signal_model_input(n, C, Q, R, F, G, x_Initial, x_trueVal, x_noiseVal, x_Estimate_Initial, P_Estimate)

    for i = 1:n
        normal_Q = normrnd(0, sqrt(Q)); %Generate normal distributed process noise centered at 0
        normal_R= normrnd(0, sqrt(R)); %Generate normal distributed meaurement noise centred at 0
   
        x_Present = F * x_Initial + G * normal_Q; %True state vector values + noise, equation 9.212 from the book
        recievedSignal = C * x_Present + normal_R; %Current measured value of postion + noise, equation 9.213

        x_prediction = F * x_Estimate_Initial;
        P_prediction = F * P_Estimate * F' + G * Q * G'; %Addition of State transition matix and control transition matrix
        rs_prediction = C * x_prediction;

        rs_error = recievedSignal - rs_prediction;
        p_til = C * P_prediction * C' + R; %denominator for Kalman Gian
        kalman_gain = P_prediction * C' * inv(p_til); %Kalman Gian calculated

        x_Estimate_Initial = F * x_Estimate_Initial + kalman_gain * rs_error;
        P_Estimate = (eye(2) - kalman_gain * C) * P_prediction  ;
        
        x_trueVal(:, i) = x_Present; 
        x_Estimate(:, i) = x_Estimate_Initial;
        k_gainValue(:, i) = kalman_gain;
        MSE(1,i) = P_Estimate(1,1);
        MSE(2,i) = P_Estimate(2,2);
        x_noiseVal(i) = recievedSignal; %Assign to x_noisy 
        x_Initial = x_Present; %Updating state vector
    end
end
