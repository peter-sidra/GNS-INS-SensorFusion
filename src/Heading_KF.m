% Load the data
DR_Data = csvread('Dead_reckoning.csv');

time = DR_Data(:,1);
mag_readings = deg2rad(DR_Data(:,7));
gyro_readings = DR_Data(:,6);
% Integrate the gyro readings
gyro_readings = cumsum(gyro_readings/2);
% Initialize gyro with the magnetometer reading
gyro_readings = gyro_readings + mag_readings(1);
Integrated_heading_solution = zeros(length(time),1);


% Initialize the state estimates and error covariance matrix
X = zeros(2,1);

P = [deg2rad(4)^2,0;...
    0,deg2rad(1)^2];

% Compute the transition matrix
tau_s = 0.5;
Phi = [1, tau_s; 0, 1];

% Compute system noise covariance matrix
S_rg = (10^-4/sqrt(3600))^2;    % gyro random noise power spectral density
% S_bgd = (deg2rad(1)/sqrt(3600))^2;
S_bgd = 2*10^-12;
Q = [S_rg*tau_s + 1/3*S_bgd*tau_s^3, 1/2*S_bgd*tau_s^2;...
    1/2*S_bgd*tau_s^2, S_bgd*tau_s];

% Formulate the measurement matrix
H = [-1, 0];

% Compute the measurement noise covariance matrix
magnetic_heading_noise_variance = deg2rad(4)^2;
R = magnetic_heading_noise_variance;

for i = 1:length(time)
    % Propagate the state estimate
    X = Phi*X;
    
    % Propagate the error covariance estimate
    P = Phi*P*Phi' + Q;
    
    % Compute the kalman gain matrix
    K = P*H'/(H*P*H' + R);
    
    % Formulate the measurement innovation vector
    delta_z = (mag_readings(i) - gyro_readings(i)) - H*X;
    
    % Update the state estimates
    X = X + K*delta_z;
    
    Integrated_heading_solution(i) = gyro_readings(i) - X(1);
end

