%%
time = GNSS_Solution(:,1);
[ndata, ~] = size(time);

% GNSS data
GNSS_lat = deg2rad(GNSS_Solution(:,1));
GNSS_lon = deg2rad(GNSS_Solution(:,2));
GNSS_height = GNSS_Solution(:,3);
GNSS_vel = GNSS_Solution(:,4:6);

[R_N_0, R_E_0] = Radii_of_curvature(GNSS_lat(1));

% Initialize state estimation error covariance matrix
P = zeros(4);
P(1,1) = 0.1^2;
P(2,2) = 0.1^2;
P(3,3) = 10^2/(R_N_0+GNSS_height(1))^2;
P(4,4) = 10^2/(R_E_0+GNSS_height(1))^2/cos(GNSS_lat(1))^2;

% Initialize the state
X = zeros(4,1);

S_DR = 0.2;
tau_s = 0.5;

% DR Data
DR_Data = csvread('Dead_reckoning.csv');
time = DR_Data(:,1);
[ndata, ~] = size(time);
avg_speed = mean(DR_Data(:,2:5),2);
heading = Integrated_heading_solution;

avg_velocity = zeros(ndata,2);
inst_velocity = zeros(ndata,2);
longtitude = zeros(ndata,1);
latitude = zeros(ndata,1);
% Initialize the DR solution using the GNSS solution
longtitude(1) = GNSS_lon(1);
latitude(1) = GNSS_lat(1);
inst_velocity(1,:) = GNSS_vel(1,1:2);

% Filtered solution
longtitude_corrected = zeros(ndata,1);
latitude_corrected = zeros(ndata,1);
inst_velocity_corrected = zeros(ndata,2);

for i = 1:ndata
    curr_heading = heading(i);
    curr_time = time(i);
    curr_longtitude = longtitude(i);
    curr_latitude = latitude(i);
    curr_height = GNSS_height(i);
    if(i==1)
        prev_time = time(i);
        prev_heading = heading(i);
        prev_longtitude = longtitude(i);
        prev_latitude = latitude(i);
%         prev_inst_vel = [avg_speed(i)*cos(curr_heading), avg_speed(i)*sin(curr_heading)];
        prev_inst_vel = inst_velocity(i,:);
        prev_height = GNSS_height(i);
    else
        prev_time = time(i-1);
        prev_heading = heading(i-1);
        prev_longtitude = longtitude(i-1);
        prev_latitude = latitude(i-1);
        prev_inst_vel = inst_velocity(i-1,:);
        prev_height = GNSS_height(i-1);
    end
    [R_N, R_E] = Radii_of_curvature(prev_latitude);
    
    % Compute the DR solution
    % ---> Calc avg and inst velcoty
    avg_velocity(i,:) = 1/2*[cos(curr_heading) + cos(prev_heading),...
        sin(curr_heading) + sin(prev_heading)]*avg_speed(i);
    
    d = 0.5;    
    inst_velocity(i,:) = (2-d)*avg_velocity(i,:) - (1-d)*prev_inst_vel;
    
    % ---> Calc long and lat
    latitude(i) = prev_latitude + avg_velocity(i,1)*(curr_time-prev_time)/(R_N + curr_height);
    longtitude(i) = prev_longtitude + avg_velocity(i,2)*(curr_time-prev_time)/(R_E + curr_height)/cos(longtitude(i));
    curr_heading = heading(i);
    curr_time = time(i);
    curr_longtitude = longtitude(i);
    curr_latitude = latitude(i);
    curr_height = GNSS_height(i);
    if(i==1)
        prev_time = time(i);
        prev_heading = heading(i);
        prev_longtitude = longtitude(i);
        prev_latitude = latitude(i);
        prev_inst_vel = [avg_speed(i)*cos(curr_heading), avg_speed(i)*sin(curr_heading)];
        prev_height = GNSS_height(i);
    end
    [R_N, R_E] = Radii_of_curvature(prev_latitude);

    % Compute the transition matrix
    Phi = eye(4);
    Phi(3,1) = tau_s/(R_N+prev_height);
    Phi(4,2) = tau_s/((R_E+prev_height)*cos(prev_latitude));
    
    % Compute the system noise covariance matrix
    Q = [S_DR*tau_s, 0, 1/2*(S_DR*tau_s^2)/(R_N+prev_height), 0;...
        0, S_DR*tau_s, 0, 1/2*(S_DR*tau_s^2)/((R_E+prev_height)*cos(prev_latitude));...
        1/2*(S_DR*tau_s^2)/(R_N+prev_height), 0, 1/3*(S_DR*tau_s^3)/(R_N+prev_height)^2,0;...
        0, 1/2*(S_DR*tau_s^2)/(R_E+prev_height)/cos(prev_latitude), 0, 1/3*(S_DR*tau_s^3)/(R_E+prev_height)^2/cos(prev_latitude)^2];
    
    % Propagate the state estimates
    X = Phi*X;
    
    % Propagate the error covariance matrix
    P = Phi*P*Phi' + Q;
    
    %Compute the measurement matrix
    H = zeros(4);
    H(1,3) = -1;
    H(2,4) = -1;
    H(3,1) = -1;
    H(4,2) = -1;
    
    % Compute the measurement noise covariance matrix
    sigma_Gr = 10;       % GNSS position measurement std deviation
    sigma_Gv = 0.05;    % GNSS velocity measurement std deviation
    R = zeros(4);
    [R_N, R_E] = Radii_of_curvature(GNSS_lat(i));
    R(1,1) = sigma_Gr^2/(R_N + curr_height)^2;
    R(2,2) = sigma_Gr^2/((R_E + curr_height)^2*cos(GNSS_lat(i))^2);
    R(3,3) = sigma_Gv^2;
    R(4,4) = sigma_Gv^2;
    
    % Compute the kalman gain matrix
    K = P*H'/(H*P*H' + R);
    
    % Formulate the measurement innovation vector
    delta_z = [GNSS_lat(i) - latitude(i);...
        GNSS_lon(i) - longtitude(i);...
        GNSS_vel(i,1:2)' - inst_velocity(i,:)'];
    delta_z = delta_z - H*X;
    
    % Update the state estimate
    X = X + K*delta_z;
    
    % Update the error covariance matrix
    P = (eye(size(K*H)) - K*H)*P;
    
%     Apply corrections via feedback to the DR processor every 50 epochs
%     and reset the state estimates
    if(mod(i,50) == 0)
        latitude(i) = latitude(i) - X(3);
        longtitude(i) = longtitude(i) - X(4);
        inst_velocity(i,:) = inst_velocity(i,:) - X(1:2)';
        
        % Reset the state estimates
        X = zeros(4,1);
    end

    % Filtered result
    latitude_corrected(i) = latitude(i) - X(3);
    longtitude_corrected(i) = longtitude(i) - X(4);
    inst_velocity_corrected(i,:) = inst_velocity(i,:) - X(1:2)';
    
end