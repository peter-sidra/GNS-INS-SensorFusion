%% Load the data
data_ranges = csvread('Pseudo_ranges.csv');
data_range_rates = csvread('Pseudo_range_rates.csv');

ranges = data_ranges(2:end,2:end);
range_rates = data_range_rates(2:end,2:end);

time = data_ranges(2:end,1);
sat_nums = data_ranges(1,2:end);
[~,num_of_sats] = size(sat_nums);

% Sat positions and velocities
sat_r_es_e = zeros(3,num_of_sats);
sat_v_es_e = zeros(3,num_of_sats);
predicted_ranges = zeros(num_of_sats,1);
predicted_range_rates = zeros(num_of_sats,1);
line_of_sight = zeros(3,num_of_sats);

GNSS_Solution = zeros(length(time),8);

% Initialize state vector and error covariance
% Initialise state estimates
GNSS_Single_Epoch_Solution
State = [  GNSS_Single_Epoch_Solution(1:3);...
    GNSS_Single_Epoch_Solution(5:7);...
    GNSS_Single_Epoch_Solution(4); GNSS_Single_Epoch_Solution(8)]; 

% Initialise error covariance matrix
P =  zeros(8);
P(1,1) = 100;
P(2,2) = 100;
P(3,3) = 100;
P(4,4) = 0.01;
P(5,5) = 0.01;
P(6,6) = 0.01;
P(7,7) = (100000)^2;
P(8,8) = (200)^2;

% Compute transition matrix
tau_s = 0.5;
Phi = [eye(3), tau_s*eye(3), zeros(3,1), zeros(3,1);...
    zeros(3), eye(3), zeros(3,1), zeros(3,1);...
    zeros(1,3), zeros(1,3), 1, tau_s;...
    zeros(1,3), zeros(1,3), 0, 1];

% Compute the system noise covariant matrix
Sa = 5;
clock_phase_PSD = 0.01;
clock_frequency_PSD = 0.04;
Q = [1/3*Sa*tau_s^3*eye(3), 1/2*Sa*tau_s^2*eye(3), zeros(3,1), zeros(3,1);...
    1/2*Sa*tau_s^2*eye(3), Sa*tau_s*eye(3), zeros(3,1), zeros(3,1);...
    zeros(1,3), zeros(1,3), clock_phase_PSD*tau_s+1/3*clock_frequency_PSD*tau_s^3, 1/2*clock_frequency_PSD*tau_s^2;...
    zeros(1,3), zeros(1,3), 1/2*clock_frequency_PSD*tau_s^2, clock_frequency_PSD*tau_s];
    
for epoch = 1:length(time)
    % Propagate the state estimates
    State = Phi*State;
    
    % Propagate the error covariance matrix
    P = Phi*P*Phi' + Q;
    
    % Predict the ranges from the approximate user position to each
    % satellite
    for i = 1:num_of_sats   % Get sat pos and vel
        [sat_r_es_e(:,i), sat_v_es_e(:,i)] = Satellite_position_and_velocity(time(epoch),sat_nums(i));
    end
    
    for i = 1:length(sat_nums)
        sagnac_matrix = eye(3);
        predicted_ranges(i) = sqrt((sagnac_matrix*sat_r_es_e(:,i)-State(1:3))'*(sagnac_matrix*sat_r_es_e(:,i)-State(1:3)));
        
        % Update the sagnac matrix
        sagnac_matrix(1,2) = omega_ie*predicted_ranges(i)/c;
        sagnac_matrix(2,1) = -sagnac_matrix(1,2);
        
        predicted_ranges(i) = sqrt((sagnac_matrix*sat_r_es_e(:,i)-State(1:3))'*(sagnac_matrix*sat_r_es_e(:,i)-State(1:3)));
        
        % Compute the line of sight vectors
        line_of_sight(:,i) = (sagnac_matrix*sat_r_es_e(:,i)-State(1:3))/predicted_ranges(i);
        
        % Predict the range rates from the approximate user position to
        % each satellite
        predicted_range_rates(i) = line_of_sight(:,i)'*(sagnac_matrix*(sat_v_es_e(:,i) + Omega_ie*sat_r_es_e(:,i))-(State(4:6) + Omega_ie*State(1:3)));
    end
    
    % Compute the measurement noise covariance matrix
    range_measurement_std = 10;
    rate_measurement_std = 0.05;
    R = [eye(num_of_sats)*range_measurement_std^2, zeros(num_of_sats);...
        zeros(num_of_sats), eye(num_of_sats)*rate_measurement_std^2];
    
    % Compute the measurement matrix
    H = [-line_of_sight', zeros(size(line_of_sight')), ones(num_of_sats,1), zeros(num_of_sats,1);...
        zeros(size(line_of_sight')), -line_of_sight', zeros(num_of_sats,1), ones(num_of_sats,1)];
    
    % Formulate the measurement innovation vector
    delta_z = [ranges(epoch,:)'-predicted_ranges-State(7); range_rates(epoch,:)'-predicted_range_rates-State(8)];
    
    % residual based outlier detection
    [num_of_measurements,~] = size(delta_z);
    range_measurement_end_idx = num_of_measurements/2;
    while(1)
        [num_of_measurements,~] = size(delta_z);
        Im = eye(num_of_measurements);
        v = (H*((H'*H)\H')-Im)*delta_z;
        Cv = (Im-H*((H'*H)\H'));
        Cv_pos = Cv(1:range_measurement_end_idx,1:range_measurement_end_idx)*range_measurement_std^2;
        Cv_vel = Cv(range_measurement_end_idx+1:end,range_measurement_end_idx+1:end)*rate_measurement_std^2;
        w_pos = v(1:range_measurement_end_idx,:)./sqrt(diag(Cv_pos));
        w_vel = v(range_measurement_end_idx+1:end,:)./sqrt(diag(Cv_vel));
        pos_threshold = 5;
        vel_threshold = 1;
        % find the largest residuals
        [largest_residual_pos, largest_residual_pos_idx] = max(abs(w_pos));
        [largest_residual_vel, largest_residual_vel_idx] = max(abs(w_vel));
        
        largest_residual_vel_idx = largest_residual_vel_idx + range_measurement_end_idx;
        
        % check if it's an outlier based on the selected threshold value
        if(largest_residual_pos > pos_threshold && largest_residual_vel > vel_threshold)
            H(largest_residual_pos_idx,:) = [];
            delta_z(largest_residual_pos_idx,:) = [];
            R(largest_residual_pos_idx,:) = [];
            R(:,largest_residual_pos_idx) = [];
            H(largest_residual_vel_idx-1,:) = [];
            delta_z(largest_residual_vel_idx-1,:) = [];
            R(largest_residual_vel_idx-1,:) = [];
            R(:,largest_residual_vel_idx-1) = [];
            range_measurement_end_idx = range_measurement_end_idx - 1;
            disp(['Here 1 epoch is ' num2str(epoch)]);
        elseif(largest_residual_pos > pos_threshold)
            H(largest_residual_pos_idx,:) = [];
            delta_z(largest_residual_pos_idx,:) = [];
            R(largest_residual_pos_idx,:) = [];
            R(:,largest_residual_pos_idx) = [];
            disp(['Here 2 epoch is ' num2str(epoch)]);
            range_measurement_end_idx = range_measurement_end_idx - 1;
        elseif(largest_residual_vel > vel_threshold)
            H(largest_residual_vel_idx,:) = [];
            delta_z(largest_residual_vel_idx,:) = [];
            R(largest_residual_vel_idx,:) = [];
            R(:,largest_residual_vel_idx) = [];
            disp(['Here 3 epoch is ' num2str(epoch)]);
        else
            break;
        end
    end
    
    % Compute the kalman gain matrix
    K = P*H'/(H*P*H' + R);
    
    % Update the state estimates
    State = State + K*delta_z;
    
    % Update the error covaraince matrix
    P = (eye(size(K*H)) - K*H)*P;
    
    % Get long and lat
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(State(1:3),State(4:6));
    GNSS_Solution(epoch,1:6) = [L_b,lambda_b,h_b,v_eb_n'];
    GNSS_Solution(epoch,1:2) = rad2deg(GNSS_Solution(epoch,1:2));
    GNSS_Solution(epoch,7:8) = State(7:8);
end
