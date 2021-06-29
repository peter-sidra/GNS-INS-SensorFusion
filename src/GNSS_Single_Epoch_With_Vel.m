% Load the data
ranges = csvread('Pseudo_ranges.csv');
range_rates = csvread('Pseudo_range_rates.csv');

time = ranges(2:end,1);
sat_nums = ranges(1,2:end);

[r_ea_e,v_ea_e] = pv_NED_to_ECEF(deg2rad(51.5074),deg2rad(0.1278),0,0);

[~,num_of_sats] = size(sat_nums);
sat_r_es_e = zeros(num_of_sats,3);
sat_v_es_e = zeros(num_of_sats,3);
predicted_ranges = zeros(num_of_sats,1);
predicted_range_rates = zeros(num_of_sats,1);
line_of_sight = zeros(num_of_sats,3);
position_solution = zeros(length(time),3);
velocity_solution = zeros(length(time),3);

% Initial guess
predicted_receiver_clock_offset = 0;
predicted_pos = [r_ea_e; predicted_receiver_clock_offset];
predicted_receiver_clock_drift = 0;
predicted_vel = [0;0;0;predicted_receiver_clock_drift];

epoch = 1;
% get sat positions and velocities at this epoch
for i = 1:num_of_sats
    [sat_r_es_e(i,:), sat_v_es_e(i,:)] = Satellite_position_and_velocity(time(epoch),sat_nums(i));
end

% Predict the ranges and range rates
for i = 1:num_of_sats
    % Predict the range
    sagnac_matrix = eye(3);
    
    temp = sagnac_matrix*transpose(sat_r_es_e(i,:)) - predicted_pos(1:3);
    predicted_ranges(i) = sqrt(transpose(temp) * temp);
    
    sagnac_matrix(1,2) = omega_ie*predicted_ranges(i)/c;
    sagnac_matrix(2,1) = -sagnac_matrix(1,2);
    
    temp = sagnac_matrix*transpose(sat_r_es_e(i,:)) - predicted_pos(1:3);
    predicted_ranges(i) = sqrt(transpose(temp) * temp);
    
    % Compute line of sight
    line_of_sight(i,:) = (sagnac_matrix*transpose(sat_r_es_e(i,:)) - predicted_pos(1:3))/predicted_ranges(i);
    
    % Predict the range rate
    predicted_range_rates(i) = line_of_sight(i,:)*(sagnac_matrix*(sat_v_es_e(i,:)'+Omega_ie*sat_r_es_e(i,:)')-...
        (predicted_vel(1:3) + Omega_ie*predicted_pos(1:3)));
end

delta_z_pos = zeros(num_of_sats,1);
H = zeros(num_of_sats, 4);
delta_z_vel = zeros(num_of_sats,1);

% Compute measurement matrix and measurement innovation vector
for i = 1:num_of_sats
    measured_sat_range = ranges(epoch+1,i+1);
    delta_z_pos(i) = measured_sat_range - predicted_ranges(i) - predicted_pos(4);
    H(i,1:3) = -line_of_sight(i,:);
    H(i,4) = 1;
    
    measured_sat_range_rate = range_rates(epoch+1,i+1);
    delta_z_vel(i) = measured_sat_range_rate - predicted_range_rates(i) - predicted_vel(4);
end

predicted_pos = predicted_pos + (H' * H) \ (H' * delta_z_pos);
predicted_vel = predicted_vel + (H' * H) \ (H' * delta_z_vel);

[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(predicted_pos(1:3,:),predicted_vel(1:3));
position_solution(epoch,:) = [rad2deg(L_b), rad2deg(lambda_b), h_b];
velocity_solution(epoch,:) = v_eb_n';

GNSS_Single_Epoch_Solution = [predicted_pos; predicted_vel];