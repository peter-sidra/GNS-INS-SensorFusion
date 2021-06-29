%% Load constants
Define_Constants;

ranges = csvread('Pseudo_ranges.csv');
range_rates = csvread('Pseudo_range_rates.csv');

time = ranges(2:end,1);
satNums = ranges(1,2:end);

[r_ea_e,v_ea_e] = pv_NED_to_ECEF(deg2rad(51.5074),deg2rad(0.1278),0,0);

[~,num_of_sats] = size(satNums);
sat_r_es_e = zeros(num_of_sats,3);
sat_v_es_e = zeros(num_of_sats,3);

% get sat positions at time 0
for i = 1:num_of_sats
    [sat_r_es_e(i,:), sat_v_es_e(i,:)] = Satellite_position_and_velocity(time(1),satNums(i));
end

r_hat_aj = zeros(num_of_sats,1);
u_aj_e = zeros(num_of_sats,3);

while(1)
    for i = 1:num_of_sats
        sagnac_matrix = eye(3);
        
        temp = sagnac_matrix*transpose(sat_r_es_e(i,:)) - r_ea_e;
        r_hat_aj(i) = sqrt(transpose(temp) * temp);
        
        sagnac_matrix(1,2) = omega_ie*r_hat_aj(i)/c;
        sagnac_matrix(2,1) = -sagnac_matrix(1,2);
        
        temp = sagnac_matrix*transpose(sat_r_es_e(i,:)) - r_ea_e;
        r_hat_aj(i) = sqrt(transpose(temp) * temp);
        
        % compute line of sight
        u_aj_e(i,:) = (sagnac_matrix*transpose(sat_r_es_e(i,:)) - r_ea_e)/r_hat_aj(i);
    end
    
    predicted_receiver_clock_offset = 0;
    x_hat_minus = [r_ea_e; predicted_receiver_clock_offset];
    delta_z = zeros(num_of_sats,1);
    H_G_e = zeros(num_of_sats, 4);
    
    for i = 1:num_of_sats
        measured_sat_range = ranges(2,i+1);
        delta_z(i) = measured_sat_range - r_hat_aj(i) - predicted_receiver_clock_offset;
        H_G_e(i,1:3) = -u_aj_e(i,:);
        H_G_e(i,4) = 1;
    end
    
    % x_hat_plus = x_hat_minus + inv(transpose(H_G_e) * H_G_e) * transpose(H_G_e) * delta_z
    x_hat_plus = x_hat_minus + (transpose(H_G_e) * H_G_e) \ (transpose(H_G_e) * delta_z);
    if(norm(x_hat_plus(1:3) - x_hat_minus(1:3))<0.1)
        break;
    end
    r_ea_e = x_hat_plus(1:3);
end

GNSS_Single_Epoch_Solution = [L_b,lambda_b,h_b,v_eb_n',x_hat_plus(4)];
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_hat_plus(1:3,:),zeros(3,1));






