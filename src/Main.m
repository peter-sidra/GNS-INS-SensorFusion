clear;

% Load the constants
Define_Constants;

% Use single epoch least squares GNSS to initialize the GNSS Kalman Filter
GNSS_Single_Epoch_With_Vel;

% Obtain GNSS solution using multi-epoch kalman filter
GNSS_Kalman_Filter_With_Outlier_Detection;

% Obtain heading solution from Gyro-Magnetometer Integration using Kalman
% Filter
Heading_KF;

% Integrate DR/GNSS using a kalman filter
DR_GNSS_Integration;

Integrated_Solution = [time, rad2deg(latitude_corrected), rad2deg(longtitude_corrected), inst_velocity_corrected, rad2deg(Integrated_heading_solution)];

% Save the solution
csvwrite('Integrated_Solution.csv',Integrated_Solution);

%%

start_point = 1;
end_point = length(Integrated_Solution);
% end_point = 300;
% Draw position
hold on;
figure(1);
geoshow(rad2deg(latitude_corrected(start_point:end_point)),rad2deg(longtitude_corrected(start_point:end_point)),'Color','b');
pbaspect([1 1 1]);
daspect([1 1 1]);
plot(rad2deg(longtitude_corrected(start_point)),rad2deg(latitude_corrected(start_point)),'g*');
plot(rad2deg(longtitude_corrected(end_point)),rad2deg(latitude_corrected(end_point)),'r*');

% Draw velocity vectors
hold on
figure(1);
quiver(rad2deg(longtitude_corrected(start_point:4:end_point)),rad2deg(latitude_corrected(start_point:4:end_point)),...
    inst_velocity_corrected(start_point:4:end_point,2),inst_velocity_corrected(start_point:4:end_point,1),0.5, 'Color', 'k');
pbaspect([1 1 1])
daspect([1 1 1])
% 
% % Draw heading arrows
l = 1;
v = l * cos(Integrated_heading_solution(start_point:4:end_point));
u = l * sin(Integrated_heading_solution(start_point:4:end_point));
quiver(rad2deg(longtitude_corrected(start_point:4:end_point)),rad2deg(latitude_corrected(start_point:4:end_point)),u,v,0.5, 'Color', 'm');

pbaspect([1 1 1])
daspect([1 1 1])

legend('Position', 'Start', 'Finish', 'Velocity', 'Heading');

% Draw GNSS_Only solution
% hold on;
% figure(1);
% geoshow(rad2deg(GNSS_lat(start_point:end_point)),rad2deg(GNSS_lon(start_point:end_point)),'Color','red');
% pbaspect([1 1 1]);
% daspect([1 1 1]);
% plot(rad2deg(GNSS_lon(start_point)),rad2deg(GNSS_lat(start_point)),'*');
% plot(rad2deg(GNSS_lon(end_point)),rad2deg(GNSS_lat(end_point)),'*');
% 
% Draw DR only position solution
% hold on;
% figure(1);
% geoshow(rad2deg(latitude(start_point:end_point)),rad2deg(longtitude(start_point:end_point)),'Color','red');
% pbaspect([1 1 1]);
% daspect([1 1 1]);
% plot(rad2deg(longtitude(start_point)),rad2deg(latitude(start_point)),'g*');
% plot(rad2deg(longtitude(end_point)),rad2deg(latitude(end_point)),'r*');
% legend;

% legend('Position', 'Start', 'Finish', 'GNSS Solution', 'GNSS Start', 'GNSS Finish');
% 
pbaspect([1 1 1])
daspect([1 1 1])

