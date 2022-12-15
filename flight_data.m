function [acc_imu,vel_gps,ang_rates,psi_mag,ang_imu,ang_hyb,gps_pose,imu_vel, dT,T] = flight_data()
%arinc = load('arinc.mat');
ingps = load('ingps.mat');

time = ingps.ingps.INGPS_time*3600;
dt = diff(time);

dT = dt';
t = time - time(1);
T = t';
%% IMU Data

g_x_imu = ingps.ingps.GMAX; 
g_y_imu = ingps.ingps.GMAY;
g_z_imu = ingps.ingps.GMAZ;

acc_imu = [g_x_imu'; g_y_imu'; g_z_imu'];    % acceleration FR?

vel_N_imu = ingps.ingps.VN_INS;           
vel_E_imu = ingps.ingps.VE_INS;
vel_V_imu = ingps.ingps.VV_INS;

imu_vel = [vel_N_imu';vel_E_imu';-vel_V_imu'];

phi_imu = deg2rad(ingps.ingps.ROLL_INS);
theta_imu = deg2rad(ingps.ingps.PITCH_INS);
psi_imu = deg2rad(ingps.ingps.THEAD_INS_c);

ang_imu = [phi_imu'; theta_imu'; psi_imu'];

p_imu = deg2rad(ingps.ingps.OMEGAX);
q_imu = deg2rad(ingps.ingps.OMEGAY);
r_imu = deg2rad(ingps.ingps.OMEGAZ);

ang_rates = [p_imu'; q_imu';r_imu'];

lat_imu = ingps.ingps.LAT_INS;
lon_imu = ingps.ingps.LON_INS;
%[pos_x_imu,pos_y_imu] = grn2eqa(lat_imu,lon_imu);

%% MAG Data

psi_mag = deg2rad(ingps.ingps.THEAD_INS_c);

%% GPS Data

vel_N_GPS = ingps.ingps.VN_SAT;
vel_E_GPS = ingps.ingps.VE_SAT;
vel_V_GPS = ingps.ingps.VV_SAT;
vel_D_GPS = -vel_V_GPS;

vel_gps = [vel_N_GPS'; vel_E_GPS'; vel_D_GPS'];

lat_gps = ingps.ingps.LAT_SAT;
lon_gps = ingps.ingps.LON_SAT;

%% Fused Data

lat_hyb = ingps.ingps.LAT_HYB;
lon_hyb = ingps.ingps.LON_HYB;

phi_hyb = deg2rad(ingps.ingps.ROLL_HYB);
theta_hyb = deg2rad(ingps.ingps.PITCH_HYB);
psi_hyb = deg2rad(ingps.ingps.THDG_HYB_c);

ang_hyb = [phi_hyb'; theta_hyb'; psi_hyb'];


alt = ingps.ingps.ALT_HYB;
origin_imu = [lat_imu(1),lon_imu(1),alt(1)];
[xNorth_imu,yEast_imu,zDown_imu] = geodetic2ned(lat_imu,lon_imu,alt,lat_imu(1),lon_imu(1),alt(1),wgs84Ellipsoid,'degrees');
[xNorth_gps,yEast_gps,zDown_gps] = geodetic2ned(lat_gps,lon_gps,alt,lat_gps(1),lon_gps(1),alt(1),wgs84Ellipsoid,'degrees');
[xNorth_hyb,yEast_hyb,zDown_hyb] = geodetic2ned(lat_hyb,lon_hyb,alt,lat_hyb(1),lon_hyb(1),alt(1),wgs84Ellipsoid,'degrees');

gps_pose = [xNorth_gps';yEast_gps';zDown_gps'];

%plot(t(1:end-100),g_x_imu(1:end-100),'b',t(1:end-100),g_y_imu(1:end-100),'k',t(1:end-100),g_z_imu(1:end-100),'r');
zUp_imu = -zDown_imu;
zUp_gps = -zDown_gps;
zUp_hyb = -zDown_hyb;

%{
figure(1) % 3d trajectory
grid on
plot3(xNorth_imu(1:end-100),yEast_imu(1:end-100),zUp_imu(1:end-100),'b');                 
hold on
grid on
plot3(xNorth_gps(1:end-100),yEast_gps(1:end-100),zUp_gps(1:end-100),'k');
hold on
grid on
plot3(xNorth_hyb(1:end-100),yEast_hyb(1:end-100),zUp_hyb(1:end-100),'r');
legend('IMU', 'GPS', 'HYB')
title('3D trajectory')
xlabel('North (m)')
ylabel('East (m)')
zlabel('Up (m)')
print(figure(1),'3D-trajectory','-depsc');
%}

%{
figure(2)
grid on
subplot(3,1,1)
plot(t(1:end-100),vel_E_imu(1:end-100),'b',t(1:end-100),vel_E_GPS(1:end-100),'k');
legend('IMU', 'GPS')
title('Velocity East')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

subplot(3,1,2)
plot(t(1:end-100),vel_N_imu(1:end-100),'b',t(1:end-100),vel_N_GPS(1:end-100),'k');
legend('IMU', 'GPS')
title('Velocity North')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

subplot(3,1,3)
plot(t(1:end-100),vel_V_imu(1:end-100),'b',t(1:end-100),vel_V_GPS(1:end-100),'k');
legend('IMU', 'GPS')
title('Velocity Up')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
%}
end