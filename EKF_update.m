function [X,P] = EKF_update(X,P,Z,H,h,R,sensor)


K = P*H'*inv(H*P*H' + R);
P = (eye(6) - K*H)*P;
X = X + K*(Z - h);
X = ang_rot(X);


end