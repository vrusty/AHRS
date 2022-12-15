function [X,P] = EKF_prediction(X,P,F_nlin,A,G,Q,dt)


Tau = G*dt;
P = A*P*A' + Tau*Q*Tau';
X = X + F_nlin*dt;
X = ang_rot(X);


end