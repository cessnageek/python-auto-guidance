close all; clear all; clc;

V2 = 2.3;%2.3; %m/s
l = 1.72; %m
d = 0.4; %m
kX = -0.13;
kPhi = 0.6;
kD = 0*0.01;
KdPhi = 0.2*0.01;%0.1;
setPoint = 0;

stepperAccel = 1.5; %rev/s^2
initSteeringAngle = 0; %rad
initSteeringVelocity = 0; %rad/s
xArray = [[]];
timeArray = [[]];
sim('Tractor_Model_Feedforward.slx');
% for x = 0:0.1:1.7
%     d = x;
%     sim('Tractor_Model.slx');
%     xArray(end+1,:) = ans.x;
%     timeArray(end+1,:) = ans.time;
% end