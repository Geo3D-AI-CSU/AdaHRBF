% name：AttitudToVector.m
% description：get the vector of attitude
% input：attitude[x,y,z,attribute] 
% return：vector[x,y,z,dx,dy,dz]
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.  
% version：ver1.0.0[2022.9.23]

function [vector]=AttitudeToVector(attitude)
    dipAngl=deg2rad(attitude(:,5));
    dipOrin=deg2rad(attitude(:,4));
    dipVct=[cos(dipAngl).*sin(dipOrin),cos(dipAngl).*cos(dipOrin),-sin(dipAngl)];
    vector=[attitude(:,1:3),dipVct];
end