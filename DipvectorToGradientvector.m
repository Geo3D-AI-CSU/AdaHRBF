% name：DipvectorToGradientvct.m
% description：get the gradient vector from attitude vector
% input：k（float）：attitude_vector[x,y,z,dx,dy,dz]
% return：gradient_vector[x,y,z,gx,gy,gz]
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.  
% version：ver1.0.0[2022.9.23]
function [gradient_vector]=DipvectorToGradientvector(attitude_vector)

    stkOrin=DipvectorToStrikvector(attitude_vector);%走向矢量
    gradient_vector=attitude_vector;
    gradient_vector(:,4:6)=cross(attitude_vector(:,4:6),stkOrin(:,4:6));%叉积求得梯度矢量
   
end
    




