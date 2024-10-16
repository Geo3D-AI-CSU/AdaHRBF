
% name：BsaeFunc.m
% description：calculating bsaefuncton value
% input：k（float）：approximate minima in calculation(useless in current version) 
%         meshgrid [X,Y,Z]
%         X_：constants points
%         Y_：constants points
%         Z_：constants points
%         type（int）：set basefunction type(useless in current version)
% return：basefunction value
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.
% version：ver1.0.0[2022.9.23]
function basefunc=BaseFunc(X,Y,Z,X_,Y_,Z_,type)

basefunc=((sqrt((X-X_).^2 + (Y-Y_).^2 +(Z-Z_).^2 )).^3);
% Gaussian
% MQ
% IMQ
end