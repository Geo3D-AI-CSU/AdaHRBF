% name：GetPartialDerivative.m
% description：calculating bsaefuncton value
% input：k（float）：approximate minima in calculation(useless in current version) 
%         meshgrid [X,Y,Z]
%         X_：constants points
%         Y_：constants points
%         Z_：constants points
%         type（int）：set basefunction type(useless in current version)
% return：x_partialDerivative,y_partialDerivative,z_partialDerivative
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.
% version：ver1.0.0[2022.9.23]
function [x_partialDerivative,y_partialDerivative,z_partialDerivative]=GetPartialDerivative(X,Y,Z,X_,Y_,Z_,k,baseFuncType)

grandg_x = @(b1,b2,b3,x,y,z)  (3*(2*b1 - 2*x).*((b1 - x).^2 + (b2 - y).^2 + (b3 - z).^2).^(1/2))./2;
grandg_y = @(b1,b2,b3,x,y,z)  (3*(2*b2 - 2*y).*((b1 - x).^2 + (b2 - y).^2 + (b3 - z).^2).^(1/2))./2;
grandg_z = @(b1,b2,b3,x,y,z)  (3*(2*b3 - 2*z).*((b1 - x).^2 + (b2 - y).^2 + (b3 - z).^2).^(1/2))./2;

x_partialDerivative=grandg_x(X,Y,Z,X_,Y_,Z_);
y_partialDerivative=grandg_y(X,Y,Z,X_,Y_,Z_);
z_partialDerivative=grandg_z(X,Y,Z,X_,Y_,Z_);

% basefunValue=BaseFunc(X,Y,Z,X_,Y_,Z_,baseFuncType);
% x_partialDerivative=(basefunValue-BaseFunc(X-k,Y,Z,X_,Y_,Z_,baseFuncType))./k;
% y_partialDerivative=(basefunValue-BaseFunc(X,Y-k,Z,X_,Y_,Z_,baseFuncType))./k;
% z_partialDerivative=(basefunValue-BaseFunc(X,Y,Z-k,X_,Y_,Z_,baseFuncType))./k;

end