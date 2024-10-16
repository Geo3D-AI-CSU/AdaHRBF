% name：GetSecPartialDerivative.m
% description：calculating bsaefuncton value
% input：k（float）：approximate minima in calculation(useless in current version) 
%         points[X,Y,Z]
%         X_：constants points
%         Y_：constants points
%         Z_：constants points
%         type（int）：set basefunction type(useless in current version)
% return：secPartialDerivativ(3*3)
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.
% version：ver1.0.0[2022.9.23]
function secPartialDerivative=GetSecPartialDerivative(x,y,z,x_,y_,z_,k,baseFuncType)

grandg_xx = @(b1,b2,b3,x,y,z)  3*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(1/2) + 3*((b1-x)^2)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);
grandg_yy = @(b1,b2,b3,x,y,z)  3*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(1/2) + 3*((b2-y)^2)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);
grandg_zz = @(b1,b2,b3,x,y,z)  3*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(1/2) + 3*((b3-z)^2)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);
grandg_xy = @(b1,b2,b3,x,y,z)  3*(b1-x)*(b2-y)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);
grandg_xz = @(b1,b2,b3,x,y,z)  3*(b1-x)*(b3-z)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);
grandg_yz = @(b1,b2,b3,x,y,z)  3*(b2-y)*(b3-z)*((b1-x)^2 + (b2-y)^2 + (b3-z)^2)^(-1/2);

e3(1,1)=grandg_xx(x,y,z,x_,y_,z_);
e3(1,2)=grandg_xy(x,y,z,x_,y_,z_);
e3(1,3)=grandg_xz(x,y,z,x_,y_,z_);
e3(2,1)=grandg_xy(x,y,z,x_,y_,z_);
e3(2,2)=grandg_yy(x,y,z,x_,y_,z_);
e3(2,3)=grandg_yz(x,y,z,x_,y_,z_);
e3(3,1)=grandg_xz(x,y,z,x_,y_,z_);
e3(3,2)=grandg_yz(x,y,z,x_,y_,z_);
e3(3,3)=grandg_zz(x,y,z,x_,y_,z_);

secPartialDerivative=e3;
% 
% k=k*(1e-2);
% basefunValue=BaseFunc(x,y,z,x_,y_,z_,baseFuncType);
% basefunValue_xk=BaseFunc(x+k,y,z,x_,y_,z_,baseFuncType);
% basefunValue_yk=BaseFunc(x,y+k,z,x_,y_,z_,baseFuncType);
% basefunValue_zk=BaseFunc(x,y,z+k,x_,y_,z_,baseFuncType);
% a=(BaseFunc(x+2.*k,y,z,x_,y_,z_,baseFuncType)-2.*basefunValue_xk+basefunValue);
% grandg_xx=(BaseFunc(x+2.*k,y,z,x_,y_,z_,baseFuncType)-2.*basefunValue_xk+basefunValue)./(k^2);
% grandg_yy=(BaseFunc(x,y+2.*k,z,x_,y_,z_,baseFuncType)-2.*basefunValue_yk+basefunValue)./(k^2);
% grandg_zz=(BaseFunc(x,y,z+2.*k,x_,y_,z_,baseFuncType)-2.*basefunValue_zk+basefunValue)./(k^2);
% grandg_xy=(BaseFunc(x+k,y+k,z,x_,y_,z_,baseFuncType)-basefunValue_xk-basefunValue_yk+basefunValue)./(k^2);
% grandg_xz=(BaseFunc(x+k,y,z+k,x_,y_,z_,baseFuncType)-basefunValue_xk-basefunValue_zk+basefunValue)./(k^2);
% grandg_yz=(BaseFunc(x,y+k,z+k,x_,y_,z_,baseFuncType)-basefunValue_yk-basefunValue_zk+basefunValue)./(k^2);
% secPartialDerivative=[grandg_xx,grandg_xy,grandg_xz;grandg_xy,grandg_yy,grandg_yz;grandg_xz,grandg_yz,grandg_zz];

end