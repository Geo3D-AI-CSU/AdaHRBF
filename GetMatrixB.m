% name：GetMatrixA.m
% description：get the divisor of coefficients solving linear system
% input：k（float）：approximate minima in calculation(useless in current version) 
%         baseFuncType（int）：set basefunction type(useless in current version) 
%         varargin{1}：[x,y,z,attribute] relatively butial depth for RBF,(ada)HRBF
%         varargin{2}：[x,y,z,gx,gy,gz] gradient magnitude for (ada)HRBF
% return：matrix B 
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.
% version：ver1.0.0[2022.9.23]
function B = GetMatrixB(varargin)

if nargin == 1
%RBF
    attribute_points = varargin{1};
    clear varargin;
elseif nargin == 2
%HRBF
    attribute_points = varargin{1};
    gradient_points = varargin{2};
    r_m=size(gradient_points,1);
    clear varargin;
else
   error('GetMatrixB:Please enter the correct number of parameters.')
end

r_num= size(attribute_points,1);
%B（1，1）
foxtra=attribute_points(1:r_num,4);

if nargin == 1
%construct matrix B in RBF
    B = zeros(r_num+4,1);
    B(1:r_num,1)=foxtra;
    return
end
%B（2，1）
golf=gradient_points(1:r_m,4:6)';

if nargin ==2
%construct matrix B in RBF
    B = zeros(r_num+3*r_m+4,1);
    B(1:r_num,1)=foxtra;
    B((r_num+1):(r_num+3*r_m),1)=golf(1:end);
    return
end
end
