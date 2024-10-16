% name：GetMatrixA.m
% description：get the dividend of coefficients solving linear system
% input：k（float）：approximate minima in calculation(useless in current version) 
%         baseFuncType（int）：set basefunction type(useless in current version) 
%         varargin{1}：[x,y,z,attribute]relatively butial depth for RBF,(ada)HRBF
%         varargin{2}：[x,y,z,gx,gy,gz]gradient magnitude for (ada)HRBF
% return：matrix A 
% author：
% version：ver1.0.0[2022.9.23]

function A = GetMatrixA(k,baseFuncType,varargin)
if nargin == 3
% RBF
    attribute_points = varargin{1};
    clear varargin;
elseif nargin == 4
% HRBF
    attribute_points = varargin{1};
    gradient_points = varargin{2};
    r_m=size(gradient_points,1);
    clear varargin;
else
   error('GetMatrixA:The number of parameters is incorrect.')
end

r_num= size(attribute_points,1);

%A（1，1）
Phi=zeros(r_num,r_num);
for i=1:r_num
    xi=ones(r_num,1)*attribute_points(i,1);
    yi=ones(r_num,1)*attribute_points(i,2);
    zi=ones(r_num,1)*attribute_points(i,3);
    Phi(i,1:r_num)=(BaseFunc(attribute_points(:,1),attribute_points(:,2),attribute_points(:,3),xi,yi,zi,baseFuncType))';
end

%A（1，3）
Charlie=ones(r_num,4);
Charlie(:,2:4)=attribute_points(1:r_num,1:3);

if nargin == 3
%construct matrix A in RBF
    A = zeros(r_num+4,r_num+4);
    A(1:r_num,1:r_num)=Phi;
    A(1:r_num,r_num+1:r_num+4)=Charlie;
    A(r_num+1:r_num+4,1:r_num)=Charlie';
    return
end

%A（1，2）
Nabla_Phi=zeros(r_num,3*r_m);
for i=1:r_m    
    xi=ones(r_num,1)*gradient_points(i,1);
    yi=ones(r_num,1)*gradient_points(i,2);
    zi=ones(r_num,1)*gradient_points(i,3);
    [Nabla_Phi(:,3*i-2),Nabla_Phi(:,3*i-1),Nabla_Phi(:,3*i)]=GetPartialDerivative(attribute_points(:,1),attribute_points(:,2),attribute_points(:,3), xi,yi,zi,k,baseFuncType);
end

%A（2，1）
Nabla_Phi_T=zeros(3*r_m,r_num);
for i=1:r_m    
    xi=ones(1,r_num)*gradient_points(i,1);
    yi=ones(1,r_num)*gradient_points(i,2);
    zi=ones(1,r_num)*gradient_points(i,3);
    [Nabla_Phi_T(3*i-2,:),Nabla_Phi_T(3*i-1,:),Nabla_Phi_T(3*i,:)]=GetPartialDerivative( xi,yi,zi,attribute_points(:,1)',attribute_points(:,2)',attribute_points(:,3)',k,baseFuncType);
end

%A（2，2）
Nabla2_Phi_T=zeros(3*r_m,3*r_m);
for i=1:r_m
    for j=1:r_m
        if i~=j
            e3=GetSecPartialDerivative(gradient_points(i,1),gradient_points(i,2),gradient_points(i,3),gradient_points(j,1),gradient_points(j,2),gradient_points(j,3),0.01,1);
        else
            e3=[0,0,0;0,0,0;0,0,0];
        end 
        Nabla2_Phi_T(3*(i-1)+(1:3),3*(j-1)+(1:3))=double(e3);
    end
end
clear e3;

%A（2，3）
Nabla_Charlie=zeros(3*r_m,4);
for i=1:3*r_m
    one_pos=mod(i-1,3);
    Nabla_Charlie(i,end+one_pos-2)=1;
end

if nargin == 4
% construct matrix A in (ada）HRBF
    A = zeros(r_num+3*r_m+4,r_num+3*r_m+4);
    A(1:r_num,1:r_num)=Phi;
    A(1:r_num,r_num+1:r_num+3*r_m)=Nabla_Phi;
    A(1:r_num,r_num+3*r_m+1:r_num+3*r_m+4)=Charlie;
    A((r_num+1):(r_num+3*r_m),1:r_num)=Nabla_Phi_T;
    A((r_num+1):(r_num+3*r_m),(r_num+1):(r_num+3*r_m))= Nabla2_Phi_T;
    A((r_num+1):(r_num+3*r_m),(r_num+3*r_m+1):(r_num+3*r_m+4))=Nabla_Charlie;
    A((r_num+3*r_m+1):(r_num+3*r_m+4),1:r_num)=Charlie';
    A((r_num+3*r_m+1):(r_num+3*r_m+4),(r_num+1):(r_num+3*r_m))=Nabla_Charlie';
    return
end

end

