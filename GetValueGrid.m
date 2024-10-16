% name:GetValueGrid.m
% description:interpolate a 3D scalar field function in southern sub-domain 
% input:k float approximate minima in calculation(useless in current version) 
%         baseFuncType int set basefunction type(useless in current version) 
%         grid_X meshgrid X
%         grid_Y meshgrid Y
%         grid_Z meshgrid Z
%         attribute_points[x,y,z,attribute] relatively butial depth for RBF,(ada)HRBF
%         gradient_points [x,y,z,gx,gy,gz] gradient magnitude for (ada)HRBF
%         varargin{1} weight coefficients a
%         varargin{2} weight coefficients b
%         varargin{3} charlie weight coefficients c
% return:valueGrid
% author:Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.��
% version:ver1.0.0[2022.9.23]
function [valueGrid]=GetValueGrid(k,baseFuncType,grid_X,grid_Y,grid_Z,attribute_points,gradient_points,varargin)

if  isnan (gradient_points)
    %RBF
    alph = varargin{1};
    charlie= varargin{3};
    %number of gradient points
    m=0;
    clear varargin;
else
    %HRBF
    alph = varargin{1};
    bravo=varargin{2};
    charlie= varargin{3};
    m=size(gradient_points,1);
    clear varargin;
end

bar = waitbar(0,'interpolating...'); % waitbar 

num=size(attribute_points,1);%%number of attribute points
firstPart=zeros(size(grid_X));%first monomial of interpolation function
for i=1:num
   Xi=ones(size(grid_X)).*attribute_points(i,1);
   Yi=ones(size(grid_Y)).*attribute_points(i,2);
   Zi=ones(size(grid_Z)).*attribute_points(i,3);
   baseFunc=BaseFunc(grid_X,grid_Y,grid_Z,Xi,Yi,Zi,1);
   firstPart=firstPart+alph(i,1).*baseFunc;
   waitbar(i/(num+m),bar,'interpolating...')    
end

%second monomial of interpolation function
secondPart=ones(size(grid_X))*charlie(1,1)+charlie(2,1)*grid_X+charlie(3,1)*grid_Y+charlie(4,1)*grid_Z;

%RBF
if isnan (gradient_points)
    valueGrid=-(firstPart+secondPart);
    close(bar)
    return
end

%thrid monomial of interpolation function
thridPart=zeros(size(grid_X));
for i=1:m
    Xi=ones(size(grid_X))*gradient_points(i,1);
    Yi=ones(size(grid_Y))*gradient_points(i,2);
    Zi=ones(size(grid_Z))*gradient_points(i,3);
    [delta_TX,delta_TY,delta_TZ]=GetPartialDerivative(grid_X,grid_Y,grid_Z,Xi,Yi,Zi,k,baseFuncType);
    betaUsed=bravo(3*(i-1)+(1:3),1);
    thridPart=thridPart+(betaUsed(1,1)*delta_TX+betaUsed(2,1)*delta_TY+betaUsed(3,1)*delta_TZ);
    waitbar((i+num)/(num+m),bar,'interpolating...')  
end
%HRBF
valueGrid=-(firstPart+thridPart+secondPart);
close(bar)
end
