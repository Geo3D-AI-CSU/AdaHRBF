% name：GetParameters.m
% description：get weight coefficients and polynomial coefficients
% input：k（float）：approximate minima in calculation(useless in current version) 
%         baseFuncType（int）：set basefunction type(useless in current version) 
%         grid_X：meshgrid X
%         grid_Y：meshgrid Y
%         grid_Z：meshgrid Z
%         varargin{1}：[x,y,z,attribute] relatively butial depth for RBF or (ada)HRBF
%         varargin{2}：[x,y,z,gx,gy,gz] gradient magnitude for (ada)HRBF
% return：alph：weight coefficients a
%         charlie：weight coefficients c
%         bravo：weight coefficients b
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.   
% version：ver1.0.0[2022.9.23]

function [alph,bravo,charlie]=GetParameters(k,baseFuncType,varargin)
% getting number of attribute points。
r_num=size(varargin{1},1);
% RBF
if nargin == 3
% constructing the RBF linear system
    A=GetMatrixA(k,baseFuncType,varargin{1});
    B=GetMatrixB(varargin{1});
    X =lsqminnorm(A,B);
    alph=X(1:r_num,1);
    charlie=X(r_num+1:r_num+4,1);
    bravo=nan;
    return
end
% get number of attribute_points。
r_m=size(varargin{2},1);

% HRBF
if nargin == 4
%  constructing the HRBF linear system
    A=GetMatrixA(k,baseFuncType,varargin{1},varargin{2});
    B=GetMatrixB(varargin{1},varargin{2});
    X =lsqminnorm(A,B);
    alph=X(1:r_num,1);
    bravo=X((r_num+1):(r_num+3*r_m),1);
    charlie=X((r_num+3*r_m+1):(r_num+3*r_m+4),1);
    return
end
error('HRBFGetPara:The number of parameters is incorrect.')
end




