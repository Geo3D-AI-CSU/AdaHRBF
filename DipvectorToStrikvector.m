% name：DipvectorToStrikvector .m
% description：get the strike vector from attitude vector
% input：attitude_vector[x,y,z,dx,dy,dz]
% return：strike_vector[x,y,z,gx,gy,gz]
% author：Linze Du, Yongqiang Tong, Baoyi Zhang, Umair Khan, Lifang Wang and Hao Deng.  
% version：ver1.0.0[2022.9.23]

function strike_vector= DipvectorToStrikvector(attitude_vector)

    dipOrin=attitude_vector(:,4:6);
    dipOrin(:,end)=0;
    stkOrin=cross(dipOrin,attitude_vector(:,4:6));
    Norm=sqrt(stkOrin(:,1).^2+stkOrin(:,2).^2);
    stkNorm=[Norm,Norm,Norm];
    stkOrin=stkOrin./stkNorm;
    strike_vector=[attitude_vector(:,1:3),stkOrin];

end