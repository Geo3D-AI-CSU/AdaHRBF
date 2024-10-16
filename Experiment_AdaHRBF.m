% loading fault attribute points
attributes_fault = xlsread('Data\FaultChange2.xls');
% loading digital elevation model
DEM = xlsread('Data\DEM100m.xlsx');
% loading stratigraphic attribute points in two sub-domains
attributes_north = xlsread('Data\NStrataPoint.xls');
attributes_south = xlsread('Data\SStrataPoint.xls');
% loading stratigraphic attribute points in two sub-domains
attitude_north = xlsread('Data\NFattitude.xlsx');
attitude_south = xlsread('Data\SFattitude.xlsx');
% approximate minima in partial derivative calculation 
k=1e-4;
% set basefunction type(cubic function) 
basefunctype=1;

% geting initial gradient magnitude from attitude points
attitude_north=AttitudeToVector(attitude_north);
attitude_south=AttitudeToVector(attitude_south);
gradient_north=DipvectorToGradientvector(attitude_north);
gradient_south=DipvectorToGradientvector(attitude_south);
% geting the scope of study area and initial the mesh grid
dataMax_n=max([attributes_north(:,[1,2,3]);attitude_north(:,[1,2,3]);attributes_south(:,[1,2,3]);attitude_south(:,[1,2,3])]);
dataMin_n=min([attributes_north(:,[1,2,3]);attitude_north(:,[1,2,3]);attributes_south(:,[1,2,3]);attitude_south(:,[1,2,3])]);
dataMax_s=max([attributes_north(:,[1,2,3]);attitude_north(:,[1,2,3]);attributes_south(:,[1,2,3]);attitude_south(:,[1,2,3])]);
dataMin_s=min([attributes_north(:,[1,2,3]);attitude_north(:,[1,2,3]);attributes_south(:,[1,2,3]);attitude_south(:,[1,2,3])]);
[meshgrid_NX,meshgrid_NY,meshgrid_NZ] = meshgrid(dataMin_n(1,1):100:dataMax_n(1,1),dataMin_n(1,2):100:dataMax_n(1,2),0:25:dataMax_n(1,3));
[meshgrid_SX,meshgrid_SY,meshgrid_SZ] = meshgrid(dataMin_s(1,1):100:dataMax_s(1,1),dataMin_s(1,2):100:dataMax_s(1,2),0:25:dataMax_s(1,3));
%%
%HRBF: interpolating an initial 3D scalar field function in northern sub-domain
[alph_n,bravo_n,charlie_n]=GetParameters(1e-4,1,attributes_north,gradient_north);
valueGrid_north=GetValueGrid(k,basefunctype,meshgrid_NX,meshgrid_NY,meshgrid_NZ,attributes_north,gradient_north,alph_n,bravo_n,charlie_n);

%HRBF: interpolating an initial 3D scalar field function in southern sub-domain 
[alphs,betas,charlies]=GetParameters(1e-4,1,attributes_south,gradient_south);
valueGrid_south=GetValueGrid(k,basefunctype,meshgrid_SX,meshgrid_SY,meshgrid_SZ,attributes_south,gradient_south,alphs,betas,charlies);
%%
%adaHRBF: interpolating an optimal 3D scalar field in northern sub-domain
[opt_gradient_north,opt_alph_n,opt_bravo_n,opt_charlie_n]=OptGradientMagnitudes(attributes_north,gradient_north,100,1e-4,300);
opt_valueGrid_north=GetValueGrid(k,basefunctype,meshgrid_NX,meshgrid_NY,meshgrid_NZ,attributes_north,opt_gradient_north,opt_alph_n,opt_bravo_n,opt_charlie_n);

%adaHRBF: interpolating an optimal 3D scalar field in southern sub-domain
[opt_gradient_south,opt_alph_s,opt_bravo_s,opt_charlie_s]=OptGradientMagnitudes(attributes_south,gradient_south,100,1e-4,300);
opt_valueGrid_south=GetValueGrid(k,basefunctype,meshgrid_SX,meshgrid_SY,meshgrid_SZ,attributes_south,opt_gradient_south,opt_alph_s,opt_bravo_s,opt_charlie_s);
%%
%RBF: interpolating a 3D scalar field for fault
[alph,~,charlie]=GetParameters(1e-4,1,attributes_fault);
valueGrid_fault=GetValueGrid(k,basefunctype,meshgrid_NX,meshgrid_NY,meshgrid_NZ,attributes_fault,nan,alph,nan,charlie);
%%
% removing unnecessary grids
index_n=find(valueGrid_fault>=0);
index_s=find(valueGrid_fault<0);

valueGrid_north(index_n)=nan;
valueGrid_south(index_s)=nan;

opt_valueGrid_north(index_n)=nan;
opt_valueGrid_south(index_s)=nan;

index=find(meshgrid_NY<=-1.21*meshgrid_NX+3379598.04 | meshgrid_NY>=-1.21*meshgrid_NX+3389045.69 | meshgrid_NY>=0.83*meshgrid_NX+2028031.65 | meshgrid_NY<=0.83*meshgrid_NX+2011701.6);
valueGrid_north(index)=nan;
valueGrid_fault(index)=nan;
opt_valueGrid_north(index)=nan;
index=find(meshgrid_SY<=-1.21*meshgrid_SX+3379598.04 | meshgrid_SY>=-1.21*meshgrid_SX+3389045.69 | meshgrid_SY>=0.83*meshgrid_SX+2028031.65 | meshgrid_SY<=0.83*meshgrid_SX+2011701.6);
valueGrid_south(index)=nan;
valueGrid_fault(index)=nan;
opt_valueGrid_south(index)=nan;

% removing grids outside the study area in the northern subdomain
bar = waitbar(0,'Removing grids outside the northern study area...'); 
for i=1:size(meshgrid_NZ,1) 
    for j=1:size(meshgrid_NZ,2)
        for k=1:size(meshgrid_NZ,3)
            high = find(DEM(:,1)==meshgrid_NX(i,j,k) & DEM(:,2) == meshgrid_NY(i,j,k));
            if (meshgrid_NZ(i,j,k)>=DEM(high,3) | meshgrid_NZ(i,j,k)<0)
                valueGrid_north(i,j,k) = nan;
                opt_valueGrid_north(i,j,k) = nan;
            end
        end
    end
    waitbar(i/(size(meshgrid_NZ,1)),bar,'Removing grids outside the northern study area...')    
end
close(bar)
%removing grids outside the study area in the southern subdomain
bar = waitbar(0,'Removing grids outside the southern study area...'); 
for i=1:size(meshgrid_SZ,1)
    for j=1:size(meshgrid_SZ,2)
        for k=1:size(meshgrid_SZ,3)
            high = find(DEM(:,1)==meshgrid_SX(i,j,k) & DEM(:,2) == meshgrid_SY(i,j,k));
            if (meshgrid_SZ(i,j,k)>=DEM(high,3) | meshgrid_SZ(i,j,k)<0)
                valueGrid_south(i,j,k) = nan;
                opt_valueGrid_south(i,j,k) = nan;
            end
        end
    end
    waitbar(i/(size(meshgrid_SZ,1)),bar,'Removing grids outside the southern study area...')    
end
close(bar)
%%
% dsiplaying
value_list=[4872,7233,8429,9145,9674,10602,11717,8429,9145,9674,10602,11142];
color_list={'#ffe7ff','#fedfff','#ecffb0','#e2e2e2','#d0cfe1','#c5c4ca','#ffdfc0','#ecffb0','#e2e2e2','#d0cfe1','#c5c4ca','#ffe7cf'};
area=[662000 677000 2568000 2582000 0 1150];
color_axis=[4000,13000];

% dsiplaying scalar field in the northern subdomain
fig1=figure('Name','HRBFscatters');
figure(fig1);
% scalar field with initial gradient magnitude
subplot(1,2,1);
caxis(color_axis)
colorbar;
axis(area);
axis equal;
set(gca,'XLim',[area(1) area(2)]);
set(gca,'YLim',[area(3) area(4)]);
set(gca,'ZLim',[area(5) area(6)]);
scatter3(meshgrid_NX(:),meshgrid_NY(:),meshgrid_NZ(:),25,valueGrid_north(:),'.');
hold on;
% scalar field with optimized gradient magnitude
subplot(1,2,2);
caxis(color_axis)
colorbar;
axis(area);
axis equal;
set(gca,'XLim',[area(1) area(2)]);
set(gca,'YLim',[area(3) area(4)]);
set(gca,'ZLim',[area(5) area(6)]);
scatter3(meshgrid_NX(:),meshgrid_NY(:),meshgrid_NZ(:),25,opt_valueGrid_north(:),'.');
hold on;
%%
% displaying stratigraphic interfaces in the northern subdomain
fig2=figure('Name','HRBFisosurfaces');         
figure(fig2);
% stratigraphic interfaces with initial gradient magnitude
subplot(1,2,1);
light('position',[-0.5,-0.5,0.5],'style','local','color','[0.5,0.5,0.5]')
grid on; box on;
axis(area);
axis equal;
set(gca,'XLim',[area(1) area(2)]);
set(gca,'YLim',[area(3) area(4)]);
set(gca,'ZLim',[area(5) area(6)]);
for i=1:7
    pic = patch(isosurface(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_north,value_list(i)));
    isonormals(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_north,pic)
    set(pic,'FaceColor',color_list{i},'EdgeColor','none');
    hold on;
end
% dsiplaying faults interfaces
pic_f = patch(isosurface(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_fault,0));
isonormals(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_fault,pic_f)
set(pic_f,'FaceColor','red','EdgeColor','none');

% stratigraphic interfaces with optimized gradient magnitude
subplot(1,2,2);
light('position',[-0.5,-0.5,0.5],'style','local','color','[0.5,0.5,0.5]')
grid on; box on;
axis(area);
axis equal;
set(gca,'XLim',[area(1) area(2)]);
set(gca,'YLim',[area(3) area(4)]);
set(gca,'ZLim',[area(5) area(6)]);
for i=1:7
    pic = patch(isosurface(meshgrid_NX,meshgrid_NY,meshgrid_NZ,opt_valueGrid_north,value_list(i)));
    isonormals(meshgrid_NX,meshgrid_NY,meshgrid_NZ,opt_valueGrid_north,pic)
    set(pic,'FaceColor',color_list{i},'EdgeColor','none');
    hold on;
end
% dsiplaying faults interfaces
pic_f = patch(isosurface(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_fault,0));
isonormals(meshgrid_NX,meshgrid_NY,meshgrid_NZ,valueGrid_fault,pic_f)
set(pic_f,'FaceColor','red','EdgeColor','none');

%%
% dsiplaying scalar field in the southern subdomain
figure(fig1);
% scalar field with initial gradient magnitude
subplot(1,2,1);
scatter3(meshgrid_SX(:),meshgrid_SY(:),meshgrid_SZ(:),25,valueGrid_south(:),'.');
% scalar field with optimized gradient magnitude
subplot(1,2,2);
scatter3(meshgrid_SX(:),meshgrid_SY(:),meshgrid_SZ(:),25,opt_valueGrid_south(:),'.');
hold off;
%%
% dsiplaying stratigraphic interfaces in the southern subdomain
figure(fig2);
% stratigraphic interfaces with initial gradient magnitude
subplot(1,2,1);
for i=8:12
    pic = patch(isosurface(meshgrid_SX,meshgrid_SY,meshgrid_SZ,valueGrid_south,value_list(i)));
    isonormals(meshgrid_SX,meshgrid_SY,meshgrid_SZ,valueGrid_south,pic)
    set(pic,'FaceColor',color_list{i},'EdgeColor','none');
    hold on;
end
%stratigraphic interfaces with optimized gradient magnitude
subplot(1,2,2);
for i=8:12
    pic = patch(isosurface(meshgrid_SX,meshgrid_SY,meshgrid_SZ,opt_valueGrid_south,value_list(i)));
    isonormals(meshgrid_SX,meshgrid_SY,meshgrid_SZ,opt_valueGrid_south,pic)
    set(pic,'FaceColor',color_list{i},'EdgeColor','none');
    hold on;
end
hold off;




 
