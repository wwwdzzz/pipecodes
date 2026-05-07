%²âÊÔlyapunov_wolf·¨
clear
global delt r b; 
delt=10;
b=8/3;
lya_pu=[];
sup_delt=1e-5;

r_test=[5:0.1:28];
lya_pu=[];
tic
r=114
% r=5
x0=1;y0=1;z0=1;
sim('lorenz_sim.mdl',50)
data1=x(3000:end);
%data1=data1+rand(length(data1),1)*0.1;
data2=y(3000:end);
data3=z(3000:end);
lya1=lyapunov(data1,length(data1),4,12,40)
% lya2=lyapunov(data2,length(data2),11,19,40)
% lya3=lyapunov(data3,length(data3),14,1,40)
toc
