%测试lyapunov_wolf法

global delt r b; 
delt=10;
b=8/3;
lya_pu=[];
sup_delt=1e-5;

r_test=[5:0.1:28];
lya_pu=[];
tic
for i=1:size(r_test,2)
r=r_test(i)
% r=5
x0=1;y0=1;z0=1;
sim('lorenz_sim.mdl',50)
data=x(3000:end);
lya=lyapunov(data,length(data),8,12,40)
lya_pu=[lya_pu,lya];
toc
end
figure
plot(r_test,lya_pu,'k');
xlabel('r')
ylabel('李雅普诺夫指数')