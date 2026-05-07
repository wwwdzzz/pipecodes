function [ln_r,ln_C]=G_P(data,N,tau,min_m,max_m,ss)
% the function is used to calculate correlation dimention with G-P algorithm
%    计算关联维数的G－P算法
% data:the time series                       时间序列
% N: the length of the time series           时间序列长度
% tau: the time delay                        时间延迟
% min_m:the least embedded dimention m       最小的嵌入维数
% max_m:the largest embedded dimention m     最大的嵌入维数
% ss:the stepsize of r                       r的步长
%skyhawk
for m=min_m:max_m
    Y=reconstitution(data,N,m,tau);%reconstitute state space
    M=N-(m-1)*tau;%the number of points in state space
    for i=1:M-1
        for j=i+1:M
            d(i,j)=max(abs(Y(:,i)-Y(:,j)));%calculate the distance of each two           
        end                                %points in state space  计算状态空间中每两点之间的距离
    end
    max_d=max(max(d));%the max distance of all points   得到所有点之间的最大距离
    d(1,1)=max_d;
    min_d=min(min(d));%the min distance of all points   得到所有点间的最短距离
    delt=(max_d-min_d)/ss;%the stepsize of r            得到r的步长
    for k=1:ss
        r=min_d+k*delt;
        C(k)=correlation_integral(Y,M,r);%calculate the correlation integral
        ln_C(m,k)=log(C(k));%lnC(r)
        ln_r(m,k)=log(r);%lnr
        fprintf('%d/%d/%d/%d\n',k,ss,m,max_m);
    end

end
fid=fopen('lnr.txt','w');
fprintf(fid,'%6.2f %6.2f\n',ln_r);
fclose(fid);
fid = fopen('lnC.txt','w');
fprintf(fid,'%6.2f %6.2f\n',ln_C);
fclose(fid);

