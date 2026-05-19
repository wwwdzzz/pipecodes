
% close all            
% clear                
% clc                    
p_test=[];
p_train=[];
t_train=[];
p_test=[];
%%  导入数据
load LSTMdata.mat;

% res1=LSTMtraindata;
% res2=LSTMdata2;
res3=LSTMdata3;
% res4=LSTMdata4;
% res5=LSTMdata_simple;
%res5=LSTMdata5;
res=res3;
%res=[res1(1:end-1,1),res1(2:end,1),res1(2:end,2),res1(2:end,3)];
%res=[res2(1:end-1,1),res2(2:end,1),res2(2:end,2),res2(2:end,3)-res2(2:end,1)];%xi xi+1 u xi+2-xi+1;
% res=[res2(2:end,1),res2(2:end,2),res2(2:end,3)-res2(2:end,1)];

datasize=size(res,1);
temp = randperm(datasize);
templ=linspace(1,datasize,datasize);
trainsize=floor(datasize*0.98);

P_train = res(templ(1: trainsize), [1 2 ])';
T_train = res(templ(1: trainsize), 3)';
M = size(P_train, 2);

P_test = res(templ(trainsize: end), [1 2  ])';
T_test = res(templ(trainsize: end), 3)';
N = size(P_test, 2);
datadim=size(P_test, 1);


[P_train, ps_input] = mapminmax(P_train, 0, 1);
P_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

P_train =  double(reshape(P_train, datadim, 1, 1, M));%
P_test  =  double(reshape(P_test , datadim, 1, 1, N));

t_train = t_train';
t_test  = t_test' ;

for i = 1 : M
    p_train{i, 1} = P_train(:, :, 1, i);
end

for i = 1 : N
    p_test{i, 1}  = P_test( :, :, 1, i);
end


layers = [
    sequenceInputLayer(datadim)               % 输入层
    
    lstmLayer(4, 'OutputMode', 'last')  % LSTM层
    reluLayer                           % Relu激活层
    
    fullyConnectedLayer(1)              % 全连接层
    regressionLayer];                   % 回归层
 

options = trainingOptions('adam', ...      
    'MaxEpochs', 50, ...                 
    'InitialLearnRate', 0.01, ...          
    'LearnRateSchedule', 'piecewise', ...  
    'LearnRateDropFactor', 0.1, ...        
    'LearnRateDropPeriod', 1200, ...       
    'Shuffle', 'every-epoch', ...        
    'Plots', 'training-progress', ...     
    'Verbose', false);

net = trainNetwork(p_train, t_train, layers, options);
%% 


t_sim1 = predict(net, p_train);
t_sim2 = predict(net, p_test );

T_sim1 = mapminmax('reverse', t_sim1, ps_output);
T_sim2 = mapminmax('reverse', t_sim2, ps_output);

error1 = sqrt(sum((T_sim1' - T_train).^2) ./ M);
error2 = sqrt(sum((T_sim2' - T_test ).^2) ./ N);

analyzeNetwork(net)

figure
plot(1: M, T_train, 'r-*', 1: M, T_sim1, 'b-o', 'LineWidth', 1)
legend('真实值', '预测值')
xlabel('预测样本')
ylabel('预测结果')
string = {'训练集预测结果对比'; ['RMSE=' num2str(error1)]};
title(string)
xlim([1, M])
grid

figure
plot(1: N, T_test, 'r-*',1:N, T_sim2, 'b-o', 'LineWidth', 1)
legend('真实值', '预测值')
xlabel('预测样本')
ylabel('预测结果')
string = {'测试集预测结果对比'; ['RMSE=' num2str(error2)]};
title(string)
% xlim([1, N])
grid

%%  相关指标计算
% R2
R1 = 1 - norm(T_train - T_sim1')^2 / norm(T_train - mean(T_train))^2;
R2 = 1 - norm(T_test  - T_sim2')^2 / norm(T_test  - mean(T_test ))^2;

disp(['训练集数据的R2为：', num2str(R1)])
disp(['测试集数据的R2为：', num2str(R2)])

% MAE
mae1 = sum(abs(T_sim1' - T_train)) ./ M ;
mae2 = sum(abs(T_sim2' - T_test )) ./ N ;

disp(['训练集数据的MAE为：', num2str(mae1)])
disp(['测试集数据的MAE为：', num2str(mae2)])

% MBE
mbe1 = sum(T_sim1' - T_train) ./ M ;
mbe2 = sum(T_sim2' - T_test ) ./ N ;

disp(['训练集数据的MBE为：', num2str(mbe1)])
disp(['测试集数据的MBE为：', num2str(mbe2)])

%%  绘制散点图
sz = 25;
c = 'b';

figure
scatter(T_train, T_sim1, sz, c)
hold on
plot(xlim, ylim, '--k')
xlabel('训练集真实值');
ylabel('训练集预测值');
xlim([min(T_train) max(T_train)])
ylim([min(T_sim1) max(T_sim1)])
title('训练集预测值 vs. 训练集真实值')

figure
scatter(T_test(1400:1800), T_sim2(1400:1800), sz, c)
hold on
 plot(xlim, ylim, '--k')
xlabel('测试集真实值');
ylabel('测试集预测值');
% xlim([min(T_test) max(T_test)])
% ylim([min(T_sim2) max(T_sim2)])
title('测试集预测值 vs. 测试集真实值')

save("LSTMsimplenet.mat",'net')
save("maxminout.mat", 'ps_output')
save("maxminin.mat", 'ps_input')