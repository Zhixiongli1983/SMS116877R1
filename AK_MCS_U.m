% SNS+主动学习Kriging程序
clc;
clear all;

% 将全部文件夹加入路径
addpath(genpath(pwd));

N_target=1e8;
load data_S_lhs_12_4var
load data_Y_lhs_12_4var
load MCS_Variations_01
load MCS_Lifetime_01
load MCS_Variations_03
load MCS_Lifetime_03

xn=[MCS_Variations_01;MCS_Variations_03];
mcs=[MCS_Lifetime_01;MCS_Lifetime_03];
% xn=MCS_Variations_01;
% mcs=MCS_Lifetime_01;
Y=Y-N_target;
mcs=mcs-N_target;

% 边界条件及载荷设置
% current=S(:,1);
N_MCS=10000; 
load MCS_Variations

history_points=[];
%初始化参数
for i=1:2000%%
theta = [10,10,10,10]; lob = [1e-7,1e-7,1e-7,1e-7]; upb = [2000,2000,2000,2000];
[dmodel, perf] = dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb);
% theta = [10,10,10,10]; lob = [1e-1,1e-1,1e-1,1e-1]; upb = [20,20,20,20];
% [dmodel, perf] = dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb);
%%指候选样本集，与下一步中的X表示的不一样
[YX MSE] = predictor(xn, dmodel);%%X指候选样本集中除去参与构件Kriging模型外的其余样本点
%% LearningFunction
UC=abs((YX)./sqrt(MSE));
   if min(UC)>2
       break
   end
   [a,b]=find(UC==min(UC));
   
    % 判断最佳点是否重合

   addpoint=find(history_points==a);
   p= isempty(addpoint);
   
   while p==0
       UC(a)=[];
       [a,b]=find(UC==min(UC));
       addpoint1=find(history_points==a); 
       p= isempty(addpoint1);
   end

   history_points(i,1)=a;   % 历史最佳点编号
   U_min(i,1)=UC(a);
   
   Xgengxin=xn(a,:);          %添加主动学习样本点
   Ygengxin = mcs(a,1);

   S=[S;Xgengxin];
   Y=[Y;Ygengxin];
   
    % 计算失效概率
   t=0;
   for n=1:N_MCS
       life_mcs=predictor(MCS_Variations(n,:), dmodel);
       if life_mcs<=0
           t=t+1;
       end
   end
   pf(i,1)=t/N_MCS;
end
