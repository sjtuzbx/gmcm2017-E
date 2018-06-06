%考虑整体暴露时间尽可能短，也要规避敌方的侦察和打击，采用适当分散机动的策略，同时还要缩短单台发射装置的最长暴露时间
%% 这是重构代码的部分
% 为了便于后面两问的结果
% 数据预处理
clear; clc; close all;
n = 2 + 6 + 60 + 62;

load Point
load edge

numOfEdges = 175; numOfF = 60; numOfPoints = size(Point, 1);
velocity_Amain = 70; % 70km/h
velocity_Bmain = 60;
velocity_Cmain = 50;

velocity_Aother = 45; % 45km/h
velocity_Bother = 35;
velocity_Cother = 30;

index_main = [175:178, 181,182,189,190,195,196,201,202,207,208,211,212,217,218,221,222,227,228,231,232,237,238,243,244,249,250,253:256,261,262];
type = zeros(2*numOfEdges, 1);
type(index_main, 1) = 1;

length = zeros(2*numOfEdges, 1);
velocity = zeros(2*numOfEdges, 3);
for i=1:2*numOfEdges
    length(i) = calculateDistance(Point(edge(i,1), :), Point(edge(i,2), :));
    velocity(i, 1) = type(i,1) * velocity_Amain + (1 - type(i, 1)) * velocity_Aother;
    velocity(i, 2) = type(i,1) * velocity_Bmain + (1 - type(i, 1)) * velocity_Bother;
    velocity(i, 3) = type(i,1) * velocity_Cmain + (1 - type(i, 1)) * velocity_Cother;
end

out = cell(1, numOfPoints);
in = cell(1, numOfPoints);
for i=1:2*numOfEdges
    temp = edge(i, :);
    out = addEdge(out, temp(1), i);
    in = addEdge(in, temp(2), i);
end

%% 整数线性规划方程组生成
[f, intcon, A, b, Aeq, beq, lb, ub] = generatePathParameters(numOfEdges, numOfF, numOfPoints,edge, length, velocity, in, out, []);
[path_x, path_val] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);

    
%% Log information
[result_path, path_id, edge_id, stagePoint] = savePathInfo(path_x, in, out, edge, [6,6,12], F, Z, numOfEdges);
[time_table, timeDF, timeFF] = saveTimeTable(edge_id, path_id, length, velocity, stagePoint);
[edge_DF, edge_FF] = divide(edge_id, stagePoint);

%% Regularization
alpha = 1;
fun = @(x) regularize(x, f, numOfEdges, alpha);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options2 = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options3 = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals', 20000,'MaxIter', 10000);

% v1 = sum(path_x(1:18*numOfEdges).^2)
% var1 = var(path_x(1:18*numOfEdges))
% [x, xval] = fmincon(fun,path_x,A,b,Aeq,beq,lb,ub,[],options3);
% %[x2, xval2] = fmincon(fun,x,A,b,Aeq,beq,lb,ub,[],options);
% v2 = sum(x(1:18*numOfEdges).^2)
% var2 = var(x(1:18*numOfEdges))

%% Information
fun = @(x) inform(x, f, numOfEdges, alpha);
v1 = sum(path_x(1:18*numOfEdges).^2)
var1 = var(path_x(1:18*numOfEdges))
[x, xval] = fmincon(fun,path_x,A,b,Aeq,beq,lb,ub,[],options3);
[x2, xval2] = fmincon(fun,x,A,b,Aeq,beq,lb,ub,[],options);
v2 = sum(x(1:18*numOfEdges).^2)
var2 = var(x(1:18*numOfEdges))