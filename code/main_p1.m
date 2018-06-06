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

%% D->F Linear Programming
[f2, A2, b2, Aeq2, beq2, lb2, ub2] = generateTimeParameters(stagePoint, edge_DF, length, velocity);
[df_path_temp, df_val] = linprog(f2,A2,b2,Aeq2,beq2,lb2,ub2);

%% zero-centered
df_path = centeredTime(df_path_temp, stagePoint);


%% D->F log information
time_table2 = saveTimeTableDF(edge_id, path_id, length, velocity, stagePoint, df_path);

%% F->F Linear Programming for initial solution
[f3,A3,b3,Aeq3,beq3,lb3,ub3] = generateTimeParametersFF(stagePoint, path_id, edge_FF, length, velocity);
znodeInfo = getNodeInfo(edge_id, path_id, stagePoint, numOfPoints);
[ff_path, ff_val] = linprog(f3,A3,b3,Aeq3,beq3,lb3,ub3); %ff_path作为方程组的初始解

%% F->F infor
time_table3 = saveTimeTableLPFF(edge_id, path_id, length, velocity, stagePoint, ff_path);

%% F->F non-linear programming
fun = @(x) f3' * x;
nonlcon = @(x) volConstraint(x, znodeInfo);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options2 = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options3 = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals', 30000,'MaxIter', 10000);

[x_non, fval_non] = fmincon(fun,ff_path,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
[x_non2, fval_non2] = fmincon(fun,x_non,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);
[x_non3, fval_non3] = fmincon(fun,x_non2,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
[x_non4, fval_non4] = fmincon(fun,x_non3,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);
[x_non5, fval_non5] = fmincon(fun,x_non4,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
[x_non6, fval_non6] = fmincon(fun,x_non5,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);

%% check
if (isempty(find((nonlcon(x_non6) > 0.001))) && (sum(find(A3*x_non6 > b3)) == 0) && (sum(abs(Aeq3 * x_non6 - beq3) < 0.0001) == size(beq3, 1)))
    disp('Not a feasible result')
end

%% F->F log information
%load main_p1
%time_table4 = saveTimeTableFF(edge_id, path_id, length, velocity, stagePoint, x_non6);
time_tablefinal = saveTimeTableAll(edge_id, path_id, length, velocity, stagePoint, df_path, x_non6);

%%
[total_time, trun, twait] = summary(x_non6, edge_id, stagePoint, df_val, fval_non6)
disp(['Best running time ',num2str(path_val)])
disp(['Our total time ', num2str(total_time)])
disp(['Our running time ', num2str(sum(trun)+df_val)])
disp(['Our waiting time ', num2str(sum(twait))])
% for i=1:8
%     if (~isempty(znodeInfo{i}))
%         number = znodeInfo{i}(1);
%         for k=1:number
%             disp([num2str(x_non6(znodeInfo{i}(1+4*k-3))),' ',num2str(x_non6(znodeInfo{i}(1+4*k-2))),' ' ,num2str(x_non6(znodeInfo{i}(1+4*k-1))),' ',num2str(x_non6(znodeInfo{i}(1+4*k)))]);
%         end
%     end
% end

