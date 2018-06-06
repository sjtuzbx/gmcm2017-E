%% Solution to p3
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
[f, intcon, A, b, Aeq, beq, lb, ub] = generatePathParameters3_1(numOfEdges, numOfF, numOfPoints,edge, length, velocity, in, out, []);
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

%% J13, J14是解出来的隐蔽待机点
[edge_single, path_single] = savePathInfoSingle(path_x, in, out,edge, [j2int(13), j2int(14)], numOfEdges, numOfF, F);

%% D->F log information
time_table2 = saveTimeTableDF(edge_id, path_id, length, velocity, stagePoint, df_path);

%% F->F Linear Programming for initial solution
%与前面用的不同，需要修正
delCars = [21, 22, 23];
[f3,A3,b3,Aeq3,beq3,lb3,ub3,edge_id_new, path_id_new] = generateTimeParametersFF3(stagePoint, path_id, edge_FF, edge_id, length, velocity, delCars, edge_single, path_single);
znodeInfo = getNodeInfo(edge_id_new, path_id_new, stagePoint, numOfPoints);
[ff_path, ff_val] = linprog(f3,A3,b3,Aeq3,beq3,lb3,ub3); %ff_path作为方程组的初始解

%%
dt = 211.0657;
for i=21:23
    stagePoint(i, 1) = 0;
end
time_table3 = saveTimeTableLPFF3(edge_id_new, path_id_new, length, velocity, stagePoint, ff_path, dt);

%% F->F non-linear programming
fun = @(x) f3' * x;
nonlcon = @(x) volConstraint(x, znodeInfo);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options2 = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
options3 = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals', 30000,'MaxIter', 10000);

[x_non, fval_non] = fmincon(fun,ff_path,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
%%
time_tablefinal = saveTimeTableAll3(edge_id_new, path_id_new, length, velocity, stagePoint, df_path, x_non);

%%
[total_time, trun, twait] = summary(x_non, edge_id, stagePoint, df_val, fval_non)
disp(['Best running time ',num2str(path_val)])
disp(['Our total time ', num2str(total_time)])
disp(['Our running time ', num2str(sum(trun)+df_val)])
disp(['Our waiting time ', num2str(sum(twait))])

%%
c = 0;
for i=F(1):F(60)
     for j=1:out{i}(1)
         c = c + x(6*numOfEdges + out{i}(1+j)*3);
         if (x(6*numOfEdges + out{i}(1+j)*3) >0)
             disp(int2label(i))
         end
     end
end
c

%%
for i=1:numOfPoints
    numOut = 0;
    for j=1:out{i}(1)
        numOut = numOut + path_x(18*numOfEdges+2*numOfF+6+out{i}(1+j));
    end
    numIn = 0;
    for j=1:in{i}(1)
        numIn = numIn + path_x(18*numOfEdges+2*numOfF+6+in{i}(1+j));
    end
    if (numOut -numIn > 0)
        disp([int2label(i), ' ', num2str(numOut -numIn)])
    end
end

%% Decrease路程时间为

