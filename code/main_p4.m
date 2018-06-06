%% Solution to p3
% 数据预处理
clear; clc; close all;
n = 2 + 6 + 60 + 62;

load Point
load edge

D = 1:2; Z = 3:8; F = 9:68; J = 69:130;

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
%%
load priority
% priority = zeros(numOfPoints, 1);
% for i=1:numOfPoints
%     c = a{i};
%     if (c(1) == 'J')
%         priority(j2int(str2num(c(2:size(c,2))))) = str2num(b{i});
%     end
% end
% save priority
    %priority(i,1) = 

%% 计算权值
w1 = 0.3; w2 = 0.2; w3 = 0.5;
weight = zeros(62, 1);

maxD = -1; minD = inf;
for i=J(1):J(62)  
  maxD = max(maxD, in{i}(1));
  minD = min(minD, in{i}(1));
end

for i=J(1):J(62)
    if ((i>=j2int(1) &&i<=j2int(11)) || (i>=j2int(12) &&i<=j2int(20)))
        type = 1;
    else
        type = 0;
    end
    
    degree = in{i}(1);
    degree = (degree - minD) / (maxD - minD);
    
    weight(i) = w1 * type + w2 * degree + w3 * (priority(i) - min(priority(69:130))) / max(priority);
end
    
%% 生成新的初始数据
% 从最大的10个里面选3个
[value, index] = sort(weight, 'descend');
comb = nchoosek(index(1:10), 3);
min_time = inf;
time_table_best = []; 
log = [];
counter = 1;
history = zeros(120, 4);
in_backup = in;
out_backup = out;

for i= 1:size(comb, 1)
    delNode = comb(i, :);
    %isValidEdge = ones(2*numOfEdges, 1);
    in = in_backup;
    out = out_backup;
    disp([num2str(i),'th iteration:'])
    for j=1:3
        for k=1:in{delNode(j)}(1)
            st = edge(in{delNode(j)}(1+k),1);
            ind = find(out{st}(2:size(out{st}, 2)) == in{delNode(j)}(1+k));
            out{st}(ind+1) = [];
            out{st}(1) = out{st}(1) - 1;
        end
        for k=1:out{delNode(j)}(1)
            st = edge(out{delNode(j)}(1+k),2);
            ind = find(in{st}(2:size(in{st}, 2)) == out{delNode(j)}(1+k));
            in{st}(ind+1) = [];
            in{st}(1) = in{st}(1) - 1;
        end
    end
    [f, intcon, A, b, Aeq, beq, lb, ub] = generatePathParameters(numOfEdges, numOfF, numOfPoints,edge, length, velocity, in, out, []);
    [path_x, path_val] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);
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
    fun = @(x) f3' * x;
    nonlcon = @(x) volConstraint(x, znodeInfo);
    options_default = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed');
    options = optimoptions(@fmincon,'Algorithm','sqp','Display','off', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
    options2 = optimoptions(@fmincon,'Algorithm','active-set','Display','off', 'FunValCheck', 'on');%,'MaxFunEvals', 300000,'MaxIter', 10000);
    options3 = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','MaxFunEvals', 30000,'MaxIter', 10000); 
    options4 = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','MaxFunEvals', 100000); 
   %% fine tuning, not always same.
    % try default first.
    %[x_non_default, fval_non_default] = fmincon(fun,ff_path,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon,options_default);

    %%
%         [x_non, fval_non] = fmincon(fun,ff_path,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options_default);
%         [x_non2, fval_non2] = fmincon(fun,x_non,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);
%         if ~(isempty(find((nonlcon(x_non) > 0.001))) && (sum(find(A3*x_non > b3)) == 0) && (sum(abs(Aeq3 * x_non - beq3) < 0.0001) == size(beq3, 1)))
%             %[x_non2, fval_non2] = fmincon(fun,x_non,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
%             [x_non3, fval_non3] = fmincon(fun,x_non2,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
%             [x_non4, fval_non4] = fmincon(fun,x_non3,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options2);
%             [x_non5, fval_non5] = fmincon(fun,x_non4,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options4);
%             x_non = x_non3;
%             log = [log,' ', num2str(i), ' ', num2str(j)];
%         end
    [x_non, fval_non] = fmincon(fun,ff_path,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
    [x_non2, fval_non2] = fmincon(fun,x_non,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);
    [x_non3, fval_non3] = fmincon(fun,x_non2,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
    [x_non4, fval_non4] = fmincon(fun,x_non3,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);
    %[x_non5, fval_non5] = fmincon(fun,x_non4,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options3);
    %[x_non6, fval_non6] = fmincon(fun,x_non5,A3,b3,Aeq3,beq3,lb3,ub3,nonlcon, options);

    %% F->F log information
    %load main_p1
    %time_table4 = saveTimeTableFF(edge_id, path_id, length, velocity, stagePoint, x_non6);
    time_tablefinal = saveTimeTableAll(edge_id, path_id, length, velocity, stagePoint, df_path, x_non);

    %%
    [total_time, trun, twait] = summary(x_non, edge_id, stagePoint, df_val, fval_non);
    disp(['Best running time ',num2str(path_val)])
    disp(['Our total time ', num2str(total_time)])
    disp(['Our running time ', num2str(sum(trun)+df_val)])
    disp(['Our waiting time ', num2str(sum(twait))])
    history(counter,:) = [path_val, total_time, sum(trun)+df_val, sum(twait)];
    counter = counter + 1;
    if (total_time < min_time)
        min_time = total_time;
        time_table_best = time_tablefinal;
        trun_best = trun+df_val;
        twait_best = twait;
    end
    save p4_solution
end

disp(time_table_best)
disp(['Our total time ', num2str(min_time)])
disp(['Our running time ', num2str(trun_best)])
disp(['Our waiting time ', num2str(twait_best)])

save p4_solution
