function [tzong, trun, twait] = summary3(x_non6, edge_id, stagePoint, df_val, fval_non6)
    numOfVar = zeros(24, 1);
    numOfVar2 = zeros(24, 1);
    znode = zeros(24, 2);
    for i=1:24
        numOfVar(i, 1) = 2*(stagePoint(i, 1)+1)-2;
        numOfVar2(i, 1) = 2*(size(edge_id{i},2)-stagePoint(i,1)+1)-2;
        znode(i, :) = [2*(stagePoint(i,2)-stagePoint(i,1)+1)-2, 2*(stagePoint(i,2)-stagePoint(i,1)+1)-1];
    end
    
    cumVar = cumsum(numOfVar);
    cumVar2 = cumsum(numOfVar(1:20)+2);
    cumVar2(21) = cumVar2(20) + numOfVar(21);
    cumVar2(22) = cumVar2(21) + numOfVar(22);
    cumVar2(23) = cumVar2(22) + numOfVar(23);
    cumVar2(24) = cumVar2(23) + numOfVar(24)+2;
    
    tzong = df_val + fval_non6;
    x_f2 = x_non6;
    s0 = 0;
    trun = 0;
    twait = 0;
    for i=1:24
       eid = edge_id{i};
       trun(i) = 0;
       twait(i) = x_f2(s0+1);
       for j=1:numOfVar2(i,1)-1
           t1 = x_f2(s0+j);
           t2 = x_f2(s0+j+1);
           if (mod(j,2) == 1)
               trun(i) = trun(i) + t2 - t1;
           else
               twait(i) = twait(i) + t2 - t1;
           end
       end
       twait(i) = twait(i) - (x_f2(s0+znode(i,2)) - x_f2(s0+numOfVar2(i,1)+1));
       s0 = cumVar2(i);
    end 
end

%%
