function [f,A,b,Aeq,beq,lb,ub] = generateTimeParametersFF(stagePoint, path_id, edge_id, length, velocity)
    numOfVar = zeros(24, 1);
    znode = zeros(24, 2);
    for i=1:24
        numOfVar(i, 1) = 2*size(edge_id{i},2);
        znode(i, :) = [2*(stagePoint(i,2)-stagePoint(i,1)), 2*(stagePoint(i,2)-stagePoint(i,1))+1];
    end

    totalVar = sum(numOfVar) + 24 * 2;
    cumVar = cumsum(numOfVar+2);

    %%
    f = zeros(totalVar, 1);
    cur = 0;
    for i=1:24
        %f(cur+1) = -1; %从0时刻计算起
        f(cur+numOfVar(i,1)) = 1;
        f(cur+numOfVar(i,1)+1) = 1; %t2
        f(cur+znode(i,2)) = -1; % 
        cur = cur + numOfVar(i,1)+2;
    end

    %%
    cur = 0;
    counter = 1;
    A = zeros(1, totalVar);
    b = zeros(1, 1);
    counter2 = 1;

    Aeq = zeros(1, totalVar);
    beq = zeros(1, 1);

    for i=1:24
       eid = edge_id{i};
       if (i<=6)
           dtype = 1;
       elseif (i<=12)
           dtype=2;
       else
           dtype=3;
       end
       for j=1:numOfVar(i,1)-1
           if (mod(j,2) == 1)
                tempx = zeros(1, totalVar);
                tempx(1, cur+j) = 1;
                tempx(1, cur+j+1) = -1;
                Aeq(counter2, :) = zeros(1, totalVar);
                Aeq(counter2, :) = tempx;
                beq(counter2, 1) = -length(eid((j+1)/2)) / velocity(eid((j+1)/2),dtype) * 60;
                counter2 = counter2 + 1; 
%                 A(counter, :) = zeros(1, totalVar);
%                 A(counter, :) = tempx;
%                 b(counter, 1) = -length(eid((j+1)/2)) / velocity(eid((j+1)/2),dtype) * 60;
%                 counter = counter + 1; 

           else
                A(counter, :) = zeros(1, totalVar);
                tempx = zeros(1, totalVar);
                tempx(1, cur+j) = 1;
                tempx(1, cur+j+1) = -1;
                A(counter, :) = tempx;
                b(counter, 1) = 0;
                counter = counter + 1; 
           end
       end

       A(counter, :) = zeros(1, totalVar);
       tempx = zeros(1, totalVar);
       tempx(1, cur+znode(i,1)) = 1;
       tempx(1, cur+numOfVar(i,1)+1) = -1;
       A(counter, :) = tempx;
       b(counter, 1) = 0;
       counter = counter + 1;

       A(counter, :) = zeros(1, totalVar);
       tempx = zeros(1, totalVar);
       tempx(1, cur+numOfVar(i,1)+1) = 1;
       tempx(1, cur+numOfVar(i,1)+2) = -1;
       A(counter, :) = tempx;
       b(counter, 1) = 0;
       counter = counter + 1;

       A(counter, :) = zeros(1, totalVar);
       tempx = zeros(1, totalVar);
       tempx(1, cur+numOfVar(i,1)+2) = 1;
       tempx(1, cur+znode(i,2)) = -1;
       A(counter, :) = tempx;
       b(counter, 1) = -10;
       counter = counter + 1;

       cur = cur + numOfVar(i,1)+2;
    end

    %%

    for i=1:23
        for j=i+1:24
            Aeq(counter2, :) = zeros(1, totalVar);
            Aeq(counter2, cumVar(i)-2) = 1;
            Aeq(counter2, cumVar(j)-2) = -1;
            beq(counter2) = 0;
            counter2 = counter2 + 1;
        end
    end
    
    lb = zeros(totalVar, 1);
    ub = [];
    
    
end