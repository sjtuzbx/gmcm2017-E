function [f, intcon, A, b, Aeq, beq, lb, up] = generatePathParameters3_1(numOfEdges, numOfF, numOfPoints,edge, length, velocity, in, out, node)
    %% 这个文件不同于generatePathParameters3
    % 在generatePathParameters3中，我们把增加的3台C装置放入原图考虑，这样子会强制它通过Z节点
    % 在本约束中，我们用另外一个图来表示
    idC = [j2int(4), j2int(6), j2int(8), j2int(13), j2int(14), j2int(15)];
    D = 1:2;
    Z = 3:8;
    F = 9:68;
    J = 69:130;
    
    Z = cat(2, Z, node);
    if (~isempty(node))
        J(find(J==node(1))) = [];
        J(find(J==node(2))) = [];
    end
     
    f = zeros(18*numOfEdges+2*numOfF+size(idC,2)+2*numOfEdges, 1); %Add 6+|E| more variables
    for i=1:6*numOfEdges
        index = i;
        while (index > 2*numOfEdges)
            index = index - 2 * numOfEdges;
        end
        f(3*i-2,1) = length(index) / velocity(index, 1) * 60;
        f(3*i-1,1) = length(index) / velocity(index, 2) * 60;
        f(3*i,1) = length(index) / velocity(index, 3) * 60;
    end
    
    for i=1:2*numOfEdges
        f(18*numOfEdges+2*numOfF+size(idC,2)+i) = length(i) / velocity(i, 3)*60; % 因为是C装置
    end

    intcon = 1:size(f,1);
    A = zeros(numOfF, size(f,1));
    for i=1:numOfF
        A(i, i+18*numOfEdges) = 1;
        A(i, i+18*numOfEdges+numOfF) = 1;
    end
    b = ones(numOfF, 1);
    
    % D -> F
    %Aeq = 0;
    counter = 1;
    counter2 = numOfF;
    % \Sigma n_i = 3
    Aeq = zeros(numOfF, size(f,1));
    for i=1:size(idC,2)
        Aeq(counter, 18*numOfEdges+2*numOfF+i) = 1;
    end
    beq(counter) = 3;
    counter = counter + 1;
    
    t = [6, 3, 3];
    for i=1:numOfPoints
        if (sum(find(i == D)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:out{i}(1)
                    Aeq(counter, out{i}(1+j)*3-k+1) = 1;
                end
                for j=1:in{i}(1)
                    Aeq(counter, in{i}(1+j)*3-k+1) = -1;
                end
                beq(counter, 1) = t(k);
                counter = counter + 1;
            end
        elseif (sum(find(i == J)) ~= 0 || sum(find(i == Z)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:out{i}(1)
                    Aeq(counter, out{i}(1+j)*3-k+1) = 1;
                end
                for j=1:in{i}(1)
                    Aeq(counter, in{i}(1+j)*3-k+1) = -1;
                end
                beq(counter, 1) = 0;
                counter = counter + 1;
            end
        elseif (sum(find(i == F)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:in{i}(1)
                    Aeq(counter, in{i}(1+j)*3-k+1) = 1;
                end
            end
            Aeq(counter, i-8+18*numOfEdges) = -1;
            beq(counter, 1) = 0;
            counter = counter + 1;
        end
    end

    % F -> Z
    for i=1:numOfPoints
        if (sum(find(i == F)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                if (k~=1) %B和A型
                    for j=1:out{i}(1)
                        Aeq(counter, 6*numOfEdges+out{i}(1+j)*3-k+1) = 1;
                    end
                    for j=1:in{i}(1)
                        Aeq(counter, in{i}(1+j)*3-k+1) = -1;
                        Aeq(counter, 6*numOfEdges+in{i}(1+j)*3-k+1) = -1; % 这里要注意Z->F的约束！！
                    end
                    beq(counter, 1) = 0;
                    counter = counter + 1;
                else %C型
                    for j=1:in{i}(1)
                        Aeq(counter, 6*numOfEdges+in{i}(1+j)*3-k+1) = 1; 
                    end
                    beq(counter, 1) = 0;
                    counter = counter + 1;
                    
                    for j=1:out{i}(1)
                        Aeq(counter, out{i}(1+j)*3-k+1) = 1;
                    end
                    beq(counter, 1) = 0;
                    counter = counter + 1;
                    
                    % chu - ru <= 0
                    A(counter2, :) = zeros(1, size(f,1));
                    for j=1:out{i}(1)
                        A(counter2, 6*numOfEdges+out{i}(1+j)*3-k+1) = 1;
                    end
                    for j=1:in{i}(1)
                        A(counter2, in{i}(1+j)*3-k+1) = -1;
                    end
                    b(counter2, 1) = 0;
                    counter2 = counter2 + 1;                   
                end
              end  
        elseif (sum(find(i == J)) ~= 0 || sum(find(i == D)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:out{i}(1)
                    Aeq(counter, 6*numOfEdges + out{i}(1+j)*3-k+1) = 1;
                end
                for j=1:in{i}(1)
                    Aeq(counter, 6*numOfEdges + in{i}(1+j)*3-k+1) = -1;
                end
                beq(counter, 1) = 0;
                counter = counter + 1;
            end
        elseif (sum(find(i == Z)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:out{i}(1)
                    Aeq(counter, 12*numOfEdges + out{i}(1+j)*3-k+1) = 1;
                    Aeq(counter, 6*numOfEdges + out{i}(1+j)*3-k+1) = 1; % new added, might influence p1
                end
                for j=1:in{i}(1)
                    Aeq(counter, 6*numOfEdges + in{i}(1+j)*3-k+1) = -1;
                    Aeq(counter, 12*numOfEdges + in{i}(1+j)*3-k+1) = -1;
                end
                beq(counter, 1) = 0;
                counter = counter + 1;
            end
        end
    end
    
    % 离开第二阶段F的C车为9
    Aeq(counter, :) = zeros(1, size(f, 1));
    for i=F(1):F(60)
         for j=1:out{i}(1)
             Aeq(counter, 6*numOfEdges + out{i}(1+j)*3) = 1;
         end
    end
    beq(counter, 1) = 9;
    counter = counter + 1;

    for k=2:3
        Aeq(counter, :) = zeros(1, size(f, 1));
        for i=Z(1):Z(6)
             for j=1:in{i}(1)
                 Aeq(counter, 6*numOfEdges + in{i}(1+j)*3-k+1) = 1;
             end
        end
        beq(counter, 1) = t(k) * 2;
        counter = counter + 1;
    end
    
    %进入Z的C车为9
    Aeq(counter, :) = zeros(1, size(f,1));
    for i=Z(1):Z(6)
         for j=1:in{i}(1)
             Aeq(counter, 6*numOfEdges + in{i}(1+j)*3) = 1;
         end
    end
    beq(counter, 1) = t(1) * 2-3 ;
    counter = counter + 1;

    % Z -> F
    for i=1:numOfPoints
        if (sum(find(i == J)) ~= 0 || sum(find(i == D)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:out{i}(1)
                    Aeq(counter, 12*numOfEdges+out{i}(1+j)*3-k+1) = 1;
                end
                for j=1:in{i}(1)
                    Aeq(counter, 12*numOfEdges+in{i}(1+j)*3-k+1) = -1;
                end
                beq(counter, 1) = 0;
                counter = counter + 1;
            end
        elseif (sum(find(i == F)) ~= 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for k = 1:3
                for j=1:in{i}(1)
                    Aeq(counter, 12*numOfEdges+in{i}(1+j)*3-k+1) = 1;
                end
            end
            for j=1:in{i}(1)
                Aeq(counter, 18*numOfEdges+2*numOfF+size(idC,2)+in{i}(1+j)) = 1;
            end
            Aeq(counter, i-8+18*numOfEdges+numOfF) = -1;
            beq(counter, 1) = 0;
            counter = counter + 1;
        end
    end
    
    %处理额外的图
    for i=1:numOfPoints
        if (sum(i==idC) == 1)
            Aeq(counter, :) = zeros(1, size(f, 1));    
            for j=1:out{i}(1)
                Aeq(counter, 18*numOfEdges+2*numOfF+size(idC,2)+out{i}(1+j)) = 1;
            end
            for j=1:in{i}(1)
                Aeq(counter, 18*numOfEdges+2*numOfF+size(idC,2)+in{i}(1+j)) = -1;
            end
            [~, ind] = find(i ==idC);
            Aeq(counter, 18*numOfEdges+2*numOfF+ind) = -1;
            beq(counter, 1) = 0;
            counter = counter + 1;        
        elseif (sum(find(i == F)) == 0)
            Aeq(counter, :) = zeros(1, size(f, 1));
            for j=1:out{i}(1)
                Aeq(counter, 18*numOfEdges+2*numOfF+size(idC,2)+out{i}(1+j)) = 1;
            end
            for j=1:in{i}(1)
                Aeq(counter, 18*numOfEdges+2*numOfF+size(idC,2)+in{i}(1+j)) = -1;
            end
            beq(counter, 1) = 0;
            counter = counter + 1;           
       end
    end

    Aeq(counter, :) = zeros(1, size(f, 1));
    for i=1:numOfF
         Aeq(counter, i+18*numOfEdges) = 1;
    end
    beq(counter, 1) = 24;
    counter = counter + 1;
    

    Aeq(counter, :) = zeros(1, size(f, 1));
    for i=1:numOfF
         Aeq(counter, i+numOfF+18*numOfEdges) = 1;
    end
    beq(counter, 1) = 24;
    counter = counter + 1;

    lb = zeros(size(f));
    up = zeros(size(f));
    for i=1:2*numOfEdges
        up(3*i-2, 1) = 6;
        up(3*i-1, 1) = 6;
        up(3*i, 1) = 12;
        up(3*i-2+6*numOfEdges, 1) = 6;
        up(3*i-1+6*numOfEdges, 1) = 6;
        up(3*i+6*numOfEdges, 1) = 12;
        up(3*i-2+12*numOfEdges, 1) = 6;
        up(3*i-1+12*numOfEdges, 1) = 6;
        up(3*i+12*numOfEdges, 1) = 12;
    end
    for i=1:numOfF
        up(18*numOfEdges+i, 1) = 1;
        up(18*numOfEdges+numOfF+i, 1) = 1;
    end
    for i=1:size(idC,2)
        up(18*numOfEdges+numOfF*2+i, 1) = 2;
    end
    for i=1:2*numOfEdges
         up(18*numOfEdges+numOfF*2+size(idC,2)+i, 1) = 3;
    end
end


