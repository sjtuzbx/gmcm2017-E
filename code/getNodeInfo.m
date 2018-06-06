function log = getNodeInfo(edge_id, path_id,stagePoint, numOfPoints)
    %% Extra Information for logging
    numOfVar = zeros(24, 1);
    znode = zeros(24, 2);
    for i=1:24
        numOfVar(i, 1) = 2*(size(edge_id{i},2)-stagePoint(i,1)+1)-2;
        znode(i, :) = [2*(stagePoint(i,2)-stagePoint(i,1)+1)-2, 2*(stagePoint(i,2)-stagePoint(i,1)+1)-1];
    end

    log = cell(numOfPoints, 1);
    cur = 0;
    for i=1:24
        for j=1:size(edge_id{i},2)-stagePoint(i,1)
           if (2*j == znode(i,1))
               %disp([num2str(t1),' ', num2str(t2), ' ', num2str(t2-t1)])
               id = path_id{i}(j+1+stagePoint(i,1));
               if isempty(log{id})
                 log{id}(1) = 1;
                 log{id}(2) = cur+znode(i,1);
                 log{id}(3) = cur+numOfVar(i,1)+1;
                 log{id}(4) = cur+numOfVar(i,1)+2;
                 log{id}(5) = cur+znode(i,2);
               else
                 log{id}(1) = log{id}(1) + 1;
                 log{id}(4*log{id}(1)-2) = cur+znode(i,1);
                 log{id}(4*log{id}(1)-1) = cur+numOfVar(i,1)+1;
                 log{id}(4*log{id}(1)) = cur+numOfVar(i,1)+2;
                 log{id}(4*log{id}(1)+1) = cur+znode(i,2);
               end
           end
        end
        cur = cur + numOfVar(i,1) + 2;
    end