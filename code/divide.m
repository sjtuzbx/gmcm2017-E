function [edge_DF, edge_FF] = divide(edge_id, stagePoint)
    edge_DF = cell(24, 1);
    edge_FF = cell(24, 1);
    for i=1:24
        eid = edge_id{i};
        for j=1:size(eid,2)
            if (j<=stagePoint(i,1))
                edge_DF{i}(j) = eid(j);
            else
                edge_FF{i}(j-stagePoint(i,1)) = eid(j);
            end
        end
    end
end