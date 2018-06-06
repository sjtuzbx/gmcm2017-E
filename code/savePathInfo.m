function [result_path, path_id, edge_id, stagePoint] = savePathInfo(x, in, out, edge, numOfEquipment, F, Z, numOfEdges)
    result_path = cell(24,1);
    path_id = cell(24,1);
    edge_id = cell(24,1);
    stagePoint = zeros(24, 2);
    [result_path(1:6, 1), t1, t2, t3] = savePath(2, x, in, out, edge, numOfEquipment(1), F, Z, numOfEdges);
    path_id(1:6,1) = t1;
    edge_id(1:6,1) = t2;
    stagePoint(1:6,:) = t3;
    [result_path(7:12, 1), t1, t2, t3] = savePath(1, x, in, out, edge, numOfEquipment(2), F, Z, numOfEdges);
    path_id(7:12,1) = t1;
    edge_id(7:12,1) = t2;
    stagePoint(7:12,:) = t3;
    [result_path(13:24, 1), t1, t2, t3] = savePath(0, x, in, out, edge, numOfEquipment(3), F, Z, numOfEdges);
    path_id(13:24,1) = t1;
    edge_id(13:24,1) = t2;
    stagePoint(13:24,:) = t3;
    
    fid = fopen('path.txt','w');
    for i=1:24
        fprintf(fid,'%s\n',result_path{i});
    end
    fclose(fid);
    
end