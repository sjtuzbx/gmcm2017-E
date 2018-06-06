function [result_path, path_id, edge_id,  stagePoint] = savePath(type, x, in, out, edge, numOfEquipment, F, Z, numOfEdges)
    % type: 2 = A, 1 = B, 0 = C
    if (type == 2)
        letter = 'A';
    elseif (type == 1)
        letter = 'B';
    else
        letter = 'C';
    end
    
    stagePoint = zeros(numOfEquipment, 2);
    temp_x = x;
    eps = 0.0001;
    result_path = cell(numOfEquipment,1);
    path_id = cell(numOfEquipment,1);
    edge_id = cell(numOfEquipment,1);
    
    for i=1:numOfEquipment
        counter = 1;
        start = nan;
        for s = 1:2 
            for k=1 : out{s}(1)
                e = out{s}(1+k);
                if (temp_x(3*e-type, 1) > eps)
                    edge_id{i}(counter) = e;
                    start = s;
                    next = edge(e, 2);
                    temp_x(3*e-type, 1)= temp_x(3*e-type, 1)-1;
                    break;
                end
            end
            if (~isnan(start))
                break;
            end
        end
        if (isnan(start))
            break;
        end
        
        path = [letter, int2str(i),':',int2label(start)];
        path_id{i}(counter) = start;
        counter = counter + 1;
        start = next;
        while (next ~= -1)
            path = [path,'->',int2label(next)];
            path_id{i}(counter) = next;
            counter = counter + 1;
            next = -1;
            for k=1:out{start}(1)
                e = out{start}(1+k);
                if (temp_x(3*e-type, 1) > eps)
                    next = edge(e, 2);
                    edge_id{i}(counter-1) = e;
                    temp_x(3*e-type, 1) = temp_x(3*e-type, 1) - 1;
                    break;
                end
            end
            if (sum(find(start == F) ~= 0) && next == -1)
                %path = [path,'->',int2label(next)];
                stagePoint(i, 1) = counter-2;
                break;
            end
            start = next;
        end
        
        next = start;
        % F->Z
        while (next ~= -1)
            next = -1;
            for k=1:out{start}(1)
                e = out{start}(1+k);
                if (temp_x(6*numOfEdges+3*e-type, 1) > eps)
                    next = edge(e, 2);
                    edge_id{i}(counter-1) = e;
                    temp_x(6*numOfEdges+3*e-type, 1) = temp_x(6*numOfEdges+3*e-type, 1) - 1;
                    break;
                end
            end
            if (sum(find(start == Z) ~= 0) && next == -1)
                %path = [path,'->',int2label(next)];
                stagePoint(i, 2) = counter-2;
                break;
            end
            start = next;
            path = [path,'->',int2label(next)];
            path_id{i}(counter) = next;
            counter = counter + 1;
        end
        
        next = start;
        % Z->F
        while (next ~= -1)
            next = -1;
            for k=1:out{start}(1)
                e = out{start}(1+k);
                if (temp_x(12*numOfEdges+3*e-type, 1) > eps)
                    next = edge(e, 2);
                    edge_id{i}(counter-1) = e;
                    temp_x(12*numOfEdges+3*e-type, 1) = temp_x(12*numOfEdges+3*e-type, 1) - 1;
                    break;
                end
            end
            if (sum(find(start == F) ~= 0) && next == -1)
                path = [path,'->',int2label(next)];
                %path_id{i}(counter) = next;
                %counter = counter + 1;
                break;
            end
            start = next;
            path = [path,'->',int2label(next)];
            path_id{i}(counter) = next;
            counter = counter + 1;
        end

        %disp(path)
        result_path{i} = path;
    end
end