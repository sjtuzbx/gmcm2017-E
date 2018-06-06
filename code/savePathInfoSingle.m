function [edge_id, path_id] = savePathInfoSingle(temp_x, in, out,edge, snode, numOfEdges, numOfF, F)
    path_id = cell(3,1);
    edge_id = cell(3,1);
    for i=1:3
        counter = 1;
        start = nan;
        for sp = 1:2
            s = snode(sp);
            for k=1 : out{s}(1)
                e = out{s}(1+k);
                if (temp_x(18*numOfEdges+2*numOfF+6+e, 1) > eps)
                    edge_id{i}(counter) = e;
                    start = s;
                    next = edge(e, 2);
                    temp_x(18*numOfEdges+2*numOfF+6+e, 1)= temp_x(18*numOfEdges+2*numOfF+6+e, 1)-1;
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

        path = [int2label(s),':',int2label(start)];
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
                if (temp_x(18*numOfEdges+2*numOfF+6+e, 1) > eps)
                    next = edge(e, 2);
                    edge_id{i}(counter-1) = e;
                    temp_x(18*numOfEdges+2*numOfF+6+e, 1) = temp_x(18*numOfEdges+2*numOfF+6+e, 1) - 1;
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
        disp(path)
        result_path{i} = path;
end