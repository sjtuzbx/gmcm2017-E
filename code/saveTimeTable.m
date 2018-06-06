function [time_table, timeDF, timeFF] = saveTimeTable(edge_id, path_id, length, velocity, stagePoint)
    time_table = cell(24,1); 
    total_time = 0; timeDF = 0; timeFF = 0;
    for i=1:24
       eid = edge_id{i};
       if (i<=6)
           dtype = 1;
           letter = 'A';
       elseif (i<=12)
           dtype=2;
           letter = 'B';
       else
           dtype=3;
           letter = 'C';
       end
       t = 0; t1 = 0; t2 = 0;
       info = [letter,int2str(i),':',int2str(t),' ',int2label(path_id{i}(1))];
       for j=1:size(eid,2)
           if (j <= stagePoint(i,1))
               t1 = t1 + length(eid(j)) / velocity(eid(j),dtype) * 60; 
           else 
               t2 = t2 + length(eid(j)) / velocity(eid(j),dtype) * 60; 
           end
           t = t + length(eid(j)) / velocity(eid(j),dtype) * 60; 
           info = [info, ' ', num2str(t), '->',int2label(path_id{i}(j+1)) ];
       end
       total_time = total_time + t; timeDF = timeDF + t1; timeFF = timeFF + t2;
       disp(info)
       time_table{i} = info;
    end

    fid = fopen('time_table.txt','w');
    for i=1:24
        fprintf(fid,'%s\n',time_table{i});
    end
    fclose(fid);
    
end