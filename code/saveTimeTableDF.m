function time_table= saveTimeTableDF(edge_id, path_id, length, velocity, stagePoint, x_f)
    numOfVar = zeros(24, 1);
    for i=1:24
        numOfVar(i, 1) = 2*(stagePoint(i, 1)+1)-2;
    end
    cumVar = cumsum(numOfVar);

    time_table = cell(24,1);
    st = 1;
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
       t = x_f(st);
       info = [letter,int2str(i),':',int2label(path_id{i}(1)),' ',num2str(t)];
       for j=1:stagePoint(i,1)
           t1 = x_f(st+2*j-1);
           if (j~=stagePoint(i,1))
               t2 = x_f(st+2*j);
               info = [info, '->', num2str(t1), ' ',int2label(path_id{i}(j+1)),' ', num2str(t2)];
           else
               info = [info, '->', num2str(t1), ' ',int2label(path_id{i}(j+1))];
           end
       end
       st = cumVar(i)+1;
       disp(info)
       time_table{i} = info;
    end
end