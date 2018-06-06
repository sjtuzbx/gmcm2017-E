function time_tablefinal  = saveTimeTableAll(edge_id, path_id, length, velocity, stagePoint, x_f, x_f2)
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

    time_table = cell(24,1);
    st = 1; st2 = 1;
             
    for i=1:24
       dt = 0;
       if (i<=6)
           letter = 'A';
       elseif (i<=12)
           letter = 'B';
       else
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
               dt = t1;
               info = [info, '->', num2str(t1), ' ',int2label(path_id{i}(j+1))];
           end
       end
       st = cumVar(i)+1;

       t2 = x_f2(st2);
       info = [info,' ',num2str(t2+dt)];
       for j=1:size(edge_id{i},2)-stagePoint(i,1)
           t1 = x_f2(st2+2*j-1);
           if (j~=size(edge_id{i},2)-stagePoint(i,1))
               t2 = x_f2(st2+2*j);
               info = [info, '->', num2str(t1+dt), ' ',int2label(path_id{i}(j+1+stagePoint(i,1))),' ', num2str(t2+dt)];
           else
               info = [info, '->', num2str(t1+dt), ' ',int2label(path_id{i}(j+1+stagePoint(i,1)))];
           end
       end
       st2 = cumVar2(i)+1;

       disp(info)
       time_tablefinal{i} = info;
    end
end