function df_path = centeredTime(df_path_temp, stagePoint)
    df_path = df_path_temp;
    numOfVar = zeros(24, 1);
    for i=1:24
        numOfVar(i, 1) = 2*(stagePoint(i, 1)+1)-2;
    end
    
    minTime = inf;
    cur = 0;
    for i=1:24
        minTime = min(df_path_temp(cur+1), minTime);
        cur = cur + numOfVar(i,1);
    end
    
    cur=0;
    for i=1:24
        df_path(cur+1) = df_path_temp(cur+1) - minTime;
        cur = cur + numOfVar(i,1);
    end
end