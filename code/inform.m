function y = inform(x, f, numOfEdges, alpha)
    t = 0;
    for i=1:18*numOfEdges
        t = t + x(i)*log(x(i));
    end
    y = f' * x + t;
end