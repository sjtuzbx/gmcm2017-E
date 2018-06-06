function y = regularize(x, f, numOfEdges, alpha)
    t = 0;
    for i=1:18*numOfEdges
        t = t + x(i)^2;
    end
    y = f' * x + t;
end