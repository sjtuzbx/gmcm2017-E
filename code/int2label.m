function out = int2label(x)
    if x <=2
        out = ['D' int2str(x)];
    elseif x <= 8
        out = ['Z' int2str(x-2)];
    elseif x <= 68
        out = ['F' int2str(x-8)];
    else
        out = ['J' int2str(x-68)];
    end
end