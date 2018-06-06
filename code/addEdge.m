function out = addEdge(graph, a, e)
    if isempty(graph{a})
        graph{a}(1) = 1;
        graph{a}(2) = e;
    else
        graph{a}(graph{a}(1)+2) = e;
        graph{a}(1) = graph{a}(1) + 1;
    end
    out = graph;
end