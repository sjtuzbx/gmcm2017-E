function [res, ceq] = volConstraint(x,log)
    n = size(log,1);
    counter = 1;
    for i=1:n
        if (~isempty(log{i}))
            number = log{i}(1);
            if (number >= 2)
                a = zeros(2, 1);
                c = zeros(2, 1);
                for p=1:number-1
                    for q=p+1:number
                         a(1) = x(log{i}(1+4*p-3)); c(1) = x(log{i}(1+4*p-1));
                         a(2) = x(log{i}(1+4*q-3)); c(2) = x(log{i}(1+4*q-1));
                         
                         [~,ind] = max(a); [~, ind2] = min(a);
                         res(counter) = -c(ind)+ 10 + c(ind2);
                         counter = counter + 1;
                    end
                end
            end
                         
            if (number >=3)
                a = zeros(3, 1);
                b = zeros(3, 1);
                c = zeros(3, 1);
                d = zeros(3, 1);
                for p=1:number-2
                    for q=p+1:number-1
                        for r=q+1:number
                            a(1) = x(log{i}(1+4*p-3)); b(1) = x(log{i}(4*p-2+1));c(1) = x(log{i}(1+4*p-1));d(1) = x(log{i}(1+4*p));
                            a(2) = x(log{i}(1+4*q-3)); b(2) = x(log{i}(4*q-2+1));c(2) = x(log{i}(1+4*q-1));d(2) = x(log{i}(1+4*q));
                            a(3) = x(log{i}(1+4*r-3)); b(3) = x(log{i}(4*r-2+1));c(3)= x(log{i}(1+4*r-1));d(3) = x(log{i}(1+4*r));
                            [value,ind] = max(b);
                            indicator = logical(ones(3,1));
                            indicator(ind) = 0;
                            res(counter) = min(d(indicator)) - b(ind);
                            counter = counter + 1;

                            res(counter) = max(c(indicator))+ 10 - c(ind);
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
    end
    ceq = [];
end