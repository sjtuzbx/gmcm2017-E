clear; clc
load p4_solution

[value, ind] = sort(history(:,2), 'descend');
comb(ind(1:3), :)
dt=7893.2
for i=1:3
    disp([int2label(comb(ind(i),1)),' ' ,int2label(comb(ind(i),2)),' ' ,int2label(comb(ind(i),3)), ' ', num2str(value(i)-7893)]);
end