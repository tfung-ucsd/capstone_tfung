% n = 10;
% test_matrix = zeros(n);
% stdev = 5;
% for i = 1:n
%     for j = (i+1):n
%         test_matrix(i,j) = randn*stdev;
%         test_matrix(j,i) = test_matrix(i,j);
%     end
% end
% test_matrix

labels={};
for i = 1:20
    labels = [labels num2str(i)];
end