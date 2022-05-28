function [D] = fake_pdist(X)
    D = zeros(1,length(X));
    count = 1;
    for i = 1:length(X)
        for j = (i+1):length(X)
            D(count) = norm(X(i,:)-X(j,:));
            count = count + 1;
        end
    end
end