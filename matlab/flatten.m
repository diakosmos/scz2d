function z=flatten(x);
%
%  flattens arrays. probably already exists, but not to my knowledge.
%
[n,m]=size(x);
z=ones(1,n*m)';
for i=1:n,
    for j=1:m,
        z((i-1)*m + j) = x(i,j);
    end
end
