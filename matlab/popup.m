function z=popup(x,n,m)
z=zeros(n,m);
for i=1:n,
    for j=1:m,
        z(i,j)=x((i-1)*m + j);
    end
end