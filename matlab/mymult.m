function c=mymult(a,b);

[na,ma]=size(a);
[nb,mb]=size(b);
bigb=zeros(na*nb,ma*mb);  % just to declare it to help with memory 
for io=1:nb,
    for jo=1:mb,      % it would be nice to do this without for loops, but i don't know how
        for ii=1:na,
            for ji=1:ma,   % for i-inner, i-outer, j-inner, j-outer
                i=(io-1)*na + ii;
                j=(jo-1)*ma + ji;
                bigb(i,j)=b(io,jo);
            end
        end
    end
end
biga = repmat(a,nb,mb);
%size(biga)
%size(bigb)
c = biga .* bigb;
