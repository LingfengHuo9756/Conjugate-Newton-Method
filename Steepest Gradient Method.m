function [z,sum,flag,fd,gd,gg]=zuisu(x,x0,h1,h2,f,n,n1)
g=sym(zeros(n,1));
for i=1:n
    g(i,1)=diff(f,x{i});
end
g0=subs(g,x,x0);
sum=0;
flag=0;
gg=zeros(14,1);
while(sum<5&&norm(g0,inf)>h2) 
        p0=(-1)*g0;
        phii=@(a)subs(f,x,x0+a*p0);
        if(norm(g0,inf)<=h1)
            a0=goldensection2(0,1,0.01,phii);
            flag=1;
        else
            a0=goldensection1(0,1,(sum+1)*n1,phii);
        end
        x0=x0+a0*p0;
        g0=subs(g,x,x0);
        sum=sum+1;
        %gg(sum)=norm(g0,inf);
end
fd=subs(f,x,x0);
gd=norm(g0,inf);
z=x0;
end