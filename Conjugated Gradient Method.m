function [z,sum1,flag,fg,gg]=gonge1(x,x0,h1,h2,f,n,n1) %,fg,gg,gg1
g=sym(zeros(n,1));
for i=1:n
    g(i,1)=diff(f,x{i});
end
g0=subs(g,x,x0);
p0=(-1)*g0;
phi=@(a)subs(f,x,x0+a*p0);
a0=goldensection1(0,1,n1,phi);
x3=x0+a0*p0;
sum1=1;
sum=1;
flag=0;
%gg1=zeros(10,1);
while(sum1<4&&norm(subs(g,x,x3),inf)>h2)%norm(subs(g,x,x3))>h2;sum1<4;
    g3=subs(g,x,x3);
    %gg1(sum1)=norm(g3,inf);
    if(mod(sum1,n)==0)
        p0=(-1)*g3;     
    else
        b=(g3'*g3)/(g0'*g0);
        p0=(-1)*g3+b*p0;
    end
        x0=x3;
        g0=g3;
        phii=@(a)subs(f,x,x0+a*p0);
        if(norm(g3)<=h1)
            a0=goldensection2(0,1,0.01,phii);
            flag=1;
        else
            if(mod(sum1,2)==0)
                sum=0;
                a0=goldensection1(0,1,(sum+1)*n1,phii);
            else
                a0=goldensection1(0,1,(sum+1)*n1,phii);
            end
        end
        x3=x0+a0*p0;
        sum1=sum1+1;sum=sum+1;
end
fg=subs(f,x,x3);
gg=norm(subs(g,x,x3),inf);
%gg1(sum1)=gg;
z=x3;
end
        
        
    

