function [z,sum,fn1,gn1,flag]=conjugated_Newton_Method(x,x0,h1,h2,f,n,n1)
g=sym(zeros(n,1));
for i=1:n
    g(i,1)=diff(f,x{i});
end
G=sym(zeros(n,n));
for i=1:n
    for j=1:n
        G(i,j)=diff(g(i,1),x{j});
    end
end
G1=inv(G);
t1=subs(G1,x,x0);
g0=subs(g,x,x0);
p0=(-1)*t1*g0;
phi=@(a)subs(f,x,x0+a*p0);
a0=goldensection1(0,1,n1,phi);
x3=x0+a0*p0;
sum=1;
flag=0;
while(norm(subs(g,x,x3),inf)>h2)%&&sum<4
    g3=subs(g,x,x3);
    t2=subs(G1,x,x3);
    if(mod(sum,n)==0)
        p0=(-1)*t2*g3;
    else
        b=(g3'*t2*g3)/(g0'*t1*g0);
        p0=(-1)*subs(G1,x,x3)*subs(g,x,x3)+b*p0;
    end
        x0=x3;g0=g3;t1=t2;
        phii=@(a)subs(f,x,x0+a*p0);
        if(norm(g3)<=h1)
            a0=goldensection2(0,1,0.001,phii);
            flag=1;
        else
            a0=goldensection1(0,1,(sum+1)*n1,phii);
        end
        x3=x0+a0*p0;
        sum=sum+1;
end
fn1=subs(f,x,x3);
gn1=norm(subs(g,x,x3),inf);
z=x3;
end
        
        
    

