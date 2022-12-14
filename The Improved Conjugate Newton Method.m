function [z,sum,ff,gn,fn]=quasi_conjugated_direction_Method(x,x0,h2,f,n,J)%,flag,det,gg
g=sym(zeros(n,1));
for i=1:n
    g(i,1)=diff(f,x{i});
end
A=J'*J;u=max(diag(subs(A,x,x0)));G1=A+u*eye(n);v=2;
g0=subs(g,x,x0);
k1=subs(G1,x,x0)\g0;
k2=g0'*k1;
p01=(-1)*k1;
p0=p01+0;
%phi=@(a)subs(f,x,x0+a*p0);
%a0=goldensection1(0,1,n1,phi);
%x3=x0+a0*p0;
x3=x0+p0;
x3_=x0+p01;
sum=1;%sum1=1;
%flag=0;
%det=zeros(8,1);
%gg=zeros(8,1);
ff=zeros(8,1);
while(sum<4&&norm(subs(g,x,x3),inf)>h2)
        ff(sum)=subs(f,x,x0);
        sigma=(ff(sum)-subs(f,x,x3_))/((p01'*(u*p01-g0))/2);
        if(sigma>0)% Update u through sigma
            u=u*max(1/3,1-(2*sigma-1)^3);
            v=2;
        else
            u=u*v;v=2*v;
        end
        k3=k2;
        G1=A+u*eye(n);% Update the Hessian Matrix
        g3=subs(g,x,x3);
        k1=subs(G1,x,x3)\g3;
        k2=g3'*k1;
        b=k2/k3;
        p01=(-1)*k1;
        if(mod(sum,n)==0)
            p0=p01;
        else
            p0=p01+b*p0;
        end
        x0=x3;
        g0=g3;
        x3=x0+p0;
        x3_=x0+p01;
        sum=sum+1;
end
g3=subs(g,x,x3);
ff(sum)=subs(f,x,x3);
gn=vpa(norm(g3,inf));
fn=vpa(ff(sum));
z=x3;
end