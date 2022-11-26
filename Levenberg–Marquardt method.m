function [z,sum,uu,ff,gl,fl]=LM(x,x0,h2,h3,f,n,J)%det,gg
g=sym(zeros(n,1));
for i=1:n
    g(i,1)=diff(f,x{i});
end
A=J'*J;u=max(diag(subs(A,x,x0)));G1=A+u*eye(n);v=2;t2=subs(G1,x,x0);
g0=subs(g,x,x0);
p01=(-1)*(t2\g0);
%phi=@(a)subs(f,x,x0+a*p0);
%a0=goldensection1(0,1,n1,phi);
%x3=x0+a0*p0;
x3=x0+p01;
sum=1;
%det=zeros(8,1);
%gg=zeros(8,1);
uu=zeros(8,1);% record the penalty parameter
uu(sum)=u;
ff=zeros(8,1);% record the objective function value
while(sum<5&&norm(subs(g,x,x3),inf)>h2) 
    if(norm(x3-x0)<=(h3*(norm(x0)+h3)))% jump out of the loop when x3 has small change
        break;
    else
        ff(sum)=subs(f,x,x3);
        if(sum==1)
            sigma=(subs(f,x,x0)-ff(sum))/((p01'*(u*p01-g0))/2);
        else
            sigma=(ff(sum-1)-ff(sum))/((p01'*(u*p01-g0))/2);
        end
        if(sigma>0)% Update U through sigma
            u=u*max(1/3,1-(2*sigma-1)^3);
            v=2;
        else
            u=u*v;v=2*v;
        end
        G1=A+u*eye(n);%¸üÐÂÄâº£Éª¾ØÕó
        t2=subs(G1,x,x3);
        g3=subs(g,x,x3);
        p01=(-1)*(t2\g3);
       %t3=norm(g3,inf);
       %det(sum)=abs(norm(g0,inf))-abs(t3);
       %gg(sum)=vpa(t3);
        x0=x3;
        g0=g3;
        x3=x0+p01;
        sum=sum+1;
        uu(sum)=u;
    end
end
g3=subs(g,x,x3);
% det(sum)=abs(norm(g0,inf))-abs(norm(g3,inf));
%gg(sum)=vpa(norm(g3,inf));
ff(sum)=subs(f,x,x3);
gl=vpa(norm(g3,inf));
fl=vpa(ff(sum));
z=x3;
end