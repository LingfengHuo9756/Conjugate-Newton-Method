f1=zeros(6,1);g1=zeros(6,1);t1=zeros(6,1);sum2=zeros(6,1);z2=zeros(6,6);
%% Conjugate Gradient Method
t0=clock;
[z,sum1,flag,fg,gg]=gonge1(x,x0,h1,h2,f,n,n1);
tg=etime(clock,t0); 
f1(1)=fg;g1(1)=gg;t1(1)=tg;sum2(1)=sum1;z2(:,1)=z;
%% Method with weighted average step
% t0=clock;
% [z2,sum,flag,flag1,fj,gj]=jiaquan(x,x0,h1,h2,f,n,n1);
% tj=etime(clock,t0);
% f1(2)=fj;g1(2)=gj;t1(2)=tj;sum2(2)=sum;
%% The Improved Conjugated Newton Method
t0=clock;
[z,sum,ff,gn,fn]=newgonge2(x,x0,h2,f,n,J);
tn=etime(clock,t0); 
f1(2)=fn;g1(2)=gn;t1(2)=tn;sum2(2)=sum;z2(:,2)=z;
%% Levenberg¨CMarquardt Method
t0=clock;
[z,sum,uu,ff,gl,fl]=LM(x,x0,h2,h3,f,n,J);
tl=etime(clock,t0); 
f1(3)=fl;g1(3)=gl;t1(3)=tl;sum2(3)=sum;z2(:,3)=z;
%% Conjugated Newton Method
% t0=clock;
% [z,sum,fn1,gn1,flag]=newgonge(x,x0,h1,h2,f,n,n1);
% tn1=etime(clock,t0);
% f1(3)=fn1;g1(3)=gn1;t1(3)=tn1;sum2(3)=sum;
%% Steepest Gradient Method
t0=clock;
[z,sum,flag,fd,gd,gg]=zuisu(x,x0,h1,h2,f,n,n1);
td=etime(clock,t0);
f1(4)=fd;g1(4)=gd;t1(4)=td;sum2(4)=sum;z2(:,4)=z;

