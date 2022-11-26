function hahaha=goldensection1(x0,det,n,f)
[a,b]=jintui1(x0,det,f);
x1=a+0.382*(b-a);
x2=a+0.618*(b-a);
sum=1;
while(sum<=n)
    f1=f(x1);
    f2=f(x2);
    if(f1>f2)
        a=x1;x1=x2;x2=a+0.618*(b-a);
    elseif(f1==f2)
        a=x1;b=x2;
    else
        b=x2;x2=x1;x1=a+0.382*(b-a);
    end
    sum=sum+1;
end
hahaha=(a+b)/2;
end
        
        