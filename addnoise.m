function noise = addnoise(N,Mean,STD)

y=randn(1,N); 
y=y/std(y); 
y=y-mean(y); 
a=Mean; 
b=STD; 
y=a+b*y; 

noise = y;
end