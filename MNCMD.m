function [IFmset,IA,smset]=MNCMD(s,fs,eIF,alpha,beta,var,tol)

[mm,nn]=size(s);
if mm>nn  
    s=s';
end
clear mm nn
[M,N]=size(s);
[kk,nn]=size(eIF);
if kk>nn  
    eIF=eIF';
end
clear kk nn
[K,N]=size(eIF);

t=(0:N-1)/fs; 
e=ones(N,1);
e2=-2*e; 
e2(1)=-1; e2(end)=-1;
oper=spdiags([e e2 e],-1:1,N,N); 
opedoub = oper'*oper;
sini=zeros(K,N);cosi=zeros(K,N);
xim=zeros(K,N,M);yim=zeros(K,N,M);
iternum=500; 
IFsetiter=zeros(K,N,iternum+1); IFsetiter(:,:,1)=eIF; 
ssetiter=zeros(K,N,M,iternum+1); 
lamuda=zeros(M,N); 
u=zeros(M,N); 
sum_x=zeros(M,N);
sum_y=zeros(M,N);
interx=zeros(M,N);
amp2=zeros(M,N); 
ybar=zeros(M,N);
xbar=zeros(M,N);
deltaIF=zeros(M,N);
tempIF=zeros(M,N);
for i=1:K
    sini(i,:)=sin(2*pi*(cumtrapz(t,eIF(i,:))));
    cosi(i,:)=cos(2*pi*(cumtrapz(t,eIF(i,:))));
    Bi=spdiags(sini(i,:)',0,N,N); Bdoubi=spdiags((sini(i,:).^2)',0,N,N); 
    Ai=spdiags(cosi(i,:)',0,N,N); Adoubi=spdiags((cosi(i,:).^2)',0,N,N); 
    for m=1:M
        xim(i,:,m)=(2/alpha*opedoub+Adoubi)\(Ai'*s(m,:)');
        yim(i,:,m)=(2/alpha*opedoub+Bdoubi)\(Bi'*s(m,:)');
        ssetiter(i,:,m,1)=xim(i,:,m).*cosi(i,:)+yim(i,:,m).*sini(i,:);
    end
end


iter=1; 
sDif=tol+1;
for m=1:M
    sum_x(m,:)=sum(xim(:,:,m).*cosi,1); 
    sum_y(m,:)=sum(yim(:,:,m).*sini,1); 
end

while (sDif>tol && iter<=iternum)
    
    betathr=10^(iter/36-10); 
    if betathr>beta
        betathr=beta;
    end
    
    for m=1:M
        u(m,:)=projec(s(m,:)-sum_x(m,:)-sum_y(m,:)-lamuda(m,:)/alpha,var);
    end
    
    for i=1:K
        
        Bi=spdiags(sini(i,:)',0,N,N); Bdoubi=spdiags((sini(i,:).^2)',0,N,N); 
        Ai=spdiags(cosi(i,:)',0,N,N); Adoubi=spdiags((cosi(i,:).^2)',0,N,N); 
        for m=1:M
            
            sum_x(m,:)=sum_x(m,:)-xim(i,:,m).*cosi(i,:); 
            xim(i,:,m)=(2/alpha*opedoub+Adoubi)\(Ai'*(s(m,:)-sum_x(m,:)-sum_y(m,:)-u(m,:)-lamuda(m,:)/alpha)');
            interx(m,:)=xim(i,:,m).*cosi(i,:); 
            sum_x(m,:)=sum_x(m,:)+interx(m,:);
            
            sum_y(m,:)=sum_y(m,:)-yim(i,:,m).*sini(i,:); 
            yim(i,:,m)=(2/alpha*opedoub+Bdoubi)\(Bi'*(s(m,:)-sum_x(m,:)-sum_y(m,:)-u(m,:)-lamuda(m,:)/alpha)');
            
            amp2(m,:)=xim(i,:,m).^2+yim(i,:,m).^2;
            
            ybar(m,:)=Differ(yim(i,:,m),1/fs); xbar(m,:)=Differ(xim(i,:,m),1/fs);
            deltaIF(m,:)=(xim(i,:,m).*ybar(m,:)-yim(i,:,m).*xbar(m,:))./(amp2(m,:))/2/pi; 
            temp=(2/betathr*opedoub+speye(N))\deltaIF(m,:)';
            tempIF(m,:)=eIF(i,:)-0.5*temp'; 
        end
        
        eIF(i,:)=sum(tempIF.*amp2,1)./sum(amp2,1);
        
       
        sini(i,:)=sin(2*pi*(cumtrapz(t,eIF(i,:))));
        cosi(i,:)=cos(2*pi*(cumtrapz(t,eIF(i,:))));
        
       
        for m=1:M
            sum_x(m,:)=sum_x(m,:)-interx(m,:)+xim(i,:,m).*cosi(i,:);
            sum_y(m,:)=sum_y(m,:)+yim(i,:,m).*sini(i,:);
            ssetiter(i,:,m,iter+1)=xim(i,:,m).*cosi(i,:)+yim(i,:,m).*sini(i,:);
        end
    end
    IFsetiter(:,:,iter+1)=eIF;
    
    
    for m=1:M
        lamuda(m,:)=lamuda(m,:)+alpha*(u(m,:)+sum_x(m,:)+sum_y(m,:)-s(m,:));
    end
    
   
    for m=1:M
        if norm(u(m,:)+sum_x(m,:)+sum_y(m,:)-s(m,:))>norm(s(m,:))
            lamuda(m,:)=zeros(1,length(t));
            for i=1:K
                Bi=spdiags(sini(i,:)',0,N,N); Bdoubi=spdiags((sini(i,:).^2)',0,N,N); % Bdoubm=Bm'*Bm
                Ai=spdiags(cosi(i,:)',0,N,N); Adoubi=spdiags((cosi(i,:).^2)',0,N,N); % Adoubm=Am'*Am
                xim(i,:,m)=(2/alpha*opedoub+Adoubi)\(Ai'*s(m,:)');
                yim(i,:,m)=(2/alpha*opedoub+Bdoubi)\(Bi'*s(m,:)');
                ssetiter(i,:,m,iter+1)=xim(i,:,m).*cosi(i,:)+yim(i,:,m).*sini(i,:);  
            end
            sum_x(m,:)=sum(xim(:,:,m).*cosi,1); 
            sum_y(m,:)=sum(yim(:,:,m).*sini,1); 
        end
    end
    
   
    sDif=0;
    for i=1:K
        for m=1:M
            sDif=sDif+(norm(ssetiter(i,:,m,iter+1)-ssetiter(i,:,m,iter))/norm(ssetiter(i,:,iter))).^2;
        end
    end
    iter=iter+1;
end

IFmset=IFsetiter(:,:,1:iter);
smset=ssetiter(:,:,:,1:iter);
IA=xim.^2+yim.^2;



















