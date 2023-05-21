function [Spec,f] = STFT(Sig,SampFreq,N,WinLen);


if (isreal(Sig))
Sig = hilbert(Sig);
end

SigLen = length(Sig);


WinLen = ceil(WinLen / 2) * 2;
t = linspace(-1,1,WinLen)';
sigma = 0.28;
WinFun = (pi*sigma^2)^(-1/4)*exp((-t.^2)/2/(sigma^2));

Lh = (WinLen - 1)/2; 

Spec = zeros(N,SigLen) ;  

for iLoop = 1:SigLen,
    
    tau = -min([round(N/2)-1,Lh,iLoop-1]):min([round(N/2)-1,Lh,SigLen-iLoop]);
    temp = floor(iLoop + tau);
    temp1 = floor(Lh+1+tau);
    rSig = Sig(temp);

    rSig = rSig .* conj(WinFun(temp1));
    Spec(1:length(rSig),iLoop) = rSig;
end;

Spec = fftshift(fft(Spec),1); 

[nLevel, SigLen] = size(Spec);

f = linspace(-SampFreq/2,SampFreq/2,nLevel);
t = (0: SigLen-1)/SampFreq;



