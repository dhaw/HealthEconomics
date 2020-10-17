function f=hcFinalSize
C=[1.7112    0.4927    0.4049    0.0293;
    0.3919    2.5527    0.4753    0.0348;
    0.3234    0.5548    0.8996    0.0728;
    0.0528    0.1904    0.3744    0.3830];
h=[0.0018,0.0030,0.0582,0.1348]';
NNbar=[4064198,12192593,36577778,13005432]';
NNsum=sum(NNbar);
A=C;%.*repmat(NNbar',4,1);

R0=(0:.01:2);
lr=length(R0);
Z=zeros(lr,1);
x0=ones(4,1);

for i=1:lr
    fun=@(z)(1-z-exp(-R0(i)*A*z));
    zi=fsolve(fun,x0);
    Z(i)=sum(NNbar.*zi);%/NNsum;%NNbar.*h.*
end

f=Z;
plot(R0,Z);
xlabel('R0')
ylabel('Hospitalisations')

    
