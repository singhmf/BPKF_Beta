function[Out]=BPKF_Sim(ooP,Start,DnSamp,tEnd,doPlot)

%% HARD-CODED for speed
W=ooP.Param{1};
D=ooP.Param{2};
C=ooP.Param{3};
S=ooP.Param{4};
V=ooP.Param{5};

nX=size(Start,1);
nDat=size(Start,2);
dt=1;
tVec=0:(dt*DnSamp):tEnd;
nT=numel(tVec);

dW=diag(ooP.rtQ)/10;

if nDat==1
    Out=nan(nX,nT);
    Out(:,1)=Start;
else
    Out=nan(nX,nT,nDat);
    Out(:,1,:)=Start;
end
if mean(dW(:)~=0)~=0
    if nDat==1
    Noise=dW.*randn(nX,tEnd/dt);
    else
    Noise=dW.*randn(nX,nDat,tEnd/dt);
    end
    
    
if nDat==1
    for i=1:(tEnd/dt)
        Start=single(W*double(tanh(S.*Start+V)))+D.*Start+C+Noise(:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        Start=W*tanh(S.*Start+V)+D.*Start+C+Noise(:,:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end

else

if nDat==1
    for i=1:(tEnd/dt)
        Start=W*tanh(S.*Start+V)+D.*Start+C;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        Start=W*tanh(S.*Start+V)+D.*Start+C;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end
end

if ~isempty(doPlot)&&strcmpi(doPlot(1),'y')
    figure
    if nDat==1
        plot(tVec,Out);
    else
        plot(tVec,squeeze(Out(:,:,1)));
    end
end
end