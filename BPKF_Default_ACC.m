function[Out,SimDat]=BPKF_Default_ACC(SimFile)

SimDat=load(SimFile);
H=SimDat.H;
[~,nX]=size(H);
MeasSet={SimDat.Meas};
Xguess={randn(nX,size(SimDat.Meas,2))/100};
Wmask=zeros(nX);
Wmask(SimDat.W~=0)=nan;


KalSpec.KFtype='E';
KalSpec.doEstack='n';
KalSpec.BFcov='Sim';%Sim';
    KalSpec.SimLength=200;
    KalSpec.nSim=40;
    KalSpec.nSaveStart=200;

KalSpec.nStep=20;
KalSpec.nRec=1;
KalSpec.minStep=5;
KalSpec.ErrMat='Mahal';
KalSpec.Pbase=eye(nX);
KalSpec.Pfix=eye(nX);
KalSpec.decP=.8;
KalSpec.decFix=.6;

ModelSpec.Hkron=1;
ModelSpec.H=SimDat.H;
ModelSpec.R=SimDat.R;
ModelSpec.Q=SimDat.Q;
ModelSpec.freeQ='n';
ModelSpec.fixC=SimDat.C;
ModelSpec.fixS=1;
ModelSpec.fixV=0;
ModelSpec.nRep=1;
if size(SimDat.D,2)~=1
    ModelSpec.fixD=diag(SimDat.D);
else
    ModelSpec.fixD=SimDat.D;
end
ModelSpec.Wmask=Wmask;


ParStr.L2err=0;
ParStr.CovErr=.1;
ParStr.BatchSz=1;
ParStr.NBatch=50000;

GradSpec.reg=.001;
GradSpec.dec1=.98;GradSpec.dec2=.95;
GradSpec.Rate=2.5*10^-4;

%GradSpec.Clip=5;
%GradSpec.ClipLength=5000;
%GradSpec.ClipType='Start';



Out=BPKF_Full(Xguess,MeasSet,ParStr,KalSpec,ModelSpec,GradSpec);
end