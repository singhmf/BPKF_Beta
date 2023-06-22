function[]=BPKF_Plot_Err(Out,SmthSpace,varargin)

if nargin==1
    disp('Assuming no smoothing kernel')
    SmthSpace=1;
end

if isfield(Out,'kalCV')
    if isfield(Out,'batchCV')
        batchCV=Out.batchCV;
    else
        disp('Assuming CV rate is 500')
        batchCV=(1:size(Out.kalCV,2))*500;
    end
    doCV=true;
else
    doCV=false;
end


if numel(varargin)==1
    varargin=repelem(varargin,1,2);
end

figure
subplot(1,2,1)
plot(SmthSpace:size(Out.kalE,2),convn(Out.kalE,ones(1,SmthSpace)/SmthSpace,'valid')');
set(gca,'ColorOrderIndex',1)
if doCV
    hold on;plot(batchCV,Out.kalCV,'.')
end
if ~isempty(varargin)&& strcmpi(varargin{1},'y')
    set(gca,'yscale','log')
end
title('Kal Error')
pbaspect([1 1 1])

subplot(1,2,2)
plot(SmthSpace:size(Out.recE,2),convn(Out.recE,ones(1,SmthSpace)/SmthSpace,'valid')')
set(gca,'ColorOrderIndex',1)
if doCV
    hold on;plot(batchCV,Out.recCV,'.')
end
if ~isempty(varargin)&& strcmpi(varargin{2},'y')
    set(gca,'yscale','log')
end
title('Rec Error')
pbaspect([1 1 1])
end