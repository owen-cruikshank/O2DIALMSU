function delta = findDelta(alphaModel,alpha,Non,Noff,rm,ts)
dr = rm(2)-rm(1);
[~,startAlt]=min(abs(5900-rm));
startAlt

delta = zeros(startAlt,length(ts));
for ii=startAlt-1:-1:1
    delta(ii,:)=(2*dr*(alphaModel(ii,:)-alpha(ii,:))+delta(ii+1,:).*(1./Non(ii+1,:)-1./Noff(ii+1,:)))./(1./Non(ii,:)-1./Noff(ii,:));
end