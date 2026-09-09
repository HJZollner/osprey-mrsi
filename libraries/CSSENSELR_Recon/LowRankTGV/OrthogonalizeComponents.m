function [ nV_tc] = OrthogonalizeComponents(V_tc)
    % apply the Graham Schmidt orthogonalization to the V_tc components

    NbComp=size(V_tc,2);
    nV_tc=V_tc;

    for a=2:NbComp
        %nV_tc(:,a)=nV_tc(:,a)-nV_tc(:,1:(a-1)) *diag(1.0./diag(nV_tc(:,1:(a-1))'*nV_tc(:,1:(a-1)))) *nV_tc(:,1:(a-1))'*V_tc(:,a);
        CoLinCoefMat= diag(1.0./diag(nV_tc(:,1:(a-1))'*nV_tc(:,1:(a-1)))) *nV_tc(:,1:(a-1))'*V_tc(:,a); %splitting the operation make it way faster !! (why????)
        nV_tc(:,a)=nV_tc(:,a)- ( nV_tc(:,1:(a-1)) * CoLinCoefMat );    
    end

end
