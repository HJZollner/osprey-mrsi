function  [New_mrsiData_ctkkk,New_Traj_as3,mrsiReconParams] = InterpFIDNbPts(mrsiData_ctkkk,mrsiReconParams,GoalPts)

[NbCoil NbTime M N Slc] = size(mrsiData_ctkkk);

mrsiData_tckkk = permute(mrsiData_ctkkk,[2,1,3,4,5]);
clear mrsiData_ctkkk
Xd=1:(NbTime);
Xm=linspace(1,(NbTime),(GoalPts));

NewNbTime=GoalPts;

fprintf(['Reducing the number of points in the FID from ',num2str(NbTime),' to ',num2str(NbTime),' ...']);

mrsiData_tckkk=interp1(Xd,mrsiData_tckkk,Xm,'spline');

New_mrsiData_ctkkk = permute(mrsiData_tckkk,[2,1,3,4,5]);

