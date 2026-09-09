function Oper_rk = ForwDFT3_Masked_Operator(ImSize,MaskK,MaskI)

% Input Sizes
N1 = ImSize(1);
N2 = ImSize(2);
N3 = ImSize(3);



[InTraj_x, InTraj_y,InTraj_z] = ndgrid((0:(N1-1)),(0:(N2-1)),(0:(N3-1)));%ndgrid( (1:N1),(1:N2),(1:N3));
InTraj_r3 = [InTraj_x(MaskI(:)>0), InTraj_y(MaskI(:)>0),InTraj_z(MaskI(:)>0)];

clear InTraj_x InTraj_y InTraj_z

[OutTraj_x, OutTraj_y,OutTraj_z] = ndgrid( (0:(N1-1))/N1,(0:(N2-1))/N2,(0:(N3-1))/N3);
OutTraj_k3 = [OutTraj_x(MaskK(:)>0), OutTraj_y(MaskK(:)>0),OutTraj_z(MaskK(:)>0)];

clear OutTraj_x OutTraj_y OutTraj_z

InTraj_r3 = single(InTraj_r3);
OutTraj_k3 = single(OutTraj_k3);

Oper_rk = exp(-2*pi*1i*InTraj_r3*transpose(OutTraj_k3));

% Normalization
Oper_rk = Oper_rk / sqrt(size(InTraj_r3,1)); 

