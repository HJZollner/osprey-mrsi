function VOIMask = MakeVOIMask(mrsiReconParams);
 
  SizeGrid =[  mrsiReconParams.mrProt.Ncol; mrsiReconParams.mrProt.Nlines;mrsiReconParams.mrProt.Nslc];
  
  SizeVOI =[  mrsiReconParams.mrProt.VOIWidth; mrsiReconParams.mrProt.VOIHeight; mrsiReconParams.mrProt.VOI3D];
  SizeFoV =[  mrsiReconParams.mrProt.FoVWidth; mrsiReconParams.mrProt.FoVHeight; mrsiReconParams.mrProt.FoV3D];
  
  VoiGrid = round(SizeVOI./SizeFoV.*SizeGrid/2)*2;
  VOIMask = 0*mrsiReconParams.ImMask;
  VOIMask ( (1+SizeGrid(1)/2-VoiGrid(1)/2 ):(SizeGrid(1)/2+VoiGrid(1)/2 ),(1+SizeGrid(2)/2-VoiGrid(2)/2 ):(SizeGrid(2)/2+VoiGrid(2)/2 ),(1+SizeGrid(3)/2-VoiGrid(3)/2 ):(SizeGrid(3)/2+VoiGrid(3)/2 ) ) = 1;
  
  end
