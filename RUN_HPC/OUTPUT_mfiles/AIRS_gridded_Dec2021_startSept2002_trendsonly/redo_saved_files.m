for JOB = 1 : 64

  JOB

  clear subjac
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_subjacLatBin' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_subjac_nostruct_LatBin' num2str(JOB,'%02i') '.mat'];
  loader = ['load ' foutsubjac];
  eval(loader);

  save(foutsubjac2,'-struct', 'subjac');
end
