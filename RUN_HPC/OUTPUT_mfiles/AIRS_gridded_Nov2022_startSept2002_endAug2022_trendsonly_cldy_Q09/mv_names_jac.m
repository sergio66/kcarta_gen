for JOB = 1 : 64
  xfoutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_20yr_CLD_Q09_'           num2str(JOB,'%02i') '.mat'];
  xfoutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(JOB,'%02i') '.mat'];

  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjacLatBin_kCARTA_ERA5_20yr_CLD_Q09_'           num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(JOB,'%02i') '.mat'];

  mver = ['!mv ' xfoutsubjac  ' ' foutsubjac];  eval(mver)
  mver = ['!mv ' xfoutsubjac2 ' ' foutsubjac2]; eval(mver)
end
