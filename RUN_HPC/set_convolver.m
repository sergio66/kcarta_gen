if ~exist('iDoConvolve')
  iDoConvolve = +1;
end

%%%%% %%%% iDoConvolve = -1;   %% fast no convolve so we can test eg DISORT vs PCLSAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iInstr')
  iInstr = 1;   % AIRS only
  iInstr = 2;   % IASI only
  iInstr = 4;   % CRIS all hi/CHIRP/CrIS lo
  iInstr = 14;  % AIRS + CRIS all hi/CHIRP/CrIS lo
  iInstr = 124; % AIRS + IASI + CRIS all hi/CHIRP/CrIS lo
  
  iDoConvolve = +1; iInstr = 14;  % AIRS + CRIS hi
  iDoConvolve = +1; iInstr = 124; % AIRS + IASI + CRIS hi/CHIRP/CrIS lo
  
  iDoConvolve = +1; iInstr = 1;   % AIRS
  iDoConvolve = +1; iInstr = 2;   % IASI     only
  iDoConvolve = +1; iInstr = 4;   % CRIS all
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iDoConvolve = +1; iInstr = 124; % AIRS + IASI + CRIS hi/CHIRP/CrIS lo
  iDoConvolve = +1; iInstr = 1;   % AIRS
end
  
