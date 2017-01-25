function [hOut,pOut,nOut] = wrapper_levelsRTP_to_levelstext_and_lblrtmTAPE5OD(iProfile,fRTP,fTXT)

%% input
%%   iProfile = which profile in the RTP to process
%%   fRTP     = name of input LEVELS rtp file
%%   fTXT [optional] = name of text file to write out for kCARTA

%{
example
fRTP = '/home/sergio/KCARTA/IP_PROFILES/junk49.ip.rtp';
iProfile = 1;
wrapper_levelsRTP_to_levelstext_and_lblrtmTAPE5OD(iProfile,fRTP);
%}

disp('WARNING : addpath to this dir, as the file fTXT and TAPE5x are printed locally here')

if nargin == 2
  fTXT = 'junk_levels.prf';
end

%% make text levels file for KCARTA
levelsRTP_to_levelstext(iProfile,fRTP,fTXT)

addpath /home/sergio/KCARTA/MATLAB/LBLRTM

%% make TAPE5 for LBLRTM
v1in = 605;
v2in = 630;
v2in = 605;
dvOut = 0.0025;
[hOut,pOut,nOut] = template_levelsRTP_to_tape5_ODonly(fRTP,iProfile,v1in,v2in,dvOut);
