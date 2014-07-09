% this script writes out this code
%      PARAMETER (kNumkCompT=188)
%      DATA kaTag        / 1,   2,   3,   4,   5,   6,   7,  8,   9   /
%      DATA kaNumkComp   / 08, 89,  09,  20,  24,  16,  20,  1,   1   /
%      DATA kaCTag       /'q', 'r', 'm', 'n', 'v', 'w', 's', 'k', 'p' /
%      DATA kaMinFr      /0500.00000,   0605.00000,    4050.00000,  5000.00000,
%     $                   8000.00000,   14000.00000,   2830.00000, 
%     $                   0200.00000,   0355.00000/
%      DATA kaMaxFr     /0605.00000,   2830.00000,    4950.00000,   8000.00000,
%     $                   14000.00000,  22000.00000,   3305.00000,
%     $                    0355.00000,   0500.00000/
%      DATA kaFrStep  /1.500000e-3,  2.500000e-3,   1.000000e-2,  1.500000e-2,
%     $                 2.500000e-2,  5.000000e-2,   2.500000e-3,
%     $                 0.500000e-3,  1.000000e-3/
%      DATA kaBlSize  /015.00000,    025.00000,     100.00000,    150.00000,
%     $                   250.00000,    500.00000,   025.00000,
%     $                   005.00000,    010.00000/
%      DATA kaFineFrStep  /3.000000d-4, 5.000000d-4,  1.000000d-3, 1.000000d-3,
%     $                     1.000000d-3, 1.000000d-3,  5.000000d-4,
%     $                     0.500000d-4, 1.000000d-4/
%      DATA kaMediumFrStep/1.000000d-1, 1.000000d-1,  1.000000d-1, 1.000000d-1,
%     $                     1.000000d-1, 1.000000d-1,  1.000000d-1,
%     $                      1.000000d-1, 1.000000d-1/
%      DATA kaCoarseFrStep/3.000000d-1, 5.000000d-1,  5.000000d-1, 5.000000d-1,
%     $                     5.000000d-1, 5.000000d-1,  5.000000d-1,
%     $                     5.000000d-1,  5.000000d-1/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% adjust these as necessary
%%%%% set in stone in Oct 2009 is "r" stands for 605-2830 cm-1 = Tag 2
%%%%% and these boundary definitions have some freedom

%%% orig????
kaTag      = [1      2     3     4     5     6     7     8     9];
kaCTag     = {'q',  'r',  'm',  'n',  'v',  'w',  's',  'k',  'p'};
kaMinFr    = [00500 00605 04050 05000 08000 14000 02830 00200 00355];
kaMaxFr    = [00605 02830 04950 08000 14000 22000 03305 00355 00505];

%%% Nov 2009
kaTag      = [15    20    30    35    40    50    25    10    12   ];
kaCTag     = {'q',  'r',  'm',  'n',  'v',  'w',  's',  'k',  'p'  };
kaMinFr    = [00500 00605 04050 05000 08000 14000 02830 00200 00355];
kaMaxFr    = [00605 02830 04950 08000 14000 22000 03305 00355 00505];

%%% early Dec 2009
kaTag      = [15    20    30    35    40    50    25    10    12   ];
kaCTag     = {'q',  'r',  'm',  'n',  'v',  'w',  's',  'k',  'p'  };
kaMinFr    = [00500 00605 04050 05000 08000 14000 02830 00140 00300];
kaMaxFr    = [00605 02830 04950 08000 14000 22000 03580 00310 00510];

%%% late Dec 2009
kaTag      = [15    20    30    35    40    50    25    10    12     08  ];
kaCTag     = {'q',  'r',  'm',  'n',  'v',  'w',  's',  'k',  'p',   'j' };
kaMinFr    = [00500 00605 04050 05000 08000 14000 02830 00140 00300  00080];
kaMaxFr    = [00605 02830 04950 08000 14000 22000 03580 00310 00510  00150];

%%% Jan 2010
%% range      <-- IR --> <SWIR> <----- NIR -----> <VIS> <UV>  <------- FIR --->
kaTag      = [15    20    25    30    35    40    50    55    12    10     08  ];
kaCTag     = {'q',  'r',  's',  'm',  'n',  'o',  'v',  'u',  'p',  'k',   'j'  };
kaMinFr    = [00500 00605 02830 03550 05550 08250 12000 25000 00300 00140  00080];
kaMaxFr    = [00605 02830 03580 05650 08400 12250 25000 44000 00510 00310  00150];

%%% Dec 2011
%% range      <-- IR --> <SWIR> <----- NIR -----> <VIS> <UV> XXXX <------- FIR ----------------------->
kaTag      = [15    20    25    30    35    40    50    55        12    10     08    06    04    02];
kaCTag     = {'q',  'r',  's',  'm',  'n',  'o',  'v',  'u',      'p',  'k',  'j',   'h',  'g',  'f'  };
kaMinFr    = [00500 00605 02830 03550 05550 08250 12000 25000     00300 00140 00080 00050 00030 00015];
kaMaxFr    = [00605 02830 03580 05650 08400 12250 25000 44000     00510 00310 00150 00080 00050 00030];

squareavg = 5;
numpts = 10000;
%% this is the "center freq" used to define the resolution
kaUseJunk = [00500 00800 02500 04000 06000 08000 14000 40000      00300 00200 00100 00075 00040 00025];

disp('---> output saved in dump_kaTag_for_predefined.param')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%         code below this does NOT need to be altered     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/SPECTRA
disp('running runXtopts_params_smart')
for ii = 1 : length(kaTag)
  topts = runXtopts_params_smart(kaUseJunk(ii));
  kaFrStep(ii) = squareavg*topts.ffin; 
  kaBlSize(ii) = squareavg*numpts*topts.ffin; 
  kaFineFrStep(ii) = topts.ffin; 
  kaMediumFrStep(ii) = topts.fmed; 
  kaCoarseFrStep(ii) = topts.fcor; 
  kaNumkComp(ii) = (kaMaxFr(ii)-kaMinFr(ii))/kaBlSize(ii);
  end

all_ints = abs(kaNumkComp-round(kaNumkComp));
for ii = 1 : length(kaTag)
  if all_ints(ii) < 1e-5
    all_ints(ii) = 0;
    end
  end
if sum(all_ints ~= 0)
  aax = [kaMinFr' kaMaxFr' kaNumkComp'];
  disp('kaMinFr kaMaxFr kaNumkComp');
  fprintf(1,'%10.4f    %10.4f     %10.6f \n',aax');
  error('whoops kaNumkComp not all integers ... better check kaMinFr,kaMaxFr,kaBlSize')
  end

kNumkCompT = sum(kaNumkComp)

[Y,I] = sort(kaMinFr);
%I = 1 : length(I);

fid = fopen('dump_kaTag_for_predefined.param','w');
fprintf(fid,'      PARAMETER (kNumkCompT = %4i ) \n',kNumkCompT);
fprintf(fid,'      DATA kaTag        /');
for jj = 1 : length(kaTag)
  if kaTag(I(jj)) < 10
    fprintf(fid,'  0%1i',kaTag(I(jj)));
  else
    fprintf(fid,'  %2i',kaTag(I(jj)));
    end
  if jj < length(kaTag);   fprintf(fid,','); end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaNumkComp   /');
for jj = 1 : length(kaNumkComp)
  if kaNumkComp(I(jj)) < 10
    fprintf(fid,'  0%1i',round(kaNumkComp(I(jj))));
  else 
    fprintf(fid,'  %2i',round(kaNumkComp(I(jj))));
    end
  if jj < length(kaTag);   fprintf(fid,','); end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaCTag       /');
for jj = 1 : length(kaCTag)
  fprintf(fid,' ''%s''',char(kaCTag(I(jj))));
  if jj < length(kaTag);   fprintf(fid,','); end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaMinFr        /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4f',kaMinFr(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaMaxFr        /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4f',kaMaxFr(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaFrStep       /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4e',kaFrStep(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaBlSize       /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4f',kaBlSize(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaFineFrStep   /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4e',kaFineFrStep(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaMediumFrStep /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4e',kaMediumFrStep(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');

fprintf(fid,'      DATA kaCoarseFrStep /');
for jj = 1 : length(kaCTag)
  fprintf(fid,'  %8.4e',kaCoarseFrStep(I(jj)));
  if jj < length(kaTag);   fprintf(fid,','); end
  if mod(jj,4) == 0
     fprintf(fid,'\n'); 
     fprintf(fid,'     $                  ');
     end
  end
fprintf(fid,'  / \n');


fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
