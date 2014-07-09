function [outwave,out_array]=run7choose(gasID,fmin,fmax,profname,topts)

% same as run7.m, except can choose which quantum number stuff to compute!
% function [outwave,out_array]=run7(gasID,fmin,fmax,profname,{topts})
%
% if you straight away want 0.0025 res with no 5 pt averaging  
% also note the coarse freq point spacing is 0.5 ==> can do quad interp  
% also note the far wing coarse inclusion is 100 cm  
% [fr,k]=run6(3,705,755,0.0025,0.1,0.5,1,1,2,100,1,0.0,0.0,LVG,-1,profname);
%     or {topts.nbox = 1, topts.xfar = 100.0}
%     [fr,k]=run7(3,705,755,profname,topts);
% if you do 0.0005 res with 5 pt averaging to give 0.0025  
% also note the coarse freq point spacing is 1.0 ==> can do linear interp  
% [fr,k]=run6(3,705,730,0.0005,0.1,0.5,1,1,2,25,5,0.0,0.0,LVG,-1,profname);
%     [fr,k]=run7(3,705,755,profname);
%
%this function computes the lines spectra for gasID
%
%  TYPE     VAR           DESCRIPTION              TYPICAL VALUE
%----------------------------------------------------------------------
%            these need to be passed in by the user
%------------------------------------------------------------------------
%integer   gasID          HITRAN gas ID                   3
%
%integer   fmin           minimum freq (cm-1)            605
%integer   fmax           maximum freq (cm-1)            630
%
%matrix   profname       this contains N x 5 matrix : layernum,p,pp,t,x
%                         where column 1 = layer number 
%                               column 2 = layer pressure (atm)
%                               column 3 = layer partial pressure (atm)
%                               column 4 = layer temperature (k)
%                               column 5 = layer gas amt (kilomoles/cm2)
%------------------------------------------------------------------------
%  TYPE     VAR           DESCRIPTION                   DEFAULT VALUE
%            these are default values; can be reset using structure {topts}
%------------------------------------------------------------------------
%real      ffin           fine point spacing (cm-1)             0.0005
%real      fmed           medium point spacing (cm-1)           0.1
%real      fcor           coarse point spacing (cm-1)           0.5
%
%real      fstep          wide mesh width size (cm-1)           1.0
%real      xnear          near wing distance(cm-1)              1.0
%real      xmed           med wing distance(cm-1)               2.0
%real      xfar           far wing distance(cm-1)               25.0
%
%integer   nbox           boxcar sum size (odd integer)            5
%
%real      str_far        min line strength for far wing lines     0.0
%real      str_near       min line strength for near wing lines    0.0
%
%char      LVG            (L)orentz,Voi(G)t,(V)anHuber             'V'
%integer   CKD            continumm no (-1)                         -1
%                         yes water : (0,21,23) 
%                         yes O2,N2 : (+1) 
%string    HITRAN         path to HITRAN database   /asl/data/hitran/h2k.by.gas
%
%real      stren_mult     multiplies ALL strengths             1.0
%real      width_mult     multiplies ALL widths                1.0
%real      tsp_mult       multiplies ALL pressure shifts       1.0

%restrictions :
%(1) xnear <= xmed <= xfar        
%(2) xnear >= fstep
%(3) xmed/fmed  xnear/fmin   fstep/fmed   fstep/ffin       are integers
%(4)fstep/(nbox*ffin)        fcor/ffin                     are integers
%(5)(fmax-fmin)/fstep        (fmax-fmin)/fcor              are integers
%
%thus ffin,fmed,fcor are all 10/2^n   n <= 5

% HISTORY : 
% 12) introduced parameters stren_mult,width_mult, tsp_mult which multiplies
%      ALL strengths, foreign/self widths, pressure shifts     Apr 03, 2002
% 11) changed it to run7(gasID,v1,v2,profile,{topts})          Feb 25, 2002
% 10) added HITRAN parameter                                   Feb 25, 2002
% 9) profile file name, which contains 5 column atmosphere profile
% 8) continuum for water, O2,N2
% 7) allows 3 different line shapes : lorentz, voigt or van huber
% 6) differsn from run4.m in that the loop over lines is now MEXED!!!
% 5) allows the user to choose ALL isotopes of a gas (0) or one of the
%    isotopes
% 4) this differs from run3.m in that the code computes the optical depths
%    of all lines in EACH layer (k*length), and then does a UNION of all the
%    lines, so that lines do not turn on and off, driving the SVD crazy
% 3) we have parameter xnear xmed xfar to define fine,med,coarse meshes
% 2) better version of run2.m (vectorised as much as possible)
% 1) same as run1WORKS except this should compute the spectra for all 100
%    layers (run1WORKS is mainly a debugger for ONE layer)

% differs from run2a.m in that we have redefined what FINE,MEDIUM,COARSE
% lines to use. take [f3 f4] as the limits of the current wide mesh
% (a) fine lines are always within  [f3-xnear  f4+xnear]     xnear ~ 1 cm-1
%     ie these lines are always within +/- 1 cm-1 of current wide mesh
% (b) medium lines are always within 
%             [f3-xmed f3-xnear] U [f4+xnear f4+xmed]      xmed ~ 2 cm-1
%     ie these lines are always 1-2 cm-1 away from current wide mesh
% (c) far lines are always within 
%             [f3-xfar f3-xmed] U [f4+xmed f4+xfar]        xfar ~ 25 cm-1
%     ie these lines are always 2-25 cm-1 away from current wide mesh
% in addition, this means that the layer loop over coarse freq grid has been
% removed, so that there is only ONE loop over layers

%
% all the voigt, vhh fcns have been MEXed
%if brd > 1e-2 can use the Martin/Puerta voigt fcn (very fast and accurate)
%if brd < 1e-2 use the Tobin's voigt fcn (slow)

% wherever code says %%%%%%%%%%%%%% USER_WARN the user can reset some params

%%%%%%%%%% this uses MEX files for 
%%%%%%%% (1) boxint2        (boxint.f is the slow MATLAB file)
%%%%%%%% (2) vhh,vhh1       (vhh is faster than tobin's vhh1)

%real      ffin           fine point spacing (cm-1)      0.0005
%real      fmed           medium point spacing (cm-1)    0.1
%real      fcor           coarse point spacing (cm-1)    0.5
%
%real      fstep          wide mesh width size (cm-1)      1.0
%real      xnear          near wing distance(cm-1)         1.0
%real      xmed           med wing distance(cm-1)          2.0
%real      xfar           far wing distance(cm-1)          25.0
%
%integer   nbox           boxcar sum size (odd integer)    5
%
%real      str_far        min line strength for far wing lines     0.0
%real      str_near       min line strength for near wing lines    0.0
%
%char      LVG            (L)orentz,Voi(G)t,(V)anHuber    'V'
%integer   CKD            continumm no (-1)               -1
%                         yes water : (0,21,23) 
%                         yes O2,N2 : (+1) 
%string    HITRAN         path to HITRAN database   /asl/data/hitran/h2k.by.gas

%real      stren_mult     multiplies ALL strengths             1.0
%real      width_mult     multiplies ALL widths                1.0
%real      tsp_mult       multiplies ALL pressure shifts       1.0

% function [outwave,out_array]=run6(gasID,fmin,fmax,ffin,fmed,fcor,fstep,...
%  xnear,xmed,xfar,nbox,strength_far,strength_near,LVG,CKD,HITRAN,profname);
ffin          = 0.0005;
fmed          = 0.1;
fcor          = 0.5;
fstep         = 1.0;
xnear         = 1.0;
xmed          = 2.0;
xfar          = 25.0;
nbox          = 5;
strength_far  = 0.0;
strength_near = 0.0;
LVG           = 'V';
CKD           = -1;
HITRAN        = '/asl/data/hitran/h2k.by.gas';
stren_mult    = 1.0;
width_mult    = 1.0;
tsp_mult      = 1.0;

allowedparams = [{'ffin'},     {'fmed'},         {'fcor'},     {'fstep'}, ...
                 {'xnear'},    {'xmed'},         {'xfar'},                ...
                 {'nbox'},     {'strength_far'}, {'strength_near'},       ...
                 {'LVG'},      {'CKD'},          {'HITRAN'},              ...
              {'stren_mult'},   {'width_mult'},  {'tsp_mult'} ];


%read in {topts} 
if nargin == 5
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(allowedparams,optvar{i})) == 1)
     eval(sprintf('%s = topts.%s;', optvar{i}, optvar{i}));
   else
     fprintf(1,'topts param not in allowed list ... \n');
     error('quitting run7');
     end
   end
 end


% Howard's latest suggestion
% option to override the above defaults with topts fields
% (other options set with kopts are passed along to rtcalc)
%if nargin == 5
%  optvar = fieldnames(topts);
%  for i = 1 : length(optvar)
%    vname = optvar{i};
%    if exist(vname, 'var')
%      eval(sprintf('%s = topts.%s;', vname, vname));
%    else
%      fprintf(1,'topts param %s not in allowed list ... \n', vname);
%      error('quitting run7');    
%    end
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (gasID == 1)         %water
  if (CKD >= 0)
    docontinuum = -1;
    fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
    fprintf(1,'\n use run6water instead ... \n');
    error('exiting out of run6');
    end
  end

allowedline=['l','L','v','V','g','G'];
okdokey=intersect(allowedline,LVG);
if (length(okdokey) ~= 1)
  error('lineshape must be L or V or G');
  end

if ((LVG == 'l') | (LVG == 'L'))
  lvgNum=-1;
elseif ((LVG == 'v') | (LVG == 'V'))
  lvgNum=0;
elseif ((LVG == 'g') | (LVG == 'G'))
  lvgNum=1;
  end

%see that all parameters makes sense
z=checkpar(fstep,ffin,fmed,nbox,fcor,fmax,fmin,xnear,xmed,xfar);
if (z ~= 1) 
  error('there is an error in your parameters')
  end

%get all lines above strengthM value
if (strength_far > strength_near)
  strengthM=strength_near;
else
  strengthM=strength_far;
  end

%%%%%%%%%%%%%% USER_WARN
% set these parameters as they determine what mininum Optical Depth 
% should be used  === k(line strengh)*max(VHH)
OptDepth_close=1e-10;
OptDepth_far=OptDepth_close*1000;

%%%%%%%%%%%%%% USER_WARN
% adjust parameter "mult" to adjust how far away from fmin,fmax we
% want to consider lines to include from HITRAN
mult=1;
low= fmin-mult*xfar;
high=fmax+mult*xfar;

%%%%%%%%%%%%%%%%%%% LOAD IN ISOTOPE MASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
load mass.dat        %get the mass of the isotopes  
[x,y]=size(mass);
%first 32 are the numbers of the isotopes  <iGasID  NumIsotopes  0.0>
isotope_num=mass(1:32,2:2);
%mext bunch are the masses of the isotopes  <iGasID  MassIsotopes Abundance>
themass=mass(33:x,2:2);
%now pick up the gas masses we need
dummy1=sum(isotope_num(1:gasID-1));
dummy2=sum(isotope_num(1:gasID));
dummy=dummy1+1:dummy2;
mass_iso=themass(dummy);
liso=length(mass_iso);

%%%%%%%%%%%%%%%%%%% IS THIS INTERACTIVE SESSION OR CLUNK THRU PROFILE %%%%%%%%%
useruser= 1;
if (useruser > 0)
  which_isotope=input('Enter which isotope : ');

  do_load=0;
  MinLayer=1; MaxLayer=1; Step=1;        %use only ONE "layer"
  NumLayers=1;
  TheLayers=MinLayer:Step:MaxLayer;

  MGC=8.314674269981136  ;  
  press=input('Enter total pressure (in atm) : ');
  partpress=input('Enter gas partial pressure (in atm) : ');
  temperature=input('Enter temperature (in K) : ');
  GasAmt=input('Enter path cell length (in cm) ');
  %change to kmoles cm-2 
  GasAmt=GasAmt*101325*partpress/1e9/MGC/temperature; %change to kmoles/cm2 
  end

%%%%%%%%%%%%%%%%%%% LOAD IN GAS PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useruser < 0
  a=['load ' profname];
  eval(a);
  %find last occurence of '/' in profname
  ii=findstr(profname,'/');ii=ii(length(ii))+1;
  profile=profname(ii:length(profname));

  press=eval([profile '(:,2)']);
  partpress=eval([profile '(:,3)']);
  temperature=eval([profile '(:,4)']);
  GasAmt=eval([profile '(:,5)']);

  subplot(2,2,1); plot(1:length(GasAmt),GasAmt); 
  subplot(2,2,2); plot(1:length(GasAmt),temperature);
  subplot(2,2,3); plot(1:length(GasAmt),press);
  subplot(2,2,4); plot(1:length(GasAmt),partpress);pause(1);
  title('Gas Profile');
  clf
  end

which_isotope=0;                             %default use all isotopes
MinLayer=1; MaxLayer=length(GasAmt); Step=1; %default step thru ALL layers
NumLayers=(MaxLayer-MinLayer+Step)/Step;
TheLayers=MinLayer:Step:MaxLayer;

%%%%%%%%%%%%%%%%%% DEFINE OUTPUT ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define output freq array, output 100 layer matrix
%output point spacing grid
output_spacing=nbox*ffin;
%eg 605:0.0025:629.9975;
fmaxUSE   = fmax-output_spacing;
outwave   = fmin:output_spacing:fmaxUSE;  
lennnn = (1:length(outwave))-1; outwave = fmin+lennnn*output_spacing;   %<-----

out_array =zeros(NumLayers,length(outwave));

%%%%%%% READ IN RELEVANT DATA FROM HITRAN DATABASE  %%%%%%%%%%%%%%%%%%%%%%%%
%initialize Q fcns (partition fcns)
[A,B,C,D,G] = qtips(gasID,liso);

%% fnamePRE='/salsify/scratch4/h96.by.gas/g';        %H96 -- old
%% fnamePRE='/salsify/scratch4/h98.by.gas/g';        %H98 -- KCARTA database 
%% fnamePRE='/salsify/scratch4/h2k.by.gas/g';        %H98 -- KCARTA database 
if (HITRAN(length(HITRAN)) == '/')
  fnamePRE = [HITRAN 'g' ];
else
  fnamePRE = [HITRAN '/g'];
  end
fnamePOST='.dat';
fnameIN=int2str(gasID);
fname=[fnamePRE fnameIN fnamePOST]

%%%%%%%%%%ZZZ
%%%%%%%%following line if we want to use latest Hitran* database
[lineORIG]=hitread(low,high,strengthM,gasID,fname);
lineORIG.stren  = lineORIG.stren  * stren_mult;
lineORIG.tsp    = lineORIG.tsp    * tsp_mult;
lineORIG.abroad = lineORIG.abroad * width_mult;
lineORIG.sbroad = lineORIG.sbroad * width_mult;

fprintf(1,'before, number of lines = %6i \n',lineORIG.linct);

iL = input('Enter lower global quantum number : ');
iU = input('Enter upper global quantum number : ');
iSiS = input('Enter isotope : ');
ii=find((lineORIG.ilsgq==iL) & (lineORIG.iusgq==iU) & (lineORIG.iso==iSiS)); 
lineORIG.iso    = lineORIG.iso(ii);
lineORIG.wnum   = lineORIG.wnum(ii);
lineORIG.stren  = lineORIG.stren(ii);
lineORIG.tprob  = lineORIG.tprob(ii);
lineORIG.abroad = lineORIG.abroad(ii);
lineORIG.sbroad = lineORIG.sbroad(ii);
lineORIG.els    = lineORIG.els(ii);
lineORIG.abcoef = lineORIG.abcoef(ii);
lineORIG.tsp    = lineORIG.tsp(ii);
lineORIG.iusgq  = lineORIG.iusgq(ii);
lineORIG.ilsgq  = lineORIG.ilsgq(ii);
lineORIG.gasid  = lineORIG.gasid(ii);
lineORIG.linct  = length(lineORIG.iso);
fprintf(1,'after, number of lines = %6i \n',lineORIG.linct);

if gasID == 19
  lineORIG=reduce19(lineORIG);
  end

fprintf(1,'\n read in hittomat \n\n');
number_of_lines = lineORIG.linct;

nwide=(fmax-fmin)/fstep;    %number of wide meshes

%number of points in wide mesh (output resolution)
%eg if fstep=1.0, nbox*ffin=0.0025, then we really should have 
% fstep/(nbox*ffin)+1 = 401 pts
%but we are letting the "last" point be the "first" point of the next wide 
%mesh, so we only have 400 pts == fstep/(nbx*ffin) 
chunklen=fstep/(nbox*ffin); 
chunk=1:chunklen;

cpptotal=cputime; 

%now do a UNION of all lines so that lines do not turn on and off, turning 
%both SVD and Fast Models crazy. Also here include only the isotope the
%user specified
if ((NumLayers > 1) & (number_of_lines > 0))
  line=findUnion(lineORIG,GasAmt,temperature,press,partpress,A,B,C,D,G,...
             gasID,mass_iso,OptDepth_close,which_isotope,TheLayers);
  number_of_lines = line.linct;  
elseif ((NumLayers == 1) & (number_of_lines > 0))
  line=lineORIG;
  number_of_lines = line.linct;  
  end

if (number_of_lines <= 0)
  fprintf(1,'found NO lines in the HITRAN file \n');
else
  fprintf(1,'found %8i lines in the HITRAN file ... proceeding\n',line.linct);
  end


%do near wing calc
%%%%%%%%%%%%%%%%%% ZZZZZZ
if (number_of_lines > 0)
  for  ii=1:nwide    %OUTER LOOP OVER WIDE MESH
  %%for  ii=1:-1    %OUTER LOOP OVER WIDE MESH

    %%%%%%%%% define the fine,medium, coarse resolution wavenumber grids %%%%
    %this is defining the temporary fine freq resolution array
    %start freq of current wide mesh
    f1=fmin+(ii-1)*fstep - ffin*(nbox-1)/2;   
    %stop freq of current wide mesh  ... note we have  the -ffin*nbox because
    % we do not want to include the "last" point here
    f2=fmin+ii*fstep - (ffin*nbox) + ffin*(nbox-1)/2; 
    fine_freq=f1:ffin:f2;
    lennnn = (1:length(fine_freq))-1; fine_freq=f1+lennnn*ffin;        %<-----

    %this is defining the temporary medium freq resolution array
    f3=fmin+(ii-1)*fstep;   %start freq of current wide mesh
    f4=fmin+(ii)*fstep;     %stop freq of current wide mesh
    med_freq=f3:fmed:f4;
    lennnn = (1:length(med_freq))-1; med_freq=f3+lennnn*fmed;          %<-----

    %this is defining the temporary coarse freq resolution array
    f5=fmin+(ii-1)*fstep;   %start freq of current wide mesh
    f6=fmin+(ii)*fstep;     %stop freq of current wide mesh
    coarse_freq=f3:fcor:f4;
    lennnn = (1:length(coarse_freq))-1; coarse_freq=f3+lennnn*fcor;    %<-----

    %%%%%%%%% define the current output resolution wavenumber grids %%%%%%%
    %at the output wavenumber resolution, this is the current wavenumber array
    output_freq=f3:nbox*ffin:f4-nbox*ffin;
    lennnn = (1:length(output_freq))-1; coarse_freq=f3+lennnn*nbox*ffin;  %<--
    %occupying the following columns in the output array outwave,out_array
    thechunk=chunk+(ii-1)*chunklen;

    %%%%%%%% sort the lines into very near, medium near, far %%%%%%%%%%%%%%%
    [very_near]=sortbins(line,f3,f4,-1,xnear);
    [medium_near]=sortbins(line,f3,f4,xnear,xmed);
    [far_wing]=sortbins(line,f3,f4,xmed,xfar);

    fprintf(1,'mesh %3i start stop freq = %8.5f %8.5f \n',ii,f3,f4);

    for jj=MinLayer:Step:MaxLayer 
      %INNER LOOP OVER LAYERS 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;

      tempr=temperature(jj);
      p=press(jj);
      ps=partpress(jj);
      %%%%%%%%%%%%%%% do the FINE GRID first (very near points)
      %fprintf(1,'\n number of VERY NEAR lines = %3i',very_near.linct);
      if very_near.linct > 0
        qfcn=q(A,B,C,D,G,very_near,tempr);
        pwr=very_near.abcoef;
        for_brd=very_near.abroad;
        self_brd=very_near.sbroad;
        freq=very_near.wnum+press(jj)*very_near.tsp;  %freq shift
        energy=very_near.els;
        s0=very_near.stren;
        brd=broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
%%        fprintf(1,'%6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e \n',freq(1),brd(1),p,ps,for_brd(1),self_brd(1),pwr(1),tempr);
        strength=find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj));

        outvect=loop(very_near.iso,mass_iso,brd,strength,freq,fine_freq,...
                   tempr,very_near.linct,length(fine_freq),lvgNum);
        scum=boxint2(outvect,nbox);
        out_array(nn,thechunk)=out_array(nn,thechunk)+scum;
        end 

      %%%%%%%%%%%%%%% do the MEDIUM GRID second  (medium near points)
      %fprintf(1,'\n number of MEDIUM lines = %3i',medium_near.linct);

      clear y
      z=zeros(size(med_freq));

      if medium_near.linct > 0
        qfcn=q(A,B,C,D,G,medium_near,tempr);
        pwr=medium_near.abcoef;
        for_brd=medium_near.abroad;
        self_brd=medium_near.sbroad;
        freq=medium_near.wnum+press(jj)*medium_near.tsp;  %freq shift
        energy=medium_near.els;
        s0=medium_near.stren;
        brd=broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
        strength=find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj));

        z=loop(medium_near.iso,mass_iso,brd,strength,freq,med_freq,...
                    tempr,medium_near.linct,length(med_freq),lvgNum);

        %do a spline interpolation of med_freq onto output_freq 
        xspline=spline(med_freq,z,output_freq);      %linear not good enough
        out_array(nn,thechunk)=out_array(nn,thechunk)+xspline;

        end

      %%%%%%%%%%%%%%% do the COARSE GRID third  (farpoints)
      %fprintf(1,'\n number of FAR lines = %3i',far_wing.linct);

      clear y
      z=zeros(size(coarse_freq));

      if far_wing.linct > 0
        qfcn=q(A,B,C,D,G,far_wing,tempr);
        pwr=far_wing.abcoef;
        for_brd=far_wing.abroad;
        self_brd=far_wing.sbroad;
        freq=far_wing.wnum+press(jj)*far_wing.tsp;  %freq shift
        energy=far_wing.els;
        s0=far_wing.stren;
        brd=broad(p,ps,1.0,for_brd,self_brd,pwr,tempr,gasID);
        strength=find_stren(qfcn,freq,tempr,energy,s0,GasAmt(jj));

        z=loop(far_wing.iso,mass_iso,brd,strength,freq,coarse_freq,...
                   tempr,far_wing.linct,length(coarse_freq),lvgNum);
        %do a spline interpolation of coarse_freq onto output_freq 

        if (length(coarse_freq) == 2)  
          xspline=interp1(coarse_freq,z,output_freq);    %linear good enough 
        elseif (length(coarse_freq) == 3)  
          xspline=interp_quad(coarse_freq,z,output_freq);  %quad good enough 
        elseif (length(coarse_freq) > 3)  
          xspline=spline(coarse_freq,z,output_freq);  %cubic good enough 
          end 

        out_array(nn,thechunk)=out_array(nn,thechunk)+xspline;
        end

      end           %inner loop over layers
    end             %outer loop over wide meshes (frequencies)
  end               %if number of lines > 0

%************************************************************************
%now do continuum .. it should be smooth enough that no need to do at
%fine mesh ... just do at output mesh
%       SUBROUTINE CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CON, CKD, whichlayer) 
%                    has been changed to
%      con = CALCON23( IDGAS, NFREQ, FREQ, FSTEP, NLAY, T, P, 
%     $    PARTP, AMNT, CKD, whichlayer)


docontinuum = -1;

if (gasID == 1)         %water
  if ((CKD == 0) | (CKD == 21) | (CKD == 23) | (CKD == 24))
    docontinuum = -1;
    fprintf(1,'\n this code WILL NOT DO water continuum ... \n');
    fprintf(1,'\n use run7water/ru76watercontinuum instead ... \n');
    error('exiting out of run7');
    end
  end

if (gasID == 7)         %o2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing O2 continuum ... \n');
    end
  end

if (gasID == 22)         %n2
  if (CKD == 1)
    docontinuum = 1;
    fprintf(1,'\n doing N2 continuum ... \n');
    end
  end

%%%%%%%%%%%%%%%%% this is the orig code .... might blow up if input vectors
%%%%%%%%%%%%%%%%  are tooo large
%if (docontinuum > 0)
%    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
%      nn=(jj-MinLayer)/Step+1;
%      scum=calcon(gasID,length(outwave),outwave,1,length(press),...
%                    temperature,press,partpress,GasAmt,CKD,jj);
%      out_array(nn,:)=out_array(nn,:)+scum;
%      end
%  end      

%******************************************************

%%%%from FORTRANFILES/max.inc
%c this is max length of arrays that can be used in the Mex Files  
%c this number came out of 
%c   200000 = max number of elements in mesh 
%c              eg (755-655)/0.0005 = 160000
%c        4 = number tacked on to arrays so boxint(y,5) can be done  
%      integer MaxLen
%      parameter(MaxLen=200010)
MaxLen=200010;

dff=ffin*nbox;
if (docontinuum > 0)
  if (length(outwave) <= MaxLen)          %can do it in one big chunk!
    for jj=MinLayer:Step:MaxLayer %INNER LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      scum=calcon(gasID,length(outwave),outwave,1,length(press),...
                  temperature,press,partpress,GasAmt,CKD,jj);
      out_array(nn,:)=out_array(nn,:)+scum;
      end

  else %have to break it into frequency intervals
    for jj=MinLayer:Step:MaxLayer %LOOP 1..100 = bottom -> top
      nn=(jj-MinLayer)/Step+1;
      index0=1:10000;
      number=floor(length(outwave)/10000);

      for kk=1:number             %LOOP OVER FREQUENCY INTERVALS
        index=index0+(kk-1)*10000;
        scum=calcon(gasID,length(index),outwave(index),1,...
            length(press),temperature,press,partpress,GasAmt,CKD,jj);
        out_array(nn,index)=out_array(nn,index)+scum;
        end

      %%see if anything left over
      if (index(length(index))+1 <= length(outwave)) 
        index=index(length(index))+1:length(outwave);
        scum=calcon(gasID,length(index),outwave(index),1,...
              length(press),temperature,press,partpress,GasAmt,CKD,jj);
        out_array(nn,index)=out_array(nn,index)+scum;
        end

      end
    end            %%%%if then else loop
  end

%******************************************************
cpptotal=cputime-cpptotal
