function kcartajacdumpx = compare_matlab_f90_dumpjacoutput(kcmixjacdump,kcartajacdump,iLay,iChan);

%% x97 = compare_matlab_f90_dumpjacoutput(kcmixjacdump,kcartajacdump,2,8220);
%% x97 has matrices from kCARTA f90 put into 1:iNumLayer formt (as opposed to 1.. 100 format,
%%    ome of which are rows of zeros

%{
kcartajacdump = readkc_all_dumpjac(955,1);
kcmixjacdump = load('/home/sergio/git/kcarta/temp_jacobian_matrices_955.mat');
x97 = compare_matlab_f90_dumpjacoutput(kcmixjacdump,kcartajacdump,2,8220);
save ../WORK/f90kcartajacobian_into_x97layers_955.mat x97
%}

if nargin == 2
  iLay = 1;
  iChan = 1;
elseif nargin == 3  
  iChan = 1;
end

[mm,nn]  = size(kcmixjacdump.raaLay2Gnd');
w = kcartajacdump.w;

iF90Lay = 100 - kcartajacdump.iNumLayer + 1;
iF90Lay = iF90Lay + (iLay-1);

kcartajacdumpx = struct;

%% for jacTG,jacQG, do 4 : 100
allF90laysJ = 100-kcartajacdump.iNumLayer+1 : 100;
kcartajacdumpx.w = kcartajacdump.w;
kcartajacdumpx.jacTG    = kcartajacdump.jacTG(:,allF90laysJ);
kcartajacdumpx.jacQG    = kcartajacdump.jacQG(:,allF90laysJ);
kcartajacdumpx.rawjacQG = kcartajacdump.rawjacQG(:,allF90laysJ);
kcartajacdumpx.rawjacTG = kcartajacdump.rawjacTG(:,allF90laysJ);
kcartajacdumpx.absc     = kcartajacdump.absc(:,allF90laysJ);

%% for all others, stick to i:iProfLayer
allF90lays0 = 1 : kcartajacdump.iNumLayer;
kcartajacdumpx.raaLay2Gnd = kcartajacdump.raaLay2Gnd(:,allF90lays0);
kcartajacdumpx.raaLay2Sp  = kcartajacdump.raaLay2Sp(:,allF90lays0);
kcartajacdumpx.raaGeneral = kcartajacdump.raaGeneral(:,allF90lays0);
kcartajacdumpx.raaOneMinusTau = kcartajacdump.raaOneMinusTau(:,allF90lays0);
kcartajacdumpx.raaTau         = kcartajacdump.raaTau(:,allF90lays0);
kcartajacdumpx.raaRad         = kcartajacdump.raaRad(:,allF90lays0);
kcartajacdumpx.raaRadDT       = kcartajacdump.raaRadDT(:,allF90lays0);
kcartajacdumpx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('The main ingredients raaLay2Gnd,raaGeneral,raaOneMinusTa,raaTau,raaRad are in form do iLay = 1 :  97     and do NO NEED TO be mapped to 4:100')

figure(1); plot(w,kcartajacdumpx.raaLay2Gnd(:,iLay),'cx-',w,kcartajacdump.raaLay2Gnd(:,iLay),'b.-',w,kcmixjacdump.raaLay2Gnd(:,iLay),'r');       title('raaLay2Gnd');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaLay2Gnd(iChan,:),1:mm,'cx-',kcartajacdump.raaLay2Gnd(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaLay2Gnd(iChan,:),1:mm,'r'); title('raaLay2Gnd');
  xlim([nanmin(kcmixjacdump.raaLay2Gnd(iChan,:)) nanmax(kcmixjacdump.raaLay2Gnd(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaLay2Gnd(:,iLay)./kcartajacdump.raaLay2Gnd(:,iLay),'b',w,kcmixjacdump.raaLay2Gnd(:,iLay)./kcartajacdump.raaLay2Gnd(:,iLay),'r'); title('f90revised/f90true and kc/f90true Lay2Gnd');
figure(4); plot(kcartajacdumpx.raaLay2Gnd(iChan,:)./kcmixjacdump.raaLay2Gnd(iChan,:),1:mm); title('f90revised/kc Lay2Gnd');
printarray([(1:10); kcartajacdumpx.raaLay2Gnd(iChan,1:10); kcartajacdump.raaLay2Gnd(iChan,1:10); kcmixjacdump.raaLay2Gnd(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaLay2Sp(:,iLay),'cx-',w,kcartajacdump.raaLay2Sp(:,iLay),'b.-',w,kcmixjacdump.raaLay2Sp(:,iLay),'r');       title('raaLay2Sp');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaLay2Sp(iChan,:),1:mm,'cx-',kcartajacdump.raaLay2Sp(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaLay2Sp(iChan,:),1:mm,'r'); title('raaLay2Sp');
  xlim([nanmin(kcmixjacdump.raaLay2Sp(iChan,:)) nanmax(kcmixjacdump.raaLay2Sp(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaLay2Sp(:,iLay)./kcartajacdump.raaLay2Sp(:,iLay),'b',w,kcmixjacdump.raaLay2Sp(:,iLay)./kcartajacdump.raaLay2Sp(:,iLay),'r'); title('f90revised/f90true and kc/f90true Lay2Sp');
figure(4); plot(kcartajacdumpx.raaLay2Sp(iChan,:)./kcmixjacdump.raaLay2Sp(iChan,:),1:mm); title('f90revised/kc Lay2Sp');
printarray([(1:10); kcartajacdumpx.raaLay2Sp(iChan,1:10); kcartajacdump.raaLay2Sp(iChan,1:10); kcmixjacdump.raaLay2Sp(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaGeneral(:,iLay),'cx-',w,kcartajacdump.raaGeneral(:,iLay),'b.-',w,kcmixjacdump.raaGeneral(:,iLay),'r');       title('General');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaGeneral(iChan,:),1:mm,'cx-',kcartajacdump.raaGeneral(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaGeneral(iChan,:),1:mm,'r'); title('General');
  xlim([nanmin(kcmixjacdump.raaGeneral(iChan,:)) nanmax(kcmixjacdump.raaGeneral(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaGeneral(:,iLay)./kcartajacdump.raaGeneral(:,iLay),'b',w,kcmixjacdump.raaGeneral(:,iLay)./kcartajacdump.raaGeneral(:,iLay),'r'); title('f90revised/f90true and kc/f90true General');
figure(4); plot(kcartajacdumpx.raaGeneral(iChan,:)./kcmixjacdump.raaGeneral(iChan,:),1:mm); title('f90revised/kc General');
printarray([(1:10); kcartajacdumpx.raaGeneral(iChan,1:10); kcartajacdump.raaGeneral(iChan,1:10); kcmixjacdump.raaGeneral(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaTau(:,iLay),'cx-',w,kcartajacdump.raaTau(:,iLay),'b.-',w,kcmixjacdump.raaTau(:,iLay),'r');       title('raaTau');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaTau(iChan,:),1:mm,'cx-',kcartajacdump.raaTau(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaTau(iChan,:),1:mm,'r'); title('raaTau');
  xlim([nanmin(kcmixjacdump.raaTau(iChan,:)) nanmax(kcmixjacdump.raaTau(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaTau(:,iLay)./kcartajacdump.raaTau(:,iLay),'b',w,kcmixjacdump.raaTau(:,iLay)./kcartajacdump.raaTau(:,iLay),'r'); title('f90revised/f90true and kc/f90true Tau');
figure(4); plot(kcartajacdumpx.raaTau(iChan,:)./kcmixjacdump.raaTau(iChan,:),1:mm); title('f90revised/kc Tau');
printarray([(1:10); kcartajacdumpx.raaTau(iChan,1:10); kcartajacdump.raaTau(iChan,1:10); kcmixjacdump.raaTau(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaOneMinusTau(:,iLay),'cx-',w,kcartajacdump.raaOneMinusTau(:,iLay),'b.-',w,kcmixjacdump.raaOneMinusTau(:,iLay),'r');       title('raaOneMinusTau');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaOneMinusTau(iChan,:),1:mm,'cx-',kcartajacdump.raaOneMinusTau(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaOneMinusTau(iChan,:),1:mm,'r'); title('raaOneMinusTau');
  xlim([nanmin(kcmixjacdump.raaOneMinusTau(iChan,:)) nanmax(kcmixjacdump.raaOneMinusTau(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaOneMinusTau(:,iLay)./kcartajacdump.raaOneMinusTau(:,iLay),'b',w,kcmixjacdump.raaOneMinusTau(:,iLay)./kcartajacdump.raaOneMinusTau(:,iLay),'r'); title('f90revised/f90true and kc/f90true OneMinusTau');
figure(4); plot(kcartajacdumpx.raaOneMinusTau(iChan,:)./kcmixjacdump.raaOneMinusTau(iChan,:),1:mm); title('f90revised/kc OneMinusTau');
printarray([(1:10); kcartajacdumpx.raaOneMinusTau(iChan,1:10); kcartajacdump.raaOneMinusTau(iChan,1:10); kcmixjacdump.raaOneMinusTau(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaRad(:,mm),'cx-',w,kcartajacdump.raaRad(:,100),'b.-',w,kcmixjacdump.raaRad(:,mm),'r');       title('raaRad at TOA');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaRad(iChan,:),1:mm,'cx-',kcartajacdump.raaRad(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaRad(iChan,:),1:mm,'r'); title('raaRad');
  xlim([nanmin(kcmixjacdump.raaRad(iChan,:)) nanmax(kcmixjacdump.raaRad(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaRad(:,iLay)./kcartajacdump.raaRad(:,iLay),'b',w,kcmixjacdump.raaRad(:,iLay)./kcartajacdump.raaRad(:,iLay),'r'); title('f90revised/f90true and kc/f90true Rad');
figure(4); plot(kcartajacdumpx.raaRad(iChan,:)./kcmixjacdump.raaRad(iChan,:),1:mm); title('f90revised/kc Rad');
printarray([(1:10); kcartajacdumpx.raaRad(iChan,1:10); kcartajacdump.raaRad(iChan,1:10); kcmixjacdump.raaRad(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.raaRadDT(:,mm),'cx-',w,kcartajacdump.raaRadDT(:,100),'b.-',w,kcmixjacdump.raaRadDT(:,mm),'r');       title('raaRadDT at TOA');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.raaRadDT(iChan,:),1:mm,'cx-',kcartajacdump.raaRadDT(iChan,allF90lays0),1:mm,'b.-',kcmixjacdump.raaRadDT(iChan,:),1:mm,'r'); title('raaRadDT');
  xlim([nanmin(kcmixjacdump.raaRadDT(iChan,:)) nanmax(kcmixjacdump.raaRadDT(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.raaRadDT(:,iLay)./kcartajacdump.raaRadDT(:,iLay),'b',w,kcmixjacdump.raaRadDT(:,iLay)./kcartajacdump.raaRadDT(:,iLay),'r'); title('f90revised/f90true and kc/f90true RadDT');
figure(4); plot(kcartajacdumpx.raaRadDT(iChan,:)./kcmixjacdump.raaRadDT(iChan,:),1:mm); title('f90revised/kc RadDT');
printarray([(1:10); kcartajacdumpx.raaRadDT(iChan,1:10); kcartajacdump.raaRadDT(iChan,1:10); kcmixjacdump.raaRadDT(iChan,1:10)]')
disp('ret to continue'); pause

disp('allF90lays0 is from 1:97')
[allF90lays0(1) allF90lays0(end)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Therefore only absc and raw jacTG,jacTQ need to be mapped do iLay = 1 :  97      iL = iLay+3, so 4 : 100')

figure(1); plot(w,kcartajacdumpx.absc(:,iLay),'cx-',w,kcartajacdump.absc(:,iF90Lay),'b.-',w,kcmixjacdump.absc(:,iLay),'r');       title('absc');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.absc(iChan,:),1:mm,'cx-',kcartajacdump.absc(iChan,allF90laysJ),1:mm,'b.-',kcmixjacdump.absc(iChan,:),1:mm,'r'); title('absc');
  xlim([nanmin(kcmixjacdump.absc(iChan,:)) nanmax(kcmixjacdump.absc(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,kcartajacdumpx.absc(:,iLay)./kcartajacdump.absc(:,iF90Lay),'b',w,kcmixjacdump.absc(:,iLay)./kcartajacdump.absc(:,iLay),'r'); title('f90revised/f90true and kc/f90true absc');
figure(4); plot(kcartajacdumpx.absc(iChan,:)./kcmixjacdump.absc(iChan,:),1:mm); title('f90revised/kc absc');
printarray([(1:10); kcartajacdumpx.absc(iChan,1:10); kcartajacdump.absc(iChan,1:10); kcmixjacdump.absc(iChan,1:10)]')
disp('ret to continue'); pause

figure(1); plot(w,kcartajacdumpx.jacTG(:,iLay),'cx-',w,kcartajacdump.jacTG(:,iF90Lay),'b.-',w,kcmixjacdump.jacTG(:,iLay),'r');       title('raw raajacTG');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(kcartajacdumpx.jacTG(iChan,:),1:mm,'cx-',kcartajacdump.jacTG(iChan,allF90laysJ),1:mm,'b.-',kcmixjacdump.jacTG(iChan,:),1:mm,'r'); title('raw raajacTG');
  xlim([nanmin(kcmixjacdump.jacTG(iChan,:)) nanmax(kcmixjacdump.jacTG(iChan,:))])
  legend('f90revised','f90true','matlab','location','best');  
figure(3); plot(w,kcartajacdumpx.jacTG(:,iLay)./kcartajacdump.jacTG(:,iF90Lay),'b',w,kcmixjacdump.jacTG(:,iLay)./kcartajacdump.jacTG(:,iF90Lay),'r'); title('f90revised/f90true and kc/f90true raw raajacTG');
figure(4); plot(kcartajacdumpx.jacTG(iChan,:)./kcmixjacdump.jacTG(iChan,:),1:mm); title('f90revised/kc jacTG');
printarray([(1:10); kcartajacdumpx.jacTG(iChan,1:10); kcartajacdump.jacTG(iChan,1:10); kcmixjacdump.jacTG(iChan,1:10)]')
disp('ret to continue'); pause

wonkF90  = kcartajacdump.jacQG;
wxonkF90 = kcartajacdumpx.jacQG;
wonkMat  = squeeze(kcmixjacdump.jacQG(1,:,:));
figure(1); plot(w,wxonkF90(:,iLay),'cx-',w,wonkF90(:,iF90Lay),'b.-',w,wonkMat(:,iLay),'r');       title('raw raajacQG');
  legend('f90revised','f90true','matlab','location','best');
figure(2); plot(wxonkF90(iChan,:),1:mm,'cx-',wonkF90(iChan,allF90laysJ),1:mm,'b.-',wonkMat(iChan,:),1:mm,'r'); title('raw raajacQG');
  legend('f90revised','f90true','matlab','location','best');
figure(3); plot(w,wxonkF90(:,iLay)./wonkF90(:,iF90Lay),'b',w,wonkMat(:,iLay)./wonkF90(:,iF90Lay),'r'); title('f90revised/f90true and kc/f90true raw raajacQG');
figure(4); plot(wxonkF90(iChan,:)./wonkMat(iChan,:),1:mm,'k'); title('f90revised/f90true jacQ');  
printarray([(1:10); wxonkF90(iChan,1:10); wonkF90(iChan,1:10); kcmixjacdump.jacTG(iChan,1:10)]')
disp('ret to continue'); pause

disp('allF90laysJ is from 4:100')
[allF90laysJ(1) allF90laysJ(end)]





