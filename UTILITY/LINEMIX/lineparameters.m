function [line] = lineparameters(idgas,f1,f2,iso,iL,iU,fnameout,fband);

%%function [line] = lineparameters(idgas,f1,f2,iso,iL,iU,fnameout,fband);
% this file takes data from HITRAN and writes it out
% to a f77 binary file, so that kCARTA can compute some Doppler lineshapes!
% at non LTE temperatures
% input : idgas : gasID
%         f1,f2 : start and stop frequencies within which to search HITRAN
%         iso   : gas isotope
%         iL,iU : global lower and upper quantum states
%         fnameout : string for output name
%         band     : identifier for output name
%
% The header info contains the following integers on one line : idgas,npts. 
%    These are the gasID, number of spectra points
% Then save line parameters : elower,freq,j_lower,p_shift,stren,w,w_s,w_temp

clf
path(path,'/home/sergio/SPECTRA');

%idgas = input('Enter gasID : ');
%f1    = input('Enter start freq for hitread : ');
%f2    = input('Enter stop freq for hitread : ');

fname  = ['/asl/data/hitran/h2k.oldiso/g' num2str(idgas) '.dat'];
fprintf(1,'reading in HITRAN data %s \n',fname);
[line] = hitread(f1,f2,0,idgas,fname)

%iso   = input('Enter gas isotope : ');
%iL    = input('Enter gas lower quantum number  : ');
%iU    = input('Enter gas upper quantum number  : ');

ii=find((line.ilsgq==iL) & (line.iusgq==iU) & (line.iso==iso)); 
len = length(ii);
if (len <= 0)
  error('not enuff lines here!!!!')
  end

semilogy(line.wnum(ii),line.stren(ii),'.'); grid 

%fnameout = input('Enter output prefix (.dat) filename : ');
%fband    = input('Enter identifying band  : ');
fname = [fnameout '_' num2str(fband) '.dat']
fid=fopen(fname,'w','ieee-le');

% Write header info
% header1 = gasid, npts, isotope

stren   = line.stren(ii);
elower  = line.els(ii);
freq    = line.wnum(ii);
j_lower = line.ilsgq(ii);
j_upper = line.iusgq(ii);
p_shift = line.tsp(ii);
w       = line.abroad(ii);
w_s     = line.sbroad(ii);
w_temp  = line.abcoef(ii);

npts = length(stren);
isotope = iso;
dd= diff(iso); dd=sum(dd);
if (abs(dd) > eps)
  error('differing isotopes found here!!!')
  end

filemark= 4 + 4 + 4;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[idgas,npts,isotope],'integer*4');
fwrite(fid,filemark,'integer*4');

%Save the needed vectors 
filemark= 8 * npts;

fwrite(fid,filemark,'integer*4');
fwrite(fid,elower,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,freq,'real*8');
fwrite(fid,filemark,'integer*4');

lowerj = j_lower * 1.0;
fwrite(fid,filemark,'integer*4');
fwrite(fid,lowerj,'real*8');
fwrite(fid,filemark,'integer*4');

upperj = j_upper * 1.0;
fwrite(fid,filemark,'integer*4');
fwrite(fid,upperj,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,p_shift,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,stren,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w_s,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w_temp,'real*8');
fwrite(fid,filemark,'integer*4');

fclose(fid);
