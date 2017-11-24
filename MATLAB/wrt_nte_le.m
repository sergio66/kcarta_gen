% program wrt_nte
% same as /home/hannon/Fit_deltaR_nonLTE/Src/wrt_nte except it is for le machines
%
% Writes out non-LTE coefficients.  For use after running fitnonLTE
%
% Created: 15 March 2005, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of output binary fortran data file to create
outname = 'setnte_oct05.le.dat';
outname = 'nonLTE7_m150.le.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code below should not require modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[nchan, ncoef] = size(coefall);
[nchan, ncoef] = size(coef);

%% need to sort the channels so kCARTA does not have heart attacks
disp('sorting the channels for output ....')
[Y,I] = sort(freq);
freq   = freq(I);
idchan = idchan(I);
coef   = coef(I,:);

% Determine which channels to write out
%amaxc = max( abs(coefall),[],2 );
%iw = find( amaxc > 1E-7 );
iw = 1 : nchan;
nchanout = length(iw);

% Open output file
fid = fopen(outname,'w','ieee-le');

% FORTRAN record marker info
ifm = 4*(1 + 1 + ncoef); % 4 bytes each * (idchan, freq, coef)

% write out data to fortran file
for ii=1:nchanout
   ic = iw(ii);
   fwrite(fid,ifm,'integer*4');
   fwrite(fid,idchan(ic),'integer*4');
   fwrite(fid,freq(ic),'real*4');
   %fwrite(fid,coefall(ic,:),'real*4');
   fwrite(fid,coef(ic,:),'real*4');
   fwrite(fid,ifm,'integer*4');
end
fclose(fid);

disp(['finished writing data for ' int2str(nchanout) ' channels to file ' ...
   outname ])

%%% end of program %%%
