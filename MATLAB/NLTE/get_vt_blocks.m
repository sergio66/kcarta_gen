function [raX] = get_blocks(iNumFullBlocks,iLeftover,fid);

%fprintf(1,'iNumFullBlocks,iLeftover = %3i %3i \n',iNumFullBlocks,iLeftover);

for ix = 1 : iNumFullBlocks
  strAll = fgetl(fid);
  xarr = strread(strAll); 
  mm = length(xarr);
  %fprintf(1,'%3i : %s : %3i \n',ix,strAll,mm);
  xarr = xarr(1:5);
  iaInd = (ix-1)*5+ (1:5);
  raX(iaInd) = xarr';
end

if iLeftover == 0
  return
end

strAll = fgetl(fid);
[x1] = strread(strAll); x1 = x1(1:end-1);
raX = [raX x1];
