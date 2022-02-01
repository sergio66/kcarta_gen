function iCnt = list_anomaly_files_to_be_made(fid,notfinish)

iCnt = 0;

ii = 1;
iP = 0;
iStart = [];
iContinuous = -1;

while ii <= length(notfinish)
  if iP >= 0
    iCnt = iCnt + 1;
    fprintf(fid,'%s',num2str(notfinish(ii)));
    iStart = notfinish(ii);
  elseif iP == -1
    fprintf(fid,'-');
  end

  if length(notfinish) > ii
    if notfinish(ii)+1 == notfinish(ii+1)
      iP = iP - 1;
      iContinuous = +1;
    else
      iP = 0;
      if iStart ~= notfinish(ii)
        fprintf(fid,'%s,',num2str(notfinish(ii)));
      else
        fprintf(fid,',');
      end
    end
  elseif length(notfinish) == ii
    if iStart ~= notfinish(ii)
      fprintf(fid,'%s',num2str(notfinish(ii)));
    end
  end
  ii = ii + 1;
end
