
% basic test of mat2sdsid and sdsid2mat

% echo on

sfile = 'sdtest2.hdf';

sd_id = hdfsd('start', sfile, 'create')
status = hdfsd('setattr', sd_id, 'file attr', 'file attr val')

A1 = rand(2,3,4);
name1 = 'Array A';
attrs1 = {{'units', 'tunas'}, {'scale', 'big'}};
dims1 = {'two', 'four'};
sds_id1 = mat2sdsid(sd_id, name1, A1, attrs1, dims1)

A2 = rand(3,2,5);
name2 = 'Array B';
attrs2 = {{'units', 'bogey'}, {'scale', 'small'}};
dims2 = {'three', 'two', 'five'};
sds_id2 = mat2sdsid(sd_id, name2, A2, attrs2, dims2)

% end access to the SDS's
status = hdfsd('endaccess',sds_id1)
status = hdfsd('endaccess',sds_id2)

% close the file
status = hdfsd('end',sd_id)


% re-open and read
sd_id = hdfsd('start', sfile, 'read')

% general info about file contents
[ndatasets, nglobal_attr, status] = hdfsd('fileinfo', sd_id)

for i = 0 : nglobal_attr - 1
  [aname,data_type,count,status] = hdfsd('attrinfo',sd_id,i)
  [aval,status] = hdfsd('readattr',sd_id,i)
end

for sds_index = 0 : ndatasets-1

  % get sds ID from index
  sds_id = hdfsd('select',sd_id, sds_index)

  [sdsname, rank, dimsizes, data_type, nattrs, status] = ...
				  	     hdfsd('getinfo', sds_id);

  [name, A, attrs, dims] = sdsid2mat(sd_id, sds_id);

  fprintf(1, 'array name = %s\n', name)
  for j = 1 : length(attrs);
    fprintf(1, '%s = %s\n', attrs{j}{1}, attrs{j}{2})
  end
  for j = 1 : length(dims);
    fprintf(1, 'd%d = %s\n', j, dims{j})
  end

  status = hdfsd('endaccess',sds_id)

end

% close the file
status = hdfsd('end',sd_id)
