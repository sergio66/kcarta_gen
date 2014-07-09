
function [name, A, attrs, dims] = sdsid2mat(sd_id, sds_id);

%
% sdsid2mat -- read a matlab array from an open HDF SDS ID
%
% INPUTS
%   sd_id   - SD handle
%   sds_id  - SDS handle
%
% OUTPUT
%   name    - SDS name for A
%   A       - numeric data array; 
%   attrs   - n x 2 cell array of {name value} pairs
%   dims    - d x 1 cell array of dimension names
%
% H. Motteler, 22 Nov 00
%

% get msc info on this SDS, by ID 
[name, rank, dimsizes, data_type, nattrs, status] = ...
				  	     hdfsd('getinfo', sds_id);
if status == -1
  error('HDF SDS get info failed')
end

% read data
start = zeros(1, rank);
stride = ones(1, rank);
edge = dimsizes;
[A, status] = hdfsd('readdata', sds_id, start, stride, edge);
if status == -1
  error('HDF SDS read data failed')
end

% read attributes
attrs = {};
for attr_index = 0 : nattrs-1
  [aname, data_type, count, status] = hdfsd('attrinfo', sds_id, attr_index);
  if status == -1
    error('HDF SDS get attribute info failed')
  end
  [aval, status] = hdfsd('readattr', sds_id, attr_index);
  attrs{attr_index + 1} = {aname, aval};
  if status == -1
    error('HDF SDS get attribute value failed')
  end
end

% read dimension names
dims = {};
for dim_num = 0 : rank-1
  dim_id = hdfsd('getdimid', sds_id, dim_num);
  if dim_id == -1
    error('HDF SDS get dimension ID failed')
  end
  [dname,count,data_type,nattrs,status] = hdfsd('diminfo',dim_id);
  if status == -1
    error('HDF SDS get dimension info failed')
  end
  % dims{rank - dim_num} = dname(length(name)+2:length(dname));
  dims{rank - dim_num} = dname;
end

