
function sds_id = mat2sdsid(sd_id, sname, A, attrs, dims)

%
% mat2sdsid -- write a matlab array as an SDS to an open HDF SD ID
%
% inputs
%
%   sd_id   - SD handle
%   sname   - SDS name for A
%   A       - numeric data array; 
%   attrs   - cell array of {name, value} attributes (optional)
%   dims    - cell array of dimension names (optional)
%
% output
%
%   sds_id  - open SDS handle
%
% notes
%
%   since the SDS is saved in row order, the HDF dimensions are flipped
%
% H. Motteler, 22 Nov 00
%

% default to empty attribute lists
if nargin < 4
  attrs = {};
end
if nargin < 5
  dims = {};
end

% set SDS "create" parameters
data_type = class(A);
dimsizes = size(A);
dimsizes = fliplr(dimsizes);
rank = length(dimsizes);

% create the SDS
data_type = htype(data_type);
sds_id = hdfsd('create', sd_id, sname, data_type, rank, dimsizes);
if sds_id == -1
  error('HDF SDS create failed');
end

% assign SDS dimension names
if length(dims) > rank
  error('dimension attribute assignments exceed array rank')
end
for i = 1 : length(dims)
  dim_num = rank - i;
  dim_id = hdfsd('getdimid', sds_id, dim_num);
  if dim_id == -1
    error('HDF SDS get dimension ID failed');
  end
  if length(dims{i}) == 0
    break
  end
  % status = hdfsd('setdimname', dim_id, [sname,'_',dims{i}]);
  status = hdfsd('setdimname', dim_id, dims{i});
  if status == -1
    sname, attrs, dims, rank, dimsizes, dim_id, i
    error('HDF SDS set dimension name failed');
  end
end

% assign SDS attributes 
for i = 1 : length(attrs)
  status = hdfsd('setattr', sds_id, attrs{i}{1}, attrs{i}{2});
  if status == -1
    error('HDF SDS set attribute failed');
  end
end

% set SDS write parameters
start = zeros(1,rank);
stride = ones(1,rank);
edges = dimsizes;

% write the SDS
status = hdfsd('writedata', sds_id, start, stride, edges, A);
if status == -1
  error('HDF SDS write failed');
end



