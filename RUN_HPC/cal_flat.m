function [radLW radMW radSW frqLW frqMW frqSW] = cal_flat(vkc,rkc,nguard,opt1)

%{

Use kc2cris from the airs_decon git repo.  This is at github and
installed at /asl/packages/airs_decon.  kc2cris and other user
functions are in "source" there.  I use hamm_app and hamm_inv for
forward and inverse hamming apodization.  There are lots of examples
using these functions in the "test" directory in the same repo, see
cris_test5 for starters.

Similarly, to convolve kcarta to iasi use kc2iasi from the iasi_decon
repo.  /asl/matlib/fconv is somewhat obsolete because the functions
there were designed to keep the FFTs tractable for memory limits of
around 15 years ago.  The little parameter files there may still be
useful for instrument specs, and there is some other good stuff like
the old paired apodization and response functions

For an example of calling kc2cris for high res CrIS, all 3 bands and
with optional guard channels and tweaks to the passband filters, see
/asl/packages/ccast/motmsc/cal_flat.m

%}

%kcarta radiances with kc2cris

%addpath ../source
addpath /asl/s1/motteler/cris/ccast/motmsc/

if nargin == 2
  nguard = 4;

  opt1 = struct;
  opt1.inst_res = 'hires3';
  opt1.user_res = 'hires';     %% default 2235 hi res chan

   %% opt1.user_res = 'lowres';  (or 'midres' for the chirp stuff, if that comes up) : Howard suggestion Jan 2020
end

if (nargin == 3) | ~exist('opt1')
  opt1 = struct;
  opt1.inst_res = 'hires3';
  opt1.user_res = 'hires';     %% default 2235 hi res chan

   %% opt1.user_res = 'lowres';  (or 'midres' for the chirp stuff, if that comes up) : Howard suggestion Jan 2020
   %% opt1.user_res = 'midres';  (or 'midres' for the chirp stuff, if that comes up) : Howard suggestion Jan 2020
end

kcdir = '/asl/s1/motteler/kctest5/kcdata';

%opt1 = struct;
%opt1.inst_res = 'hires3';
%opt1.user_res = 'hires';

wlaser = 773.1307;
[instLW, userLW] = inst_params('LW', wlaser, opt1);
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[instSW, userSW] = inst_params('SW', wlaser, opt1);

%optLW = struct;  optLW.ng = 0;
%optMW = struct;  optMW.ng = 0;
%optSW = struct;  optSW.ng = 0;

% decided we want 4 guard chans
% "Note that cal_flat is just a wrapper for kc2cris.  You can modify it
% for any number of guard channels by setting the "ng" field in the
% options struct.  For example, for two guard channels for each band
% edge, set optLW.ng = 2; optMW.ng = 2, and optSW.ng = 2;"
%optLW = struct;  optLW.ng = 4;
%optMW = struct;  optMW.ng = 4;
%optSW = struct;  optSW.ng = 4;

optLW = struct;  optLW.ng = nguard;
optMW = struct;  optMW.ng = nguard;
optSW = struct;  optSW.ng = nguard;

% e5 filters
optLW.pL =  650; optLW.pH = 1100; optLW.rL = 15; optLW.rH = 20;
optMW.pL = 1200; optMW.pH = 1760; optMW.rL = 30; optMW.rH = 30;
optSW.pL = 2145; optSW.pH = 2560; optSW.rL = 30; optSW.rH = 30;

radLW = []; radMW = []; radSW = [];

for i = 1

%  kcmat = fullfile(kcdir, sprintf('kc%04d.mat', i));
%  d1 = load(kcmat);  
%  rkc = d1.radkc; vkc = d1.freqkc; clear d1

  [rtmp, frqLW] = kc2cris(userLW, rkc, vkc, optLW);
  radLW = [radLW, rtmp];

  [rtmp, frqMW] = kc2cris(userMW, rkc, vkc, optMW);
  radMW = [radMW, rtmp];

  [rtmp, frqSW] = kc2cris(userSW, rkc, vkc, optSW);
  radSW = [radSW, rtmp];

  if mod(i, 100) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n')

%save cal_flat radLW radMW radSW frqLW frqMW frqSW ...
%              userLW userMW userSW optLW optMW optSW
