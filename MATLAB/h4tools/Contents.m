% 
%                   Matlab HDF4 Tools and Demos
%                   ---------------------------
% 
%                         Version 2.2
% 
%                         H. Motteler
%                          30 Oct 02
% 
% 
% USER INTERFACE 
% ---------------
% 
% rtpread.m, rtpwrite.m     -- read and write HDF4 RTP data as a
%                              structure of arrays
% 
% rtpread2.m, rtpwrite2.m   -- read and write HDF4 RTP data as an
%                              array of structures
% 
% sdload.m, sdsave.m        -- read and write HDF4 SD's and NetCDF
%                              files as a structure of arrays
% 
% h4sdread.m, h4sdwrite.m   -- read and write HDF4 SD's and NetCDF
%                              files as a cell list of arrays
% 
% h4vsread.m, h4vswrite.m   -- read and write HDF4 vdatas as a
%                              structure of arrays
% 
% h4sgread.m, h4sgwrite.m   -- read and write HDF4 vgroups of SD's as
%                              an array of structures
% 
% gpro2rtp.m, rtp2gpro.m    -- translate between RTP and GENLN2 format
%                              profiles
% 
% DEMOS and TESTS
% ----------------
% 
% rtptest4.m, rtptest5.m    -- basic demos of rtpread and rtpwrite
% 
% srfdemo.m                 -- how to use h4vsread.m and h4vswrite.m 
%                              to read and write AIRS SRFs
% 
% sdtest1.m                 -- tests h4sgread.m and h4sgwrite.m 
% sdtest2.m                 -- tests mat2sdsid.m and sdsid2mat.m
% sdtest3.m                 -- tests h4sdread.m and h4sdwrite.m
% vstest1.m                 -- tests mat2vsfid.m and vsfid2mat.m
% vstest2.m                 -- tests h4vsread.m and h4vswrite.m
% 
% 
% SELECTED UTILITIES
% -------------------
% 
% sructcmp.m                -- compare structures, allowing for different 
%                              field sets and values within some tolerance
% 
% stransp1.m, stransp2.m    -- translate between structure of arrays and
%                              an array of structures
% 
% gproread.m, gprowrite.m   -- read and write GENLN2 format user profiles
%
% mat2vsfid.m, vsfid2mat.m  -- write and read a matlab structure array
%                              as an HDF4 vdata, to an open HDF4 file ID
% 
% mat2sdsid.m, sdsid2mat.m  -- write and read a matlab array as an HDF4
%                              SDS, to an open HDF4 SD ID
% 

