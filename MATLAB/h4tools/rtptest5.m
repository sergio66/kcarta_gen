
% test 5 -- translate GENLN2 user profiles to RTP and back to GENLN2
% 


% The fitting profiles, in GENLN input (level) format
gdir = '/home/motteler/asl/profiles/fitpro_july00/';

gpro = {'myp1', 'myp2', 'myp3', 'myp4', 'myp5', 'myp6'};

% hdf output filename
rtpfile = 'test5.rtp';

gpro2rtp(gdir, gpro, rtpfile);

rtp2gpro(rtpfile);

