this subdir has a bunch of files that allows the user to read in
kcarta output (standard and jacobian) without using fseeks!!!!!!!

the output from these readers can be easily read in, as the data files
saved here are in the form 

  num1  num2  M   
  raCFreq        1xN
  raaCData       MxN

where N = num2-num1+1 is the number of frequency points and 
      M = number of paths/mixed paths/radiances, or jacobians, that
          are output

the array of data is saved in the format
raCFreq
for ii=1:M
  raData
  end

where raCFreq, raData are arrays of data, length N

this makes the readers slower, but should work on all platforms.

sergio 1/27/99

