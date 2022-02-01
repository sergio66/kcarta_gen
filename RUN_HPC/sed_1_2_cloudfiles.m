pjunk = hdfread(use_this_rtp0,'profiles','Fields',{'cfrac','cfrac2','cfrac12','ctype','ctype2','cngwat','cngwat2','cprtop','cprtop2','cpsize','cpsize2'});
cfrac   = pjunk{1}(iiBin);
cfrac2  = pjunk{2}(iiBin);  
cfrac12 = pjunk{3}(iiBin);
ctype   = pjunk{4}(iiBin);
ctype2  = pjunk{5}(iiBin);
cngwat   = pjunk{6}(iiBin);
cngwat2  = pjunk{7}(iiBin);
cprtop   = pjunk{8}(iiBin);
cprtop2  = pjunk{9}(iiBin);
cpsize   = pjunk{10}(iiBin);
cpsize2  = pjunk{11}(iiBin);
clear pjunk

fprintf(1,' cloud1 : ctype,cfrac,cngwat,cprtop,cpsize = %2i %8.6f %8.6e %8.6f %8.6f \n',ctype,cfrac,cngwat,cprtop,cpsize);
fprintf(1,' cloud2 : ctype,cfrac,cngwat,cprtop,cpsize = %2i %8.6f %8.6e %8.6f %8.6f \n',ctype2,cfrac2,cngwat2,cprtop2,cpsize2);

%{
if (ctype > 100 & ctype2 > 100)
  iNclouds = 2;
  iCloudType1 = ctype;
  if ctype == 101
    strCloudType1 = strWaterCloud
  elseif ctype == 102
    strCloudType1 = strIceCloud;
  end
  iCloudType2 = ctype2;
  if ctype2 == 101
    strCloudType2 = strWaterCloud
  else
    strCloudType2 = strIceCloud;
  end
  sedder = [sedder ' -e "s/iCloudType1/'     num2str(ctype) '/g"  -e "s/strCloudType1/'  strCloudType1 '/g"'];
  sedder = [sedder ' -e "s/iCloudType2/'    num2str(ctype2) '/g"  -e "s/strCloudType2/'  strCloudType2 '/g"'];    
  sedder = [sedder ' ' type2nml '  > ' outnml];
 
elseif (ctype < 100 & ctype2 < 100)
  iNclouds = 0;
  sedder = [sedder ' template_Qrad.nml  > ' outnml];        
else
  iNclouds = 1;
  if ctype == 101
    strCloudType1 = strWaterCloud;
  else
    strCloudType1 = strIceCloud;
  end
  sedder = [sedder ' -e "s/iCloudType1/'     num2str(ctype) '/g"  -e "s/strCloudType1/'  strCloudType1 '/g"'];
  sedder = [sedder ' ' type1nml '  > ' outnml];    
end
%}

%% kCARTA automatically reads in ctype1, ctype2 and then figures out which cloudtype to assign to what
iNclouds = 2;
iCloudType1 = ctype;
iCloudType1 = 101;
if iCloudType1 == 101
  strCloudType1 = strWaterCloud;
elseif iCloudType1 == 201
  strCloudType1 = strIceCloud;
end
iCloudType2 = ctype2;
iCloudType2 = 201;
if iCloudType2 == 101
  strCloudType2 = strWaterCloud;
elseif iCloudType2 == 201    
  strCloudType2 = strIceCloud;
end
sedder = [sedder ' -e "s/iCloudType1/'     num2str(ctype) '/g"  -e "s/strCloudType1/'  strCloudType1 '/g"'];
sedder = [sedder ' -e "s/iCloudType2/'    num2str(ctype2) '/g"  -e "s/strCloudType2/'  strCloudType2 '/g"'];    
sedder = [sedder ' ' type2nml '  > ' outnml];
