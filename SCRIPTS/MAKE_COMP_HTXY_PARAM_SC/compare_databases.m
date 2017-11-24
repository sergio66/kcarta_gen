function [] = compare_databases(glist,database1,database2);
%% eg compare_databases(1:40,'PARAM_HT2012/comp_ir605_2830.param','PARAM_HT2016/comp_ir605_2830.param')

data1 = load(database1);
data2 = load(database2);

if data1(1,4) ~= data2(1,4)
  disp('whoops did you load in like-for-liek database eg IR vs IR or did you do eg IR vs FIR')
end

for gg = 1 : length(glist)
  gid = glist(gg);
  gid1 = find(data1(:,1) == gid);
  gid2 = find(data2(:,1) == gid);
  clear freq1 freq2
  if length(gid1 >= 1) | length(gid2 >= 1)
    fprintf(1,'gasID = %2i \n',gid)
    freq1Start = data1(gid1,2);     freq1Stop = data1(gid1,3);
    freq2Start = data2(gid2,2);     freq2Stop = data2(gid2,3);
    figure(1); clf
    for ii = 1:length(freq1Start)
      hold on; line([freq1Start(ii) freq1Stop(ii)],[0.99 0.99],'color','b','linewidth',2);
      plot([freq1Start(ii) freq1Stop(ii)],[0.99 0.99],'bo')
    end
    for ii = 1:length(freq2Start)
      hold on; line([freq2Start(ii) freq2Stop(ii)],[1.01 1.01],'color','r','linewidth',2);
      plot([freq2Start(ii) freq2Stop(ii)],[1.01 1.01],'ro')      
    end
    title(num2str(gid))
    ax = axis; axis([ax(1) ax(2) 0.98 1.02]); grid
    hold off
    disp('ret to continue'); pause    
  end
end  