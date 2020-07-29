%this version works for conjunction data
%load('conj_informative.mat');
mode=7; %3 sampleloc?4 samplefeature, 5 locmatch, 6 featurematch, 7 allmatch
population=2; %1 means CS+LMS %2 means NMS
load('conj_sampledelay_CS.mat');
load('conj_sample_CS.mat');
if mode<5
all_CS=[save_sigsample_CS{1};save_sigsample_CS{2};save_sigsample_CS{3}];
else
all_CS=[save_sigsampledelay_CS{1};save_sigsampledelay_CS{2};save_sigsampledelay_CS{3}];
end

load('conj_sampledelay_LMS.mat');
load('conj_sample_LMS.mat');
if mode<5
all_LMS=[save_sigsample_LMS{1};save_sigsample_LMS{2};save_sigsample_LMS{3};save_sigsample_LMS{4}];
else
all_LMS=[save_sigsampledelay_LMS{1};save_sigsampledelay_LMS{2};save_sigsampledelay_LMS{3};save_sigsampledelay_LMS{4}];
end

load('conj_sampledelay_NMS.mat');
load('conj_sample_NMS.mat');
if mode<5
all_NMS=[find_sigsample_NMS];
else
all_NMS=[find_sigsampledelay_NMS];
end

if population==1
find_informative=unique([all_CS;all_LMS]);
else
find_informative=unique([all_NMS]);
end
 
for r=1:5
    disp(r);
[x,y,crv_cell,population_response]=SVM_decoder(0.4,0.1,find_informative,mode);
accuracy_mat=zeros(length(x),length(x));

match_label=[];
locmatch_label=[];
featurematch_label=[];
cueloc_label=[];
sampleloc_label=[];
cuefeature_label=[];
samplefeature_label=[];
cueloc_lookup=[1,2,5,6];
cuefeature_lookup=[1,2,3,4];
match_lookup=[1,3,6,8];
sampleloc_table=[2,1];
samplefeature_table=[2,1];
for i=1:8
   if ismember(i,cueloc_lookup)   
    cueloc_label=[cueloc_label,ones(1,12)];
   else
    cueloc_label=[cueloc_label,2*ones(1,12)];
   end
   
   if ismember(i,cuefeature_lookup)   
    cuefeature_label=[cuefeature_label,ones(1,12)];
   else
    cuefeature_label=[cuefeature_label,2*ones(1,12)];
   end
   
   if ismember(i,match_lookup)   
    match_label=[match_label,[1,1,1,1,1,1,0,0,0,0,0,0]];
   else
    match_label=[match_label,[0,0,0,0,0,0,0,0,0,0,0,0]];
   end
   
   if ismember(i,match_lookup)   
    featurematch_label=[featurematch_label,[1,1,1,1,1,1,1,1,1,1,1,1]];
   else
    featurematch_label=[featurematch_label,[0,0,0,0,0,0,0,0,0,0,0,0]];
   end
    
   if ismember(i,match_lookup)   
    locmatch_label=[locmatch_label,[1,1,1,1,1,1,0,0,0,0,0,0]];
   else
    locmatch_label=[locmatch_label,[1,1,1,1,1,1,0,0,0,0,0,0]];
   end
   
   sampleloc_label=cueloc_label.*locmatch_label+sampleloc_table(cueloc_label).*(1-locmatch_label);
   samplefeature_label=cuefeature_label.*featurematch_label+samplefeature_table(cuefeature_label).*(1-featurematch_label);

end

for i=1:length(y) %loop through training time point
  
    temp_crv=crv_cell{i};
    temp_partition=temp_crv.Partition;
    for j=1:length(y)  %loop through testing time point 
     %  if i==j
     %      accuracy_mat(i,j)=y(i);
     %  else
           for p=1:10  % 10 fold cross validation
               fold_decoder=temp_crv.Trained{p};
               fold_test_index=find(training(temp_partition,p)==0);
               label=predict(temp_crv.Trained{p},population_response(:,fold_test_index,j)');
               if mode==1
               correct_count(p)=length(find(label==cueloc_label(fold_test_index)')); %change what to decode as needed
               elseif mode==2
               correct_count(p)=length(find(label==cuefeature_label(fold_test_index)')); %change what to decode as needed
               elseif mode==3
               correct_count(p)=length(find(label==sampleloc_label(fold_test_index)')); %change what to decode as needed
               elseif mode==4
               correct_count(p)=length(find(label==samplefeature_label(fold_test_index)')); 
               elseif mode==5
               correct_count(p)=length(find(label==locmatch_label(fold_test_index)'));
               elseif mode==6
               correct_count(p)=length(find(label==featurematch_label(fold_test_index)'));
               else 
               correct_count(p)=length(find(label==match_label(fold_test_index)'));
               end
               all_count(p)=length(fold_test_index);
           end
           accuracy_mat(i,j)=sum(correct_count)/sum(all_count);
    %   end
    end
end
store_mat(r,:,:)=accuracy_mat;
end
plot_mat=squeeze(mean(store_mat,1));
save sampledelayallmatch_nonlinear.mat store_mat
[X,Y] = meshgrid(-0.8:0.1:4.8);
[Xq,Yq] = meshgrid(-0.8:0.02:4.8);
Vq = interp2(X,Y,plot_mat,Xq,Yq);
imagesc(flipud(Vq));
colormap jet;
hold on;
plot([40,40],[0,281],'k','LineWidth',2);
plot([65,65],[0,281],'k','LineWidth',2);
plot([140,140],[0,281],'k','LineWidth',2);
plot([165,165],[0,281],'k','LineWidth',2);
plot([0,281],[240,240],'k','LineWidth',2);
plot([0,281],[215,215],'k','LineWidth',2);
plot([0,281],[140,140],'k','LineWidth',2);
plot([0,281],[115,115],'k','LineWidth',2);
yticks([40 90 140 190 240]);
yticklabels({'4s','3s','2s','1s','0s'});
xticks([40 90 140 190 240]);
xticklabels({'0s','1s','2s','3s','4s'});
%ylabel('training time');
%xlabel('testing time');
caxis([0.5,1]);