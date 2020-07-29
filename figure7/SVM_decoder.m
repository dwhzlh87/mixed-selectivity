function [x_time,time_perform,crv_cell,population_response]=SVM_decoder(win_size,step_size,temp_index, mode)
%this version works with conjunction data
% win_size and step_size to compute slide windown spike count,temp_index
% select cells to use
%mode=1;    %for testing
%temp_index=find_informative; 
%step_size=0.1;
%win_size=0.4;
load('new_conj_data.mat');
informative_select=temp_index;
edge1=[-1:step_size:5-win_size];   %6 second of time 
edge2=[-1+win_size:step_size:5];
cell_id=[];
cell_cueloc=[];
cell_cuefeature=[];
cell_sampleloc=[];
cell_samplefeature=[];
cell_ifmatch=[];
cell_locmatch=[];
cell_featurematch=[];

cell_bincount=[];
trial_count=0;
class_cueloc=[ones(1,8),2*ones(1,8),ones(1,8),2*ones(1,8)];
class_cuefeature=[ones(1,16),2*ones(1,16)];
class_sampleloc=repmat([1,1,2,2],1,8);
class_samplefeature=repmat([1,1,1,1,2,2,2,2],1,4);
class_locmatch=double(class_cueloc==class_sampleloc);
class_featurematch=double(class_cuefeature==class_samplefeature);
class_allmatch=double(class_locmatch & class_featurematch);

cueloc_lookup=[1,2,5,6];
cuefeature_lookup=[1,2,3,4];
match_lookup=[1,3,6,8];
sampleloc_table=[2,1];
samplefeature_table=[2,1];
%data extraction
for i=1:length(informative_select) %loop cell

    temp_cell_data=new_conj_cell(informative_select(i),:);    
    for j=1:8  %loop cuecalss
        temp_cueclass_data=temp_cell_data{j};
        cell_id=[cell_id,i*ones(1,length(temp_cueclass_data))];
        cueclass_trialclass=[temp_cueclass_data.trialclass];
        
        temp_cueloc=class_cueloc(cueclass_trialclass);
        cell_cueloc=[cell_cueloc,temp_cueloc]; 
        
        temp_cuefeature=class_cuefeature(cueclass_trialclass);
        cell_cuefeature=[cell_cuefeature,temp_cuefeature];
        
        temp_ismatch=[temp_cueclass_data.IsMatch];
        cell_ifmatch=[cell_ifmatch,temp_ismatch];
        
        temp_samplefeature=class_samplefeature(cueclass_trialclass);
        cell_samplefeature=[cell_samplefeature,temp_samplefeature];
        
        temp_sampleloc=class_sampleloc(cueclass_trialclass);
        cell_sampleloc=[cell_sampleloc,temp_sampleloc];
        
        temp_locmatch=class_locmatch(cueclass_trialclass);
        cell_locmatch=[cell_locmatch,temp_locmatch];
        
        temp_featurematch=class_featurematch(cueclass_trialclass);
        cell_featurematch=[cell_featurematch,temp_featurematch];
        
        cueclass_ontime=[temp_cueclass_data.Cue_onT];
        cueclass_spiketimes={temp_cueclass_data.TS};
        for p=1:length(temp_cueclass_data) %loop through trial
            trial_count=trial_count+1;
            temp_TS=cueclass_spiketimes{p}-cueclass_ontime(p);
            for q=1:length(edge1)  %loop through timebin            
            cell_bincount(trial_count,q)=length(find(temp_TS>edge1(q) & temp_TS<edge2(q)));
            end
        end
    end
end
match_label=[];
locmatch_label=[];
featurematch_label=[];
cueloc_label=[];
sampleloc_label=[];
cuefeature_label=[];
samplefeature_label=[];
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


%build pseudopopulation
for i=1:length(informative_select) %loop through cell
  
    temp_cellindex=find(cell_id==i);
    ind_cell_bincount=cell_bincount(temp_cellindex,:);
    ind_cell_cueloc=cell_cueloc(temp_cellindex);
    ind_cell_cuefeature=cell_cuefeature(temp_cellindex);
    ind_cell_locmatch=cell_locmatch(temp_cellindex);
    ind_cell_featurematch=cell_featurematch(temp_cellindex);
  %  ind_cell_ifmatch=cell_ifmatch(temp_cellindex);
    %loop through cueclasses
    lookup_table1=[1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2]; %cueloc
    lookup_table2=[1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]; %cuefeature
    lookup_table3=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]; %locmatch
    lookup_table4=[1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1]; %featurematch
    for j=1:16
    temp_index=find((ind_cell_cueloc==lookup_table1(j) & ind_cell_cuefeature==lookup_table2(j)) & (ind_cell_locmatch==lookup_table3(j) & ind_cell_featurematch==lookup_table4(j)));
    rand_select=randperm(length(temp_index));
    select_index=rand_select(1:6);
    population_response(i,1+(j-1)*6:6*j,:)=ind_cell_bincount(temp_index(select_index),:);
    end

end
t=templateSVM('KernelFunction','linear');
%svm for each time point
for i=1:length(edge1)
   % disp(i);
    if mode==1
        decode_label=cueloc_label;
    elseif mode==2
        decode_label=cuefeature_label;
    elseif mode==3
        decode_label=sampleloc_label;
    elseif mode==4
        decode_label=samplefeature_label;
    elseif mode==5
        decode_label=locmatch_label;
    elseif mode==6
        decode_label=featurematch_label;
    elseif mode==7
        decode_label=match_label;
    end
    OBJ=fitcecoc(population_response(:,:,i)',decode_label,'Coding','onevsall','Learners',t);
    cvecoc = crossval(OBJ,'KFold',10);
   % cvecoc = crossval(OBJ,'leaveout','on');
    Yhat = kfoldPredict(cvecoc);
    confuse_mat=confusionmat(decode_label,Yhat);
    num_correct=diag(confuse_mat);
    correct_pertrue_class=num_correct./sum(confuse_mat,2);
    mean_correct_linearSVM=mean(correct_pertrue_class);
    time_perform(i)=mean_correct_linearSVM;
    crv_cell{i}=cvecoc;
end
x_time=linspace(-1,5,length(edge1))+0.5*win_size;
end
