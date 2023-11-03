




load('DATA1');load('DATA2');load('DATA3');
data=[DATA1;DATA2];
labels=[1:51,1:51];
t_data=DATA3;
t_labels=1:51;

[predicted_labels,nn_index,accuracy,ed,ind]=KNN_(3,data,labels,t_data,t_labels);
predicted_labels
nn_index
accuracy