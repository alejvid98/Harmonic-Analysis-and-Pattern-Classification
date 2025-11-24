function label = predictKNN(input_data)
    mdlStruct = load('knnModel.mat', 'mdl');
    mdl = mdlStruct.mdl;
    
% data = load('nuevos_datos/arritmia_prueba5.mat'); %1
%data =  load('nuevos_datos/sana_prueba1.mat'); %0
%data =  load('nuevos_datos/inf_prueba1.mat'); %2
%data =  load('nuevos_datos/insuf_prueba2.mat'); %3



input_data_numeric = input_data;

 %   input_data_numeric = data.P1; 

 input_data_numeric = input_data_numeric(:,1:63);    
 label = predict(mdl, input_data_numeric);

end

