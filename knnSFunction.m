function knnSFunction(block)
    setup(block);
end

function setup(block)
    % Register number of input and output ports
    block.NumInputPorts = 1;
    block.NumOutputPorts = 1;

    % Setup functional port properties
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    % Set input port properties
    block.InputPort(1).Complexity = 'Real';
    block.InputPort(1).DataTypeId = 0;
   %  block.InputPort(1).Dimensions = [];
    block.InputPort(1).Dimensions = [1,69];
%block.InputPort(1).DimensionsMode = 'Variable';
    block.InputPort(1).DirectFeedthrough = true;

    % Set output port properties
    block.OutputPort(1).Complexity = 'Real';
    block.OutputPort(1).DataTypeId = 0;
    block.OutputPort(1).Dimensions = 1;

    % Register methods
    block.RegBlockMethod('Outputs', @outputs);    
%block.NumDialogPrms =1;
%block.SampleTimes = [0 0];
end

function outputs(block)
    input_data = block.InputPort(1).Data;
    label = predictKNN(input_data);
    block.OutputPort(1).Data = label;
end



