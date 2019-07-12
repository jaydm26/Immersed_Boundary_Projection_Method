function [x,y] = E_operation(DataType,source)
    
    global body_map
    switch DataType
        case "node"
            x = zeros(length(body_map(:,1)),1);
            for k = 1:length(body_map(:,1))
                x(k,1) = E_op(body_map(k,1),body_map(k,2),source);
            end
        case "cell"
            x = zeros(length(body_map(:,1)),1);
            for k = 1:length(body_map(:,1))
                x(k,1) = E_op(body_map(k,1),body_map(k,2),source);
            end
        case "edge"
            x = zeros(length(body_map(:,1)),1);
            y = zeros(length(body_map(:,2)),1);
            for k = 1:length(body_map(:,1))
                x(k,1) = E_op(body_map(k,1),body_map(k,2),source,1);
                y(k,1) = E_op(body_map(k,1),body_map(k,2),source,2);
            end
    end
end
            
    