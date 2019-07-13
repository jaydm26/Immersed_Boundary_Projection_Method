function [x,y] = ECL_inv(gamma)
    
    sf = L_inv_operation(gamma.x,"node");
     
    q = curl_2(sf);
    
    [x,y] = E_operation("edge",q);
end