function gamma = CTH(Fx,Fy)
    
    q = H_operation("edge",Fx,Fy);
    
    gamma = curl_2(q);
end