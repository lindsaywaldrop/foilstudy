function test_adaptive_step()

lam = 0.05

for i=1:2000
    
    if mod(i,250)==0
        
        lam = 0.55*lam;
       
        fprintf('n=%d && lam=%.4f\n',i,lam);
        
    end
end