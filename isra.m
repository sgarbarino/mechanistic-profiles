function x = isra(A, At, y, x, opt)
    
    iter = 0;
    Aty = At*y;
%     x_all = zeros(opt.maxiter,numel(x));
    
    while(iter < opt.maxiter)  
%           xold = x;
        x = x .* (Aty) ./ (At*(A*x));
        iter = iter + 1;
        if mod(iter,100) == 0
%             disp(iter);
%             disp(x)
        end
    end
end