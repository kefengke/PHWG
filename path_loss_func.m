function L = path_loss_func(x)

assert(isempty(find(x<0)))
alpha=4;

L = min(1,x.^(-alpha));
% L = zeros(length(x),1);d=0.5;
% for ii = 1:length(x)
%     if x(ii)<= d
%         L(ii) = 100;
%     end
% %L = (1+x).^(-alpha);
 end

