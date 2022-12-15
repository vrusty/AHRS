function X = ang_rot(X)

 X(3) = mod(X(3),2*pi);
% X(2) = 
% if sum(X(1:3)>=2*pi)
%     X(X(1:3)>=2*pi) = X(X(1:3)>=2*pi) - 2*pi;
%     
% elseif sum(X(1:3)<-2*pi)
%     X(X(1:3)<-2*pi) = X(X(1:3)<-2*pi) + 2*pi;
%     
% end

end