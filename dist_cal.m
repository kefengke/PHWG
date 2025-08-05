function L = dist_cal(a,b,W)
    %assert(size(a,2) ==2 && size(b,2) ==2);
    t1 = abs(a(:,1)-b(:,1));
    t2 = abs(a(:,2)-b(:,2));
    L = sqrt((min(t1,2*W - t1)).^2 + (min(t2,2*W - t2)).^2);
    %assert(max(L)<=sqrt(2)*W);
end

