function Y = wrap_points_torus(X,win)
% function: wrap points within [-2W,2W]^2 to [-W,W]^2
    Y=X;
    assert(size(X,2) ==2);
    for ii = 1:size(X,1)
        assert(abs(X(ii,1))<=2*win && abs(X(ii,2))<=2*win);
        vec = [X(ii,1)-2*win,X(ii,1)+2*win,X(ii,1)];
        [~,ind]=min(abs(vec));
        Y(ii,1) = vec(ind);
        vec = [X(ii,2)-2*win,X(ii,2)+2*win,X(ii,2)];
        [~,ind]=min(abs(vec));
        Y(ii,2) = vec(ind);
        assert(abs(Y(ii,1))<=win && abs(Y(ii,2))<=win);
    end
end
