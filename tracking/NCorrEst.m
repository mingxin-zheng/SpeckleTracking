function data = NCorrEst(tpl,tar,st,kn)
    % tpl: template of 2D rf data, matrix, size: # of Time Samples by # RF
    % Lines (e.g. 1870 x 128)
    
    % tar: target of 2D rf data, matrix, same size as template.
    
    % st: upper left coodinates (x,y) of the object in the template
    
    % kn: kernel which have settings, such as the search range (sx,sy), the
    % width of templates (rx, cy) and the extra samples to correlate
    % (hx,hy)
    
    % Update History
    % note: add point grid, remove usage of xSkip and ySkip
    
    x1 = st.x1; y1=st.y1;
    
    % kn: kernel
    % kn.rx & kn.cy specifies the width and height of the object's image in
    % the template
    
    rx = kn.rx; cy = kn.cy; 
    
    u = zeros(rx,cy);
    v = u;
    res = u;
    %update kn fields
    [locsx,locsy,kn] = locs(kn);
    
   
    hx = kn.hx; hy = kn.hy;
    
    % check if it is out of bound
    % kn.sx & kn.sy: specifies the search range of image between frames,
    % and the displacement of the object should not be larger than [sx,sy]
    sx1 = min([kn.sx, x1+locsx(1)-hx-1]);
    sx1 = max(sx1,0);
    
    sx2 = min([kn.sx, size(tpl,1)-(x1+locsx(end))-hx]);
    sx2 = max(sx2,0);
    
    sy1 = min([kn.sy, y1+locsy(1)-hy-1]);
    sy1 = max(sy1,0);
    
    sy2 = min([kn.sy, size(tpl,2)-(y1+locsy(end))-hy]);
    sy2 = max(sy2,0);
    
    
    mask = accelerate(tpl,kn);
    
    for k = 1:rx
       for p = 1:cy
           ix = locsx(k);
           iy = locsy(p);
           val = (x1+ix)>0 & (x1+ix)<=size(mask,1) & (y1+iy)>0 & (y1+iy)<=size(mask,2);
           if val&&(mask(x1+ix,y1+iy)|~any([rx-1,cy-1])) == 1

           ref    = tpl(x1+ix-hx   :x1+ix+hx   , ...
                        y1+iy-hy   :y1+iy+hy);
           target = tar(x1+ix-hx-sx1:x1+ix+hx+sx2, ...
                        y1+iy-hy-sy1:y1+iy+hy+sy2);
           c = normxcorr2(ref,target);
           [res(k,p),ind] = max(c(:));
           [mr,~] = size(c);
           
           % suppress warning
           
           v_subest = 0;
           if (ind>mr) && (ind+mr)<numel(c)
               y = [c(ind-mr) c(ind) c(ind+mr)];
               v_subest = QuadFit(y);
           end
           u(k,p) = mod(ind,mr)-(2*hx+1)-sx1;
           v(k,p) = ceil(ind/mr)-(2*hy+1)-sy1+v_subest;
           
           end
           
       end
    end
    
    data.u = u;
    data.v = v;
    data.res = res;
end

function mask = accelerate(rf,kn)
    % reject data patch that has low signal intensity
    env = abs(hilbert(double(rf)));
    win= ones(kn.hx*2+1,kn.hy*2+1);
    env_fil = imfilter(env,win/sum(win(:)));
    mask = env_fil>0; % was 150 
    mask(:,1:kn.hy) = 0;
    mask(:,end-kn.hy:end) = 0;
    mask(1:kn.hx,:) = 0;
    mask(end-kn.hx:end,:) = 0;
end

function [x,y,kn] = locs(kn)
    rx = kn.rx;
    cy = kn.cy;
    w = kn.w;
    h = kn.h;
    
    if rx>h || cy>w
        error('number of elements in rx/cy needs to be larger than h/w');
    end
    
    if rx ~= round(rx) || cy ~=round(cy)
        error('kn containts number that is not integer');
    end
    
    if rx == 1
            x = round(h/2);
            kn.hx = max(kn.hx,floor(kn.h/2));
    else
            x = round(linspace(0,h-1,rx));
            kn.hx = max(kn.hx,round(mean(diff(x))/2));
    end
    
    if cy == 1
            y = round(w/2);
            kn.hy = max(kn.hy,floor(kn.w/2));
    else
            y = round(linspace(0,w-1,cy));
            kn.hy = max(kn.hy,round(mean(diff(y))/2));
    end
    
end

function [ output ] = QuadFit( y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
	temp1 = y(1) - y(3);
	temp2 = y(1) + y(3) - 2*y(2);
    if  abs(temp1)<1e-6
		output = 0;
    else
		answer = temp1/(temp2*2.0);
		
        if answer<1
            output = answer;
        else
            output = 0.0;
        end
    end

end

