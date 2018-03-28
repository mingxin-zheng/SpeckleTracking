function img = b_map(rf,option)

label_option = true;

if isstruct(rf)

    fs = rf.zone.transmitFrequency(1);
    if isfield(rf,'d_x')
        d_x = rf.d_x;
    else
        if isfield(rf,'resolution')
            d_x = rf.resolution(1);
            fs = 1540/(rf.resolution(2)*2);
        else
            d_x = 0.3/1000; 
        end
    end
    rf_lines = rf.rf;
else

    fs = 36e6;
    d_x = 0.3/1000;
    rf_lines = rf;
end

if ~isa(rf_lines,'double')
    rf_lines = double(rf_lines);
end
env=abs(hilbert(rf_lines));    % take Hilbert transfrom of real signal
dispFigure = 1;
switch nargin
    case 1   
        new_env=2*round((env.^0.5)/(1e6^0.4)*250);
    case 2
        switch option  
            case 'log'
                log_env=env(1:end,:)/max(env(:));
                log_env=20*log10(log_env);
                new_env=127/60*(log_env+60);
            case 'no_show' %3
                new_env=2*round((env.^0.5)/(1e6^0.4)*250);
                dispFigure = 0;
            case 'no_show_logmap'
                log_env=env(1:end,:)/max(env(:));
                log_env=20*log10(log_env);
                new_env=127/60*(log_env+60);
                dispFigure = 0;
            case 'no_label' % create b_map from rf line data w/o labels
                new_env=2*round((env.^0.5)/(1e6^0.4)*250);
                label_option = false;
            case 'envelop' % create b_map from abs(hilbert(rf)) data.
                new_env = 2*round((rf.^0.5)/(1e6^0.4)*250);
                label_option = true;
            case 'envelop_no_label'
                new_env = 2*round((rf.^0.5)/(1e6^0.4)*250);
                label_option = false;
        end

    otherwise
        new_env=2*round((env.^0.5)/(1e6^0.4)*250);
end
    if dispFigure
        [n,m]=size(new_env);

        image(((1:m-1)*d_x-m*d_x/2)*1000,(1:n)/fs*1540/2*1000,new_env)
        axis('image')
        
        if label_option
        
        xlabel('Lateral distance [mm]', 'FontSize', 18)
        

        % Display Element with axial depth
        % image(1:m,(1:n)/fs*1540/2*1000,new_env)
        % xlabel('Lateral Element Number [mm]')
        % axis('normal')
        ylabel('Axial distance [mm]', 'FontSize', 18)
        else
            axis off;
        end
        colormap(gray(127))
    end
    
    img = new_env;
    grid(gca,'minor')
    set(gca,'FontSize',18);
end



%the multiplier 2 is added 5/3/2016
