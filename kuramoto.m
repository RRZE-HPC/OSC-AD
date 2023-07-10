% We focus on closing the gap and come-up with analytical modelling by analyzing 
% at the mesoscopic level what parameters shape (de-re)synchronization of the 
% parallel programs based on a system of Kuramoto-type coupled phase oscillators.  

function Hout = kuramoto(n,kappa,beta)
    if nargin == 1
        preset = n;
    else
        preset = 0;
    end
    if nargin == 0
        n = 5;
        kappa = rand;
        beta = .660*kappa;
    elseif nargin < 3
        n = [];       % Number of oscillators. real scalars
        kappa = [];   % Coupling coefficient. real scalars
        beta = [];    % Half-width of omega interval.
    else % nargin == 3
        % n,kappa,beta are input
    end
    
    width = [];   % Standard deviation of display ring.
    step = [];    % Display step size.
    pot = [];     % Graph potential or order parameter.
    uni = [];     % Uniform or normal distribution of omega.
    quit = [];    % 0 <= t <= quit.
    H = [];       % Optional output
    wanth = nargout > 0;
    delta = [];
    
    init_vars(preset)
    
    % omega and theta are real column vectors of length n
    point = []; flag = [];  animal = []; titl = [];
    theta = []; omega = []; colors = []; r = [];
    rotate = []; txt = []; ax2 = []; psi = [];

    % topology_ij the elements of the adjacenct/connective topology matrix of size NxN
    % topology_ij = 1 if there is a connection between ossilator i and j, or topology_ij = 0, otherwise
    % periodicity/non-periodicity can make this symmetric / asymmetric.

       topology = zeros(n);
       
       % all-to-all
       % topology = ones(n);
        
       % 1. Direct short distance communication
       % uni--directional non-periodic          
%         topology(2:n+1:end) = 1; % lower first subdiaginal 
       
        % no propogation as sudden hit by boundary
%         d=1;
%         topology(d*n+1:n+1:end) = 1; % upper first subdiagonal
               
        %uni--directional periodic               
%         d=1;
%         topology(d*n+1:n+1:end) = 1; % upper first and lower (n-1)th subdiagonal
%         topology(n,1) = 1;

        % bi-directional non-periodic % both upper and lower 1st subdiagonal
%         topology(2:n+1:end) = 1;
%         d=1;
%         topology(d*n+1:n+1:end) = 1;

        % bi-directional periodic   % both upper and lower 1st and (n-1)th subdiagonals   
%         topology(2:n+1:end) = 1;
%         d=1;
%         topology(d*n+1:n+1:end) = 1; 
%         d=n-1;
%         topology(d*n+1:n+1:end) = 1; 
%         topology(n,1) = 1;

       %2. Indirect long distance communication
       % bi-directional non-periodic next-to-next ossilator only (d=+-2) % both upper and lower 2nd subdiagonal
%         d=2;
%         topology(d*n+1:n+1:end) = 1;
%         topology(3:n+1:end) = 1;
%         topology(1,n) = 0;

       % bi-directional periodic next-to-next ossilator only (d=+-2) % both upper and lower 2nd and (n-2)th subdiagonals
%         d=2;
%         topology(d*n+1:n+1:end) = 1;
%         topology(3:n+1:end) = 1;
%         topology(1,n) = 0;
%         d=n-2;
%         topology(d*n+1:n+1:end) = 1;
%         topology(n-1,1) = 1;
%         topology(n,2) = 1; 

for i = 1:n
    for j = 1:n        
        L = j-i;
        
        d = 1;                           % next oscillator communication only
        %d = 2;                          % next-to-next oscillator communication only   
        
        %if L == d                       % uni--directional non-periodic: lower dth subdiagonal 
        %if -L == d                      % uni--directional non-periodic: upper dth subdiagonal 
        %if (L == d) | (L == d-1)        % uni--directional non-periodic compact: all upper subdiagonals till d (here d=2) 
        %if (-L == d) | (L == n-d)       % uni--directional periodic: upper dth and lower (n-d)th subdiagonals 
        if abs(L) == d                   % bi-directional non-periodic: both upper and lower dth subdiagonals
        %if (abs(L) == d) | (abs(L) == n-d)  % bi-directional periodic: both upper and lower dth and (n-d)th subdiagonals    
            topology(j,i) = 1;
        end
    end
end

    % rng(0)
    stop = 0;
    while ~stop
        init_controls
        init_figure
        init_plot
        shg
        t = 0;      
        loop = 1;
        flag = 0;
        while loop
            s = max(step,.01);
            options = odeset('outputfcn',@outputfcn, ...
                'maxstep',s,'initialstep',s);         
            if isinf(quit)
                tspan = t:s:t+delta; % s=timestep, delta=number of timesteps
            else
                tspan = t:s:t+quit;  % infinite run i.e., number of timesteps=inf
            end
            [t,theta] = ode45(@ode,tspan,theta,options);
            hold on;
            t = t(end);
            theta = theta(end,:)';
        end
    end
    if wanth
        Hout = H;
    end
    close(gcf)
   
%-------------------------------------------------------------
%% Kuramoto's system of odes
    function theta_dot = ode(~,theta)
        % The nonlinear term is the gradient of the potential, partial(v)/partial(theta). 
        theta_dot = omega  - kappa/n*gradv(theta);
        
          % noise injection for certain time  
%         if t < 700           %t = 0:s/100:5000
%             %random noise with a standard deviation of 2
%              theta_dot = theta_dot + (5*theta_dot/100).*rand(n,1); %right way
%             % theta_dot = theta_dot + 0.2*sqrt(t).*rand(n,1); 
%              % theta_dot = theta_dot + 500*sin(2*pi*0.1*t) .* rand(n,1);  %  dependent on time step 
%             % Gaussian Pulse
%             %theta_dot = theta_dot + 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01))); 
%             %theta_dot = theta_dot + 2*sin(2*pi*15*t) + 2.5*gallery('normaldata',size(t),4)
%         end

    end

    % Cases: 
        % 1. sync/unstable desync/phase-locked stable desyn
            % potential: sin / tanh
        % 2. sync/unstable desync (never get phase-locked stable desync system)
            % cos/-cos/sec/-sec/-sech
        % 5. sync/phase-locked stable desync
            % sinh
        % 3. unstable desync/phase-locked stable desync system
            % -sin, tan, -tan, -tanh, sech
        % 4. interaction with strongest possible force, i.e., immediate speed-up or slow-down
        % depending on direction (modify frequency of ossilator by a strong non differential step)
            % sign, kroneckerDelta, dirac, sign(sin), heaviside, square, rectpuls,
            % rectangularPulse, tripuls, sawtooth, triangularPulse, ceil, floor
            % x = square(t) generates a square wave with period 2π for the elements of the time array t. 
            % square is similar to the sine function but creates a square wave with values of –1 and 1.
            % x = square(t,duty) generates a square wave with specified duty cycle. 
            % The duty cycle (0-100) is the percent of the signal period in which the square wave is positive.

% tanh(round(theta-theta'))): lock-phased stable desync (bottleneck
% applications with WF_amplitude=round_factor) 
% intitial injection of -3 on second ossilotor   -> pull third ossilator by -0.5 (second ossilator=-2.5) 
%                                                -> pull forth ossilator by -0.5 (second ossilator=-1.5)
%                                                -> pull fifth ossilator by -0.5 (second ossilator=-1)
%                                                   second ossilator=-0.5)

% tanh(theta-theta'): sync (decay of idle waves in applications)
% intitial injection of -3 on second ossilotor   -> pull all ossilators by
% different forces (more stongly nearest ossilator and more weakly distant neighbour)

% discrete functions like sign / square wave etc. there is no interation b/w
% ossilators without decay
    
    function g = gradv(theta)
        % theta-theta' is a matrix with elements theta(j)-theta(k).
        % The sum is by colums and produces a column vector.   
        x=theta-theta';
        
        % random unsuccessful tries
        % x=abs(x);
        % x=sign(x);
        % x=round(x); 
        % x=eq(x,zeros(n)); % let the ossilators interact that has zero phase-differences 
        
        % engineered wavefront by using unitstep with everything shifted to the right by 0.5
        % interaction length for phases (= initial_injection/n)
       %  x(x>=-0.2) = 0;  
       %  g = dot(topology,tanh(x),2);
        % or
        % f = @(t) sin(t).*(t<-0.5);         %t = linspace(-4*pi,4*pi); plot(t, f(t)) axis([4*pi  4*pi    -1.5  1.5])
        % g = dot(topology,f(x),2);
        alpha=1.5;
        sigma=(pi/alpha)*1.5;
        f = @(t) -sin(alpha*t).*(abs(t)<sigma) + 1.*(t>sigma) -1.*(t<-sigma);

        % Half-Wave Rectified Sine Function
         %f = @(t) sin(t).*(sin(t)>=0);        
        %t = linspace(-4*pi,4*pi); plot(t, f(t)) axis([4*pi  4*pi    -1.5  1.5])                                                                            
         g = dot(topology,f(x),2); 
        
        % original global potential / interaction matrix (any non-zero value change frequencies)
        % g = sum(sin(x),2);         
        
        % modified with topology
         %g = dot(topology,-sin(2*x),2);  
        
        % modified for scalable programs with tanh()
         %g = dot(topology,tanh(x),2);  
        
        % modified for botlleneck programs with -tanh()* gaussian
          % gaussian = @(t,amp,mu,sigma) amp*exp(-(((t-mu).^2)/(2*sigma.^2))); 
          % g = dot(topology, tanh(x) * gaussian(t,1,0,40), 2); %t, amp= 1 / 1/(4*0.1*sqrt(2*pi)), mu=0, sigma=0.1
           %expo = @(t,interaction_length) exp(-interaction_length * t);
           %g = dot(topology, tanh(x) * expo(x,0.9), 2); %t, amp= 1 / 1/(4*0.1*sqrt(2*pi)), mu=0, sigma=0.1
    end

%     function v = potential(theta)
%         v = 0;
%         for k = 1:n
%             j = k+1:n;
%             v = v + sum(sin((theta(j)-theta(k))/2).^2);       
%         end
%         v = (4/n^2)*v;
%     end
    
         
%% 
%-------------------------------------------------------------

    function status = outputfcn(t,theta,odeflag)
        % Called after each successful step of the ode solver.
        
        if isequal(odeflag,'init') && wanth
            % assert(t(1) == 0)
            H.t = 0; 
            H.theta = theta;
            H.pot = potential(theta);
            H.gradv = gradv(theta);
            H.psi = psi;
            H.order = 0-0;
        end
        if isempty(odeflag)  % Not 'init' or 'last'.
            fileID = fopen('open_d1_d2.txt','a+'); %'w'
            for j = 1:length(t)
                if loop == 0
                    break
                end
                % Order parameter.  
                z = 1/n*sum(exp(1i*theta(:,j)));
                psi = 0;
                if get(rotate,'value') == 1
                    % Rotating frame of reference.
                    psi = angle(z);
                    z = abs(z);
                end
                for k = 1:n
                    set(point(k), ...
                        'xdata',r(k)*cos(theta(k,j)-psi), ...
                        'ydata',r(k)*sin(theta(k,j)-psi))
                end

                % Length of arrow is order parameter.
                arrow(0,z);

                if pot
                    % Potential
%                     v = potential(theta(:,j));
%                     set(titl,'string', ...
%                         ['potential (' sprintf('%6.3f',v) ')'])
%                     addpoints(animal,t(j),v)
                    v=sum(dot(topology,tanh(abs(theta(:,j)-theta(:,j)')),2));
                    
                    %v = sum(abs(gradv(theta))); % sum(abs(dot(topology,tanh(theta-theta'),2)))
                    g = dot(topology,theta(:,j)-theta(:,j)',2)
                    %if g(1,1) > 0.004
                    if abs(v) > 0.004
                        format long
                        fprintf(fileID,'%6.3f %6.3f  %6.3f \n',t,g(1,1),v);
                        %fprintf(fileID,'%6.3f %6.3f \n',t,v);
                        %save('potential','t','v');
                    end
                    set(titl,'string', ...
                        ['timeline of phase differences, potential (=' sprintf('%6.3f',v) ')'])
                    for k = 1:n
                    addpoints(animal,t(j),g(k))
                    end
                else
                    % Order parameter
                    set(titl,'string', ...
                        ['order parameter (' sprintf('%6.3f',abs(z)) ')'])
                    addpoints(animal,t(j),abs(z))
                end
                
                if wanth
                    nt = length(H.t) + 1;
                    H.t(nt) = t(j);
                    H.theta(:,nt) = theta(:,j);
                    H.pot(nt) = potential(theta(:,j));
                    H.gradv(:,nt) = gradv(theta(:,j)); 
                    H.psi(nt) = psi;
                    H.order(nt) = abs(z);
                end

                if isinf(quit), qt = delta; else, qt = quit; end
                if rem(t(j),qt) < 1 && t(j) > 10 || t(j) >= quit
                    clearpoints(animal)
                    set(ax2,'xlim',[t(j)-1 t(j)+qt])
                    if t(j) >= quit
                        loop = 0;
                        stop = 1;
                    end
                end
            end
            fclose(fileID);
        end
        status = flag + stop;
        drawnow limitrate
    end

    function init_vars(preset)
        width = 0;   % Standard deviation of display ring.
        step = 0.1;  % Display step size.
        pot = true;  % Graph potential or order parameter.
        uni = true;  % Uniform distribution of omega.
        quit = inf;  % 0 <= t <= quit.
        %delta = 200*(1+double(preset==3));
         delta = 100;
        if wanth
            quit = delta;
        end
        
%-------------------------------------------------------------
%% presets 
        if preset > 0
            n = 5;
        end
        switch preset
            case 0              % strong coupling 
                kappa = 1;        %d_avg/Tcomp+Tcomm
                beta = 0;
                step = 0.1;
                n = 8;
            case 1              % free ossilations 
                kappa = 0;
                beta = 0;
                step = 1;
            case 2              % strong coupling 
                kappa = .75;  
                beta = 0;
                step = 1;
            case 3
                kappa = .36;    % strong coupling compare to frequency spread of oscillators
                beta = .24;
                step = 0.5;
            case 4
                kappa = .36;    % weak coupling compare to frequency spread of oscillators
                beta = .23;
                step = 0.5;
            case 5
                uni = false;    %100 oscillators
                n = 100;  
                kappa = .10;
                beta = .05;
                step = 0.5;
            otherwise
                display(['unknown preset:' int2str(preset)])
                scream
        end
%-------------------------------------------------------------
%% 
        if wanth
            H.n = n;
            H.kappa = kappa;
            H.beta = beta;
            H.t = [];
            H.theta = [];
            H.pot = [];
            H.gradv = [];
            H.psi = [];
            H.order = [];
        end
        loop = 0;
    end

    function txtval(v,vmin,vmax,fmt,cb,k)
        % Slider with text.
        txt(k) = uicontrol('style','text', ...
            'string',sprintf(fmt,v), ...
            'units','normal', ...
            'position',[.04 0.92-.10*k .18 .05], ...
            'background',[.94 .94 .94], ...
            'fontsize',get(0,'defaultuicontrolfontsize')+2, ...
            'horiz','left', ...
            'background','w');
        uicontrol('style','slider', ...
            'units','normal', ...
            'position',[.04 0.88-.10*k .18 .05], ...
            'min',vmin, ...
            'max',vmax, ...
            'value',v, ...
            'callback',cb);
    end

    function t = toggle(str,v,cb,k,lr)
        % Toggle switches.
        switch lr
            case 'l'  % left side
                x = .04;
                dx = .18;
                y = .34;
                dy = .07;
            case 'r'  % right side
                x = .90;
                dx = .08;
                y = .98;
                dy = .08;
        end
        t = uicontrol('style','toggle', ...
            'units','normal', ...
            'position',[x y-k*dy dx dy-.02], ...
            'string',str, ...
            'value',v, ...
            'callback',cb);
    end
        
    function init_controls
        % Initialize buttons, sliders and toggles.
        shg
        % Preset buttons
        uicontrol('style','text', ...
            'string','preset', ...
            'units','normal', ...
            'horiz','left', ...
            'background','w', ...
            'position', [.04 .92 .10 .04])
        for k = 1:5
            uicontrol('style','radio', ...
                'units','normal', ...
                'position',[.04*k .88 .04 .04], ...
                'background','w', ...
                'value',double(k == preset), ...
                'callback',@radiocb)
        end
        
        % Sliders
        txt = zeros(5,1);
        txtval(n,1,100,'n = %3d',@ncb,1);
        txtval(kappa,0,1.0,'kappa = %5.3f',@kappacb,2);
        txtval(beta,0,1.0,'beta = %5.3f',@betacb,3);
        txtval(step,0,1,'step = %5.3f',@stepcb,4);
        txtval(width,0,0.2,'width = %5.3f',@widthcb,5);
                
        % Toggles
        toggle('restart',0,@restartcb,1,'l');
        rotate = toggle('rotate',1,[],2,'l');
        if uni
            toggle('uniform / random',0,@unicb,3,'l');
        else
            toggle('random / uniform',1,@unicb,3,'l');
        end
        if pot
            toggle('potential / order',1,@potcb,4,'l');
        else
            toggle('order / potential',0,@potcb,4,'l');
        end
        toggle('exit',0,@stopcb,1,'r');
        toggle('help',0,@helpcb,2,'r');
        toggle('blog',0,@blogcb,3,'r');
        flag = 0;
        stop = 0;
    end

    function init_figure(~)
        % Initialize figure window.
        set(gcf,'menubar','none','numbertitle','off', ...
             'name','kuramoto','color','white')
        
        ax1 = axes('position',[.30 .34 .60 .60]);
        circle = exp((0:.01:1)*2*pi*1i);
        line(real(circle),imag(circle),'color',grey)
        axis(1.2*[-1 1 -1 1]) % change circle radius for figure 1
        axis square
        set(gca,'xtick',[],'ytick',[])
        box on
        
        ax2 = axes('position',[.3 .07 .6 .2]); % change dimentions for potential figure 2
        %animal = animatedline('linewidth',2, ...
            %'color',cyan);
        animal = animatedline('lineStyle','none','Marker','o','MarkerSize',1,'color',cyan); 
        if isinf(quit)
            %axis([0 delta 0 1.2]) %  change axis of potential figure 2
            axis([0 delta -2*pi 2*pi]) 
        else
            %axis([0 quit 0 5.2])
            axis([0 quit -2*pi 2*pi])
        end
        titl = title('');
        box on
        
        axes(ax1)
    end

%-------------------------------------------------------------
%% Initial condition   
    function init_plot
        % Initialize plot.
        
        % Oscillators initially uniform around circle.
        %theta = zeros(n,1);               % lock-step symmetry 
        %theta = (1:n)'/n*2*pi;            % uniform distribution on all ossilators in [0,2π] 
        %theta = 2*pi.*rand(n,1) -pi;      % uniform distribution  on all ossilators in [-π,π]
        theta = 2*pi*rand(n,1);           % random distribution  on all ossilators in [0,2π]      
        %theta = double((1:n)' <1.1) * 4; %-3*pi/2;    % disturbance on a single ossilator
        omega = omegas(uni);       
        
        % Plot radii.
        r = ones(n,1) + width*randn(n,1);
        
        % Parula, blue is fast, yellow is slow.
        colors = flipud(parula(ceil(n)));
        
        point = zeros(n,1);
        for k = 1:n
           point(k) = line(r(k)*cos(theta(k)),r(k)*sin(theta(k)), ...
                'linestyle','none', ...
                'marker','o', ...
                'markersize',6, ...
                'markeredgecolor','k', ...
                'markerfacecolor',colors(mod(k-1,length(colors))+1,:));
        end
        
        title(sprintf('kappa = %5.3f, beta = %5.3f',kappa,beta))
    end

    function omega = omegas(uni)
        % Intrinsic freqencies in interval of width beta centered at 1.
        if uni
            omega = ones(n,1) + beta*(-1:2/(n-1):1)';
        else
            omega = ones(n,1) + beta*(2*rand(n,1)-1);
        end
    end

%%  
    function ncb(arg,~)
        % Number of oscillators.
        n = round(get(arg,'value'));
        set(arg,'value',n)
        set(txt(1),'string',sprintf('n = %d',n))
        loop = 0;
    end

    function kappacb(arg,~)
        % Coupling parameter.
        kappa = get(arg,'value');
        set(txt(2),'string',sprintf('kappa = %5.3f',kappa))
    end

    function betacb(arg,~)
        % Spread of intrinsic frequencies.
        beta = get(arg,'value');
        set(txt(3),'string',sprintf('beta = %5.3f',beta))
        omega = omegas(uni);
    end

    function radiocb(arg,~)
        pos = get(arg,'position');
        preset = pos(1)/.04;
        init_vars(preset)
    end

    function stepcb(arg,~) 
        % Step size for display.
        step = get(arg,'value');
        set(txt(4),'string',sprintf('step = %5.3f',step))
        flag = 1;
    end 

    function widthcb(arg,~)
        % Standard deviation of random radii.
        oldw = width;
        width = get(arg,'value');
        set(txt(5),'string',sprintf('width = %5.3f',width))
        r = (r - 1)*width/(oldw+realmin) + 1;
        flag = 1;
    end

    function restartcb(~,~)
        % Restart with current parameters.
        clearpoints(animal)
        drawnow
        loop = 0;
        %flag = 1;
    end

    function unicb(arg,~)
        uni = get(arg,'value') == 1;
        if uni
            set(arg,'string','uniform / random');
        else
            set(arg,'string','random / uniform');
        end
        omega = omegas(uni);
    end

    function potcb(arg,~)
        pot = get(arg,'value') == 1;
        if pot
            set(arg,'string','potential / order')
        else
            set(arg,'string','order / potential')
        end
    end

    function stopcb(~,~)
        loop = 0;
        stop = 1;
    end

    function helpcb(~,~)
        doc('kuramoto')
        toggle('help',0,@helpcb,2,'r');
    end

    function blogcb(~,~)
        web(['https://blogs.mathworks.com/cleve/2019/10/30' ...
             '/stability-of-kuramoto-oscillators'])
        toggle('blog',0,@blogcb,3,'r');
    end

    function arrow(z0,z1)
        delete(findobj('tag','arrow_shaft'))
        delete(findobj('tag','arrow_head'))
        rho = angle(z1-z0);
        x = real(z1);
        y = imag(z1);
        u = [0 -.08 -.05 -.08 0];
        v = [0 -.05 0 +.05 0];
        s = u;
        u = u*cos(rho) - v*sin(rho) + x;
        v = s*sin(rho) + v*cos(rho) + y;
        line([real(z0) real(z1)],[imag(z0) imag(z1)], ...
            'linewidth',1.5, ...
            'color',cyan, ...
            'tag','arrow_shaft');
        patch(u,v,cyan, ...
            'tag','arrow_head');
    end
            
    function c = grey
        c = [.8 .8 .8];
    end
            
    function c = cyan
        c = [0 .6 .6];
    end
end