clear, clc, close all

% name animated gif
filename = 'f_number.gif';

% define F/#
F = 8;

% octagonal pupil? else circular
oct = true;

% parameterize simulation
N = 2048;           % grid points per side [#]
wvl = 10.6e-6;      % optical wavelength [m]
k = 2*pi/wvl;       % angular wavenumber [rad/m]
f = 50e-3;          % effective focal length [m]
S = sqrt(N*wvl*f);  % side length [m]
D = f/F;            % entrance-pupil diameter [m]
delta = S/N;        % grid spacing [m]
[x,y] = meshgrid((-N/2:N/2-1)*delta);

% define analytic function
syms chat(p)
chat(p) = real(2/pi*(acos(abs(p))-abs(p)*sqrt(1-abs(p)^2)));

% create storage vectors
tiltVector = 0:0.1:atand(1/(2*F));
powerVector = zeros(1,length(tiltVector));

if oct % octagonal pupil
    H = D/2*sqrt(1/2*(1+sqrt(2))*pi);
    square = (abs(x)<=H/2).*(abs(y)<=H/2);
    square45 = imrotate(square,45,'crop');
    amp_lns = square.*square45;
else % circular pupil
    amp_lns = x.^2 + y.^2 <= (D/2)^2;
end
phs_lns = -k*((x.^2+y.^2)/(2*f)); % thin-lens phase

% tricks to prevent wraparound...
mask0 = x.^2 + y.^2 <= (1.25*D/2)^2;
flag0 = false;
idx0 = 0;

for iter = 1:length(tiltVector)
    
    % define collimated field incident upon entrance pupil and add tilt
    fld_pup1 = amp_lns.*exp(1j*phs_lns);
    fld_pup1 = fld_pup1.*exp(1j*k*x*tand(tiltVector(iter)));

    % propagate transmitted field from entrance-pupil plane to focal plane
    fld_foc = propF_TF(fld_pup1,S,wvl,f);

    % compute irradiance in focal plane
    irr_foc = abs(fld_foc).^2;

    % propagate retroreflection from focal plane to entrance-pupil plane
    fld_pup2 = propF_TF(fld_foc,S,wvl,f);
    fld_pup2 = fld_pup2.*amp_lns;

    % this block just crudely prevents wraparound...
    % feel free to comment out and see what happens at low F-numbers
    mask = circshift(mask0,round(2*f*tand(tiltVector(iter))/delta),2);
    if mask(N/2,1) ~= 0
        flag0 = true;
        idx0 = (find(ischange(find(mask(N/2,:)),'linear'))-1);
        mask(:,1:idx0) = 0;
    else
        flag0 = false;
    end
    if ~flag0 && idx0 ~= 0
        mask = zeros(size(mask));
    end
    fld_pup2 = fld_pup2.*mask;

    % compute irradiance in entrance-pupil plane
    irr_pup = abs(fld_pup2).^2;

    % normalize and store total power in entrance-pupil plane
    powerVector(iter) = sum(sum(irr_pup));
    if iter == 1
        norm_factor = powerVector(iter);
    end
    powerVector(iter) = powerVector(iter)/norm_factor;

    % plot results
    h = figure(1);
    h.Position = [0 0 980 735];
    clf
    tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact')

    nexttile
    imagesc((irr_foc((-round(N/F)+1:round(N/F))+N/2,(-round(N/F)+1:round(N/F))+N/2).^(1/2)));
    axis image
    axis off
    hold on
    colormap(hot)
    if iter == 1
        pk1 = max(max(irr_foc.^(1/3)));
    end
    clim([0 pk1])
    title(['$\theta =$ ' sprintf('%.1f',tiltVector(iter)) '$^{\circ}$'],'Interpreter','latex','FontSize',16)

    nexttile
    imagesc(irr_pup);
    axis image
    axis off
    colormap(hot)
    if iter == 1
        pk2 = max(max(irr_pup));
    end
    clim([0 pk2])
    title(['ground truth: $F/$' sprintf('%.1f',F)],'Interpreter','latex','FontSize',16)

    nexttile
    fplot(@(theta) chat(2*F*tan(pi*theta/180)).*cos(pi*theta/180),[-atand(1/(2*F)) atand(1/(2*F))],'k--','LineWidth',2)
    hold on
    plot(tiltVector(1:iter),powerVector(1:iter),'LineWidth',2,'Color',"#D95319")
    xlim([0 atand(1/(2*F))])
    ylim([0 1])
    axis square
    set(gca, 'LineWidth', 1.5, 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', 'none', 'Xcolor', [0, 0, 0], 'YColor', [0, 0, 0], ...
        'TickLabelInterpreter', 'latex','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XTick',[-atand(1/(2*F)) 0 atand(1/(2*F))],'XTickLabel',{'$-\arctan \left( \frac{1}{2 F/\#} \right)$','0','$\arctan \left( \frac{1}{2 F/\#} \right)$'});
    xlabel('$\theta$ [$^\circ$]', 'Interpreter', 'latex')
    ylabel('$P/P_0$', 'Interpreter', 'latex')
    drawnow

    % initialize or append to animated gif
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if iter == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
    end

end

% compute F-number by fitting result to curve and display on plot
fitfun0 = fittype( @(F,b,x) double(chat(2*F*tan((x-b)*pi/180)).*cos((x-b)*pi/180)) );
opts = fitoptions( 'Method','NonlinearLeastSquares','StartPoint',[1 1]);
[f0,~] = fit(tiltVector',powerVector',fitfun0,opts);
title(['curve fit: $F/$' sprintf('%.1f',abs(f0.F))],'Interpreter','latex','FontSize',16)