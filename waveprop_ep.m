clear, clc, close all

% name animated gif
filename = 'entrance_pupil.gif';

% define F/#
F = 5.6;

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
syms edge(w)
edge(w) = real(1/pi*(acos(w)-w*sqrt(1-w^2)));

% create storage vectors
scanVector = -D/2:0.1e-3:D/2;
powerVector = zeros(1,length(scanVector));

if oct % octagonal pupil
    H = D/2*sqrt(1/2*(1+sqrt(2))*pi);
    square = (abs(x)<=H/2).*(abs(y)<=H/2);
    square45 = imrotate(square,45,'crop');
    amp_lns = square.*square45;
else % circular pupil
    amp_lns = x.^2 + y.^2 <= (D/2)^2;
end
phs_lns = -k*((x.^2+y.^2)/(2*f)); % thin-lens phase

for iter = 1:length(scanVector)

    % define collimated field incident upon entrance pupil and scan
    fld_pup = amp_lns.*exp(1j*phs_lns);
    fld_pup = fld_pup.*(x>=scanVector(iter));

    % compute irradiance in entrance-pupil plane
    irr_pup = abs(fld_pup).^2;

    % propagate transmitted field from entrance-pupil plane to focal plane
    fld_foc = propF_TF(fld_pup,S,wvl,f);

    % compute irradiance in focal plane
    irr_foc = abs(fld_foc).^2;

    % normalize and store total power in focal plane
    powerVector(iter) = sum(sum(irr_foc));
    if iter == 1
        norm_factor = powerVector(iter);
    end
    powerVector(iter) = powerVector(iter)/norm_factor;

    % manipulate entrance-pupil plane purely for visualization purposes
    irr_pup(irr_pup==0) = 3;
    irr_pup(x<=scanVector(iter)) = 0;

    % plot results
    h = figure(1);
    h.Position = [0 0 560*1.75 420*1.75];
    clf
    tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact')

    nexttile
    imagesc(irr_pup);
    axis image
    axis off
    title(['$\Delta x$ = ' sprintf('%.1f',(scanVector(iter)*1e3)) ' mm'],'Interpreter','latex','FontSize',16)
    if iter == 1
        pk1 = max(max(irr_pup));
    end
    clim([0 pk1])
    hold on
    xline((scanVector(iter)+S/2)/delta,'k-','LineWidth',2)

    nexttile
    imagesc((irr_foc((-N/64+1:N/64)+N/2,(-N/64+1:N/64)+N/2)).^(1/2));
    if iter == 1
        pk2 = max(max(irr_foc.^(1/2)));
    end
    clim([0 pk2])
    axis image
    axis off
    colormap(gray)
    title(['ground truth: $D_\mathrm{EP} =$ ' sprintf('%.1f',D*1e3) ' mm'],'Interpreter','latex','FontSize',16)

    nexttile
    fplot(@(x) edge(2*x*1e-3/D),[-D/2*1e3 D/2*1e3],'k--','LineWidth',2)
    hold on
    plot(scanVector(1:iter)*1e3,powerVector(1:iter),'LineWidth',2,'Color',"#0072BD")
    xlim([-D/2*1e3 D/2*1e3])
    ylim([0 1])
    axis square
    set(gca, 'LineWidth', 1.5, 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', 'none', 'Xcolor', [0, 0, 0], 'YColor', [0, 0, 0], ...
        'TickLabelInterpreter', 'latex','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off','XTick',[-D/2*1e3 0 D/2*1e3],'XTickLabel',{'$-\frac{D_\mathrm{EP}}{2}$','0','$\frac{D_\mathrm{EP}}{2}$'});
    xlabel('$\Delta x$ [mm]', 'Interpreter', 'latex')
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

% compute entrance-pupil diameter by fitting result to curve and display on plot
fitfun0 = fittype( @(D,b,x) double(edge(2*(x-b)/D)) );
opts = fitoptions( 'Method','NonlinearLeastSquares','StartPoint',[1e-3 1e-3]);
[f0,~] = fit(scanVector',powerVector',fitfun0,opts);
title(['curve fit: $D_\mathrm{EP} =$ ' sprintf('%.1f',f0.D*1e3) ' mm'],'Interpreter','latex','FontSize',16)