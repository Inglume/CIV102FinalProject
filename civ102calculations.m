%% 0. Initialize Parameters
L = 1280; % Length of bridge
n = L + 1; % Number of locations to evaluate bridge failure
x = linspace(0, L, n); % Define x coordinate
SFD = zeros(1, n); % Initialize SFD(x)
BMD = zeros(1, n); % Initialize SFD(x)

%% 1. Point Loading Analysis (SFD, BMD)

P = 1; % temporary, just to help calculate Pfail for two point loads


[SFDTrain, BMDTrain] = ApplyTrainLoad(x);
[SFD2L, BMD2L] = ApplyTwoLoads(1, x, SFD, BMD);

PlotTrain(x, L, SFDTrain, BMDTrain)
figure()
PlotDiagrams(x, L, SFD2L, BMD2L)

SFDTrain = max(SFDTrain(1, :), SFDTrain(2, :));
BMDTrain = max(BMDTrain(1, :), BMDTrain(2, :));

%% 2. Define cross-sections

% xc : Location, x, of cross-section change
% bft : Top Flange Width
% tft : Top Flange Thickness
% hw : Web Height
% tw : Web Thickness (Assuming 2 separate webs)
% bfb : Bottom Flange Width
% tfb : Bottom Flange Thickness
% a : Diaphragm Spacing

% Design 0

%% UNCOMMENT THIS IF YOU WANNA SEE HOW THE OUTPUT SHOULD WORK
% (each of the rows will have the same properties though since only a
% changes over the distance)

GeometricInputs0 = [];

GeometricInputs0(end + 1, :) = [0, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs0(end + 1, :) = [30, 100, 1.27, 72.46, 1.27, 80, 1.27, 520];
GeometricInputs0(end + 1, :) = [550, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs0(end + 1, :) = [580, 100, 1.27, 72.46, 1.27, 80, 1.27, 480];
GeometricInputs0(end + 1, :) = [1060, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs0(end + 1, :) = [1090, 100, 1.27, 72.46, 1.27, 80, 1.27, 160];
GeometricInputs0(end + 1, :) = [1280, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs0(end + 1, :) = [L, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];

CrossSectionProperties0 = SectionProperties(GeometricInputs0, n);

% Design 1.0000003

GeometricInputs = [];

GeometricInputs(end + 1, :) = [0, 100, 2.54, 80, 1.27, 75, 1.27, 5]; % CS #1 with plate
GeometricInputs(end + 1, :) = [5, 100, 2.54, 80, 1.27, 75, 1.27, 20];
GeometricInputs(end + 1, :) = [25, 100, 2.54, 80, 1.27, 75, 1.27, 5];
GeometricInputs(end + 1, :) = [30, 100, 2.54, 80, 1.27, 75, 0, 250]; % CS #1
GeometricInputs(end + 1, :) = [280, 100, 2.54, 80, 1.27 * 3 / 2, 75, 0, 275]; % CS #2
GeometricInputs(end + 1, :) = [555, 100, 2.54, 80, 1.27 * 3 / 2, 75, 0, 20];
GeometricInputs(end + 1, :) = [575, 100, 2.54, 80, 1.27 * 3 / 2, 75, 0, 113];
GeometricInputs(end + 1, :) = [688, 100, 2.54, 80, 1.27, 75, 0, 115]; % CS #1
GeometricInputs(end + 1, :) = [803, 100, 2.54, 80 - 1.27, 1.27 * 3 / 2, 75, 2.54, 262]; % CS #3
GeometricInputs(end + 1, :) = [1065, 100, 2.54, 80 - 1.27, 1.27 * 3 / 2, 75, 2.54, 20];
GeometricInputs(end + 1, :) = [1085, 100, 2.54, 80 - 1.27, 1.27 * 3 / 2, 75, 2.54, 170];
GeometricInputs(end + 1, :) = [1255, 100, 2.54, 80 - 1.27, 1.27, 75, 2.54, 20]; % #4
GeometricInputs(end + 1, :) = [1275, 100, 2.54, 80 - 1.27, 1.27, 75, 2.54, 5];
GeometricInputs(end + 1, :) = [L, 100, 2.54, 80 - 1.27, 1.27, 75, 2.54, 5];

CrossSectionProperties = SectionProperties(GeometricInputs, n);

%% 3. Define Material Properties
SigT = 30; % tensile stress
SigC = 6; % compressive stress
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

for i = 1 : size(GeometricInputs0, 1) - 1
    cs = CrossSectionProperties0(GeometricInputs0(i, 1) + 1, :);
    sprintf("Cross Section @ %d mm - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Qmax: %.3g Qglue: %.3g ", GeometricInputs0(i, 1), cs(8:11), cs(16))
end


%% 4. Calculate Failure Moments and Shear Forces

Fails0 = GetFails(CrossSectionProperties0, TauU, TauG, E, mu, SigT, SigC, BMD2L);
Fails = GetFails(CrossSectionProperties, TauU, TauG, E, mu, SigT, SigC, BMD2L);

%% 4.7 Calculate Failure Load
Pf0 = FailLoad(SFD2L, BMD2L, Fails0)
%Pf0 = FailLoad(SFDTrain, BMDTrain, Fails0)
Pf = FailLoad(SFD2L, BMD2L, Fails)
%Pf = FailLoad(SFDTrain, BMDTrain, Fails)

Plot2L(x, L, Pf0, CrossSectionProperties0, TauU, TauG, E, mu, SigT, SigC)
Plot2L(x, L, Pf, CrossSectionProperties, TauU, TauG, E, mu, SigT, SigC)

%% 5. Curvature, Slope, Deflections
%Defls = Deflections(x, BMDTrain, GeometricInputs(10), E);

%% Functions

%{
function [SFD, BMD] = ApplyPL(xP, P, x, SFD, BMD) % don't need this
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous SFD can be entered as input to
% construct SFD of multiple point loads
%   Output: SFD, BMD both 1-D arrays of length n
% need to account for support reaction forces too
    xA = 15; % location of support A
    xB = 1075; % location of support B
    By = P * (xP - xA) / (xB - xA);
    Ay = P - By; % force equilibrium on x
    
    [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD); % add left support point load
    [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD);% add right support point load
    [SFD, BMD] = UpdateDiagrams(xP, -P, x, SFD, BMD);% add applied force point load
end
%}

function [SFD, BMD] = ApplyTwoLoads(P, x, SFD, BMD) % P is force of each individual load
    xA = 15; % location of support A
    xB = xA + 1060; % location of support B
    xP1 = xA + 550; % applied force 1
    xP2 = xB + 190; % applied force 1
    By = (P * (xP1 - xA) + P * (xP2 - xA)) / (xB - xA); % support B located at 1050mm from support A, support A located at x = 0
    Ay = 2 * P - By; % force equilibrium on x
    
    [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD); % add left support force
    [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD); % add right support force
    [SFD, BMD] = UpdateDiagrams(xP1, -P, x, SFD, BMD); % add first applied force
    [SFD, BMD] = UpdateDiagrams(xP2, -P, x, SFD, BMD); % add second applied force
end

function [trainSFD, trainBMD] = ApplyTrainLoad(x)
% Constructs SFD and BMD from application of train Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous SFD can be entered as input to
% construct SFD of multiple point loads
%   Output: SFD, BMD both 2-D arrays of length n and height i have no idea
% need to account for support reaction forces too
    P = 0.4;
    xA = 15; % location of support A
    xB = 1075; % location of support B
    trainLength = 856; % not including length past wheels at end
    wheelSpacing = [0, 176, 340, 516, 680, 856]; % location of wheels relative to backmost ones
    trainSFD = zeros(2, size(x, 2));
    trainBMD = zeros(2, size(x, 2));
    position = [117, 424]; % two cases of train loading
    
    for i = 1 : 2 % do for two cases, not this ppopoo
        SFD = zeros(1, size(x, 2));
        BMD = zeros(1, size(x, 2));
        By = 0;
        for j = 1 : 6 % for every set of wheels
            By = By + (P / 6 * (wheelSpacing(j) + position(i) - xA));
        end
        By = By / (xB - xA);
        Ay = P - By;
        [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD);
        [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD);
        for j = 1 : 6 % for every set of wheels
            [SFD, BMD] = UpdateDiagrams(wheelSpacing(j) + position(i), -P / 6, x, SFD, BMD);
        end
        trainSFD(i, :) = SFD;
        trainBMD(i, :) = BMD;
    end
end

function [SFD, BMD] = UpdateDiagrams(xP, P, x, SFD, BMD)
    SFD(xP + 1) = SFD(xP + 1) + P;  % x location must be increased by one since matrices are 1-indexed
    
    for i = xP + 2 : length(x)
        SFD(i) = SFD(i) + P;
        BMD(i) = BMD(i - 1) + SFD(i - 1); % start changing BMD at (xP + 1) since moment at xP = 0 anyways relative to shear force
    end
end

function PlotDiagrams(x, L, SFD, BMD) % include curvature diagram
    subplot(2, 1, 1) % SFD
    plot(x, SFD)
    xlim([0 L])
    title("Shear Force over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("V (N)")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 1, 2) % BMD
    plot(x, BMD)
    xlim([0 L])
    title("Moment over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("M (N mm)")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference
    
    set(gcf, 'Name', 'SFD, BMD of Load Case #2 (two point load)') % name of window
end

function PlotTrain (x, L, trainSFD, trainBMD) % include curvature diagram
    trainLength = size(trainSFD, 2);
    subplot(2, 1, 1) % SFD
    for i = 1 : size(trainSFD, 1)
        plot(x, trainSFD)
        hold on
    end
    hold off
    %plot(x(15 : 15 + trainLength - 1), SFD)
    xlim([0 L])
    title("Shear Force over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("V (N)")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 1, 2) % BMD
    for i = 1 : size(trainSFD, 1)
        plot(x, trainBMD)
        hold on
    end
    hold off
    %plot(x(15 : 15 + trainLength - 1), BMD)
    xlim([0 L])
    title("Moment over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("M (N mm)")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference
    
    set(gcf, 'Name', 'SFD, BMD of Load Case #1 (train)') % name of window
end

% check this function
function CrossSectionProperties = SectionProperties(GeometricInputs, n) % include y plate 1, 2, 3, 4, 5, 6
% Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
%   Input: Geometric Inputs. Format will depend on user
%   Output: Sectional Properties at every value of x. Each property is a 1-D array of length n
    sect = zeros(n, 9); % columns: ybot, ytop, I, Qmax (yfla, yweb, )
    geom = zeros(n, size(GeometricInputs, 2) - 1);
    for i = 1 : (size(GeometricInputs, 1) - 1) % get row # (for each cross section)
        areas = zeros(1, 3); % length three (top, web, bottom)
        distances = zeros(1, 3); % length three (distance from bot)
        secondMomInert = zeros(1, 3); % second moment of inertias for each part
        areas(1) = GeometricInputs(i, 2) * GeometricInputs(i, 3);
        areas(2) = GeometricInputs(i, 5) * GeometricInputs(i, 4);
        areas(3) = GeometricInputs(i, 6) * GeometricInputs(i, 7);
        distances(3) = GeometricInputs(i, 7) / 2;
        distances(2) = GeometricInputs(i, 4) / 2 + GeometricInputs(i, 7);
        distances(1) = GeometricInputs(i, 3) / 2 + GeometricInputs(i, 4) + GeometricInputs(i, 7);
        secondMomInert(1) = (GeometricInputs(i, 2) * GeometricInputs(i, 3) ^ 3) / 12;
        secondMomInert(2) = (GeometricInputs(i, 5) * GeometricInputs(i, 4) ^ 3) / 12;
        secondMomInert(3) = (GeometricInputs(i, 6) * GeometricInputs(i, 7) ^ 3) / 12;
        ybot = 0; % y bar from bottom
        ytop = 0; % y bar from top
        I = 0; % sum of individual i plus area times distance squared
        Qmax = 0; % Q value from y bar
        yfla = 0; % distance from global centroid to centroid of flange above ybar
        yweb = 0; % distance from global centroid to centroid of web above ybar
        ytab = 0; % distance from global centroid to centroid of tab above ybar
        yfgl = 0; % distance from glue to flange
        Qglue = 0;% Q value from point of glue (on top)
        for k = 1 : 3
            ybot = ybot + areas(k) * distances(k);
            if k == 2
                ybot = ybot + areas(k) * distances(k); % im lazy so yea
                ybot = ybot + 20 * 1.27 * (distances(1) - GeometricInputs(i , 3) / 2 - 1.27 / 2);
            end
        end
        ybot = ybot / (sum(areas) + areas(2) + 20 * 1.27); % hardcode glue tabs cause I am in pain
        ytop = GeometricInputs(i, 3) + GeometricInputs(i, 4) + GeometricInputs(i, 7) - ybot;
        yfla = distances(1) - ybot;
        yweb = (GeometricInputs(i, 4) - ybot + GeometricInputs(i, 7)) / 2;
        ytab = yfla - (GeometricInputs(i, 3) / 2) - (1.27 / 2);
        yfgl = GeometricInputs(i, 3) / 2;
        for k = 1 : 3
            I = I + secondMomInert(k) + areas(k) * (distances(k) - ybot) ^ 2; % prbaobly somethign wrong with this calc
            if k == 2
                I = I + secondMomInert(k) + areas(k) * (distances(2) - ybot) ^ 2; % lazy * 2
                I = I + (20 * 1.27 ^ 3) / 12 + 20 * 1.27 * ytab ^ 2; % hardcode glue tabs again
            end
        end
        Qmax = areas(1) * yfla;
        Qmax = Qmax + GeometricInputs(i, 5) * 4 * yweb ^ 2; % multiply 2 times thickness of web by height (yweb * 2), multiply by distance (yweb)
        Qmax = Qmax + 20 * 1.27 * ytab;
        Qglue = areas(1) * yfla;
        for j = (GeometricInputs(i, 1) + 1) : GeometricInputs(i + 1, 1) + 1 % beginning of current cross section to beginning of next cross section
            geom(j, :) = GeometricInputs(i, 2 : end);
            sect(j, 1) = ybot;
            sect(j, 2) = ytop;
            sect(j, 3) = I;
            sect(j, 4) = Qmax;
            sect(j, 5) = yfla;
            sect(j, 6) = yweb;
            sect(j, 7) = ytab;
            sect(j, 8) = yfgl;
            sect(j, 9) = Qglue;
        end
    end
    CrossSectionProperties = [geom sect];
end

function [V_fail] = Vfail(CrossSectionProperties, TauU)
% Calculates shear forces at every value of x that would cause a matboard shear failure
%   Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
%   Output: V_fail a 1-D array of length n
    %V_fail = zeroes(1, length(x));
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 4) * 2;
    Qcent = CrossSectionProperties(:, 11);
    V_fail = (I .* b ./ Qcent * TauU)';
end

function [V_failGlue] = VfailGlue(CrossSectionProperties, TauG)
% Calculates shear forces at every value of x that would cause a glue shear failure
%   Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
%   Output: V_fail a 1-D array of length n
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 4) * 2 + 20;
    Qglue = CrossSectionProperties(:, 16);
    V_failGlue = (I .* b ./ Qglue * TauG)';
end

function [V_Buck] = VfailBuck(CrossSectionProperties, E, mu)  
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property)
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 4) * 2;
    Qcent = CrossSectionProperties(:, 11);
    h = CrossSectionProperties(:, 2) + CrossSectionProperties(:, 3) + CrossSectionProperties(:, 6);
    a = CrossSectionProperties(:, 7);
    t = CrossSectionProperties(:, 4);

    TauCrit = ((5 .* (pi^2) .* E) ./ (12 .* (1 - (mu^2)))) .* ((t ./ h) .^ 2 + (t ./ a) .^ 2);
    V_Buck = (TauCrit .* I .* b ./ Qcent)';
end

function [M_MatT] = MfailMatT(CrossSectionProperties, SigT, BMD)
% Calculates bending moments at every value of x that would cause a matboard tension failure
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array)
% Output: M_MatT a 1-D array of length n
    M_MatT = zeros(1, size(BMD, 2));
    I = CrossSectionProperties(:, 10);
    ybot = CrossSectionProperties(:, 8);
    ytop = CrossSectionProperties(:, 9);

    for i = 1 : size(BMD, 2)
        if BMD(i) >= 0 % If the moment is positive, the tension failure will be at the bottom
            M_MatT(i) = SigT * I(i) / ybot(i);
        else % If the moment is negative, the tension failure will be at the top
            M_MatT(i) = -SigT * I(i) / ytop(i);
        end
    end
end

function [M_MatC] = MfailMatC(CrossSectionProperties, SigC, BMD) % Similar to MfailMatT
    M_MatC = zeros(1, size(BMD, 2));
    I = CrossSectionProperties(:, 10);
    ybot = CrossSectionProperties(:, 8);
    ytop = CrossSectionProperties(:, 9);

    for i = 1 : size(BMD, 2)
        if BMD(i) >= 0 % If the moment is positive, the compression failure will be at the top
            M_MatC(i) = SigC * I(i) / ytop(i);
        else % If the moment is negative, the compression failure will be at the bottom
            M_MatC(i) = -SigC * I(i) / ybot(i);
        end
    end
end

function [M_Buck] = MfailBuck1(CrossSectionProperties, E, mu, BMD) % case 1, middle part of flange
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    M_Buck = zeros(1, size(BMD, 2)); % initialize list
    I = CrossSectionProperties(:, 10);
    b = (75 - CrossSectionProperties(:, 4) .* 2); % width of middle part of flange
    yfla = CrossSectionProperties(:, 9);
    t = CrossSectionProperties(:, 4);

    fail = zeros(1, size(BMD, 2)); % initialize array
    for i = 1 : size(BMD, 2)
        fail(i) = ((4*(pi^2)*E)/(12*(1-(mu^2))))*((t(i)/b(i))^2);
        if BMD(i) >= 0 % if top of cross section is under compression
            M_Buck(i) = fail(i) * I(i) / yfla(i);
        else
            M_Buck(i) = 0;
        end
    end
end

function [M_Buck] = MfailBuck2(CrossSectionProperties, E, mu, BMD) % case 2, ends of flange
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    I = CrossSectionProperties(:, 10);
    b = 10 * ones(1, size(BMD, 2)); % HARDCOODED
    yfla = CrossSectionProperties(:, 9);
    t = CrossSectionProperties(:, 4);

    fail = zeros(1, size(BMD, 2)); % initialize array
    for i = 1 : size(BMD, 2)
        fail(i) = ((0.425*(pi^2)*(E))/(12*(1-(mu^2))))*((t(i)./b(i))^2);
        if BMD(i) >= 0 % if top of cross section is under compression
            M_Buck(i) = fail(i) * I(i) / yfla(i);
        else
            M_Buck(i) = 0;
        end
    end
end

function [M_Buck] = MfailBuck3(CrossSectionProperties, E, mu, BMD) % case 3, webs
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 9) - CrossSectionProperties(:, 2);
    yweb = CrossSectionProperties(:, 9) - CrossSectionProperties(:, 2);
    t = CrossSectionProperties(:, 4);

    fail = zeros(1, size(BMD, 2)); % initialize array
    for i = 1 : size(BMD, 2)
        fail(i) = ((6*(pi^2)*(E))/(12*(1-(mu^2))))*((t(i)/b(i))^2);
        if BMD(i) >= 0 % if top of cross section is under compression
            M_Buck(i) = fail(i) * I(i) / yweb(i);
        else
            M_Buck(i) = 0;
        end
    end
end

function [M_Buck] = MfailBuck4 (CrossSectionProperties, E, mu, BMD) % case 3, webs
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 8) - CrossSectionProperties(:, 6); % y bot - bottom flange thickness
    t = CrossSectionProperties(:, 4);

    fail = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        fail(i) = ((6*(pi^2)*(E))/(12*(1-(mu^2))))*((t(i)/b(i))^2);
        if BMD(i) <= 0 % if bottom of cross section is under compression
            M_Buck(i) = -fail(i) * I(i) / b(i);
        else
            M_Buck(i) = 0;
        end
    end
end

function [M_Buck] = MfailBuck5(CrossSectionProperties, E, mu, BMD) % case 3, webs
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 5) - 2 * CrossSectionProperties(:, 4);
    ybot = CrossSectionProperties(:, 8);
    t = CrossSectionProperties(:, 6);

    fail = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        fail(i) = ((4*(pi^2)*(E))/(12*(1-(mu^2))))*((t(i)/b(i))^2);
        if BMD(i) <= 0 % if bottom of cross section is under compression
            M_Buck(i) = -fail(i) * I(i) / ybot(i);
        else
            M_Buck(i) = 0;
        end
    end
end

function [Fails] = GetFails(CrossSectionProperties, TauU, TauG, E, mu, SigT, SigC, BMD)
    % Adds all failure types into one array
    Fails = zeros(10, size(BMD, 2));
    Fails(1, :) = Vfail(CrossSectionProperties, TauU); % V_fail
    Fails(2, :) = VfailGlue(CrossSectionProperties, TauG); % V_failGlue
    Fails(3, :) = VfailBuck(CrossSectionProperties, E, mu); % V_buck
    Fails(4, :) = MfailMatT(CrossSectionProperties, SigT, BMD); % M_MatT
    Fails(5, :) = MfailMatC(CrossSectionProperties, SigC, BMD); % M_MatC
    Fails(6, :) = MfailBuck1(CrossSectionProperties, E, mu, BMD); % M_Buck1
    Fails(7, :) = MfailBuck2(CrossSectionProperties, E, mu, BMD); % M_Buck2
    Fails(8, :) = MfailBuck3(CrossSectionProperties, E, mu, BMD); % M_Buck3
    Fails(9, :) = MfailBuck4(CrossSectionProperties, E, mu, BMD); % M_Buck4
    Fails(10, :) = MfailBuck5(CrossSectionProperties, E, mu, BMD); % M_Buck5
end

function [Pfail] = FailLoad(SFD, BMD, Fails)
    mins = zeros(1, 10);
    
    for i = 1 : 3
        for j = 1 : size(SFD, 2)
            Fails(i, j) = Fails(i, j) / SFD(j);
        end
    end
    for i = 4 : 10
        for j = 1 : size(SFD, 2)
            Fails(i, j) = Fails(i, j) / BMD(j);
        end
    end

    Fails(Fails == 0) = 999999999; % i'm so cool

    for i = 1 : 10
        mins(i) = min(abs(Fails(i, :)));
    end
    
    mins;

    [Pfail, aaa] = min(mins);
end

function [Defls] = Deflections(x, BMD, I, E)
% Calculates deflections
% Input: I(1-D arrays), E (material property), BMD (1-D array)
% Outputx: Deflection for every value of x (1-D array) or for the midspan only
    CurvD = BMD ./ I ./ E;
    Thing = zeros(1, size(x, 2));
    for i = 2 : size(x, 2)
        Thing(i) = Thing(i - 1) + (CurvD(i - 1) * (x(i - 1) - 15)); % start changing BMD at (xP + 1) since moment at xP = 0 anyways relative to shear force
    end
    PlotDiagrams(x, size(x, 2) - 1, BMD, Thing)
    %Defls = (Thing(1075) / 2 - Thing(545)); % hardcode cause i'm cool

    IntegralCurvD = BMD ./ I ./ E;
    for i = 2 : size(BMD, 2)
        IntegralCurvD(i) = IntegralCurvD(i - 1) + (CurvD(i - 1)); % start changing BMD at (xP + 1) since moment at xP = 0 anyways relative to shear force
    end
    %PlotDiagrams(x, size(x, 2) - 1, BMD, IntegralCurvD);
    Defls = (Thing(1075) / IntegralCurvD(1075) / 2 - Thing(545) / IntegralCurvD(545)); % hardcode cause i'm cool
    %}
    % ^^^ NEED TO ACCOUTN FOR CENTROID OF AREAS (idk how to do this
    % might have to tdo thing - 15 cause it's offset from suport A
    Thing(1075);
    Thing(545);
end

function Plot2L(x, L, Pf, CrossSectionProperties, TauU, TauG, E, mu, SigT, SigC) % ADD LEGENDS
    SFD = zeros(1, size(x, 2));
    BMD = zeros(1, size(x, 2));
    [SFD, BMD] = ApplyTwoLoads(Pf, x, SFD, BMD);
    Fails = GetFails(CrossSectionProperties, TauU, TauG, E, mu, SigT, SigC, BMD);

    figure()
    subplot(2, 3, 1) % SFD
    plot(x, SFD, "k")
    xlim([0 L])
    title("Shear Force over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("V (N)")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 3, 2) % SFD w/ material and glue shear failure
    plot(x, SFD, "k")
    hold on
    plot(x, Fails(1, :), "r")
    plot(x, Fails(2, :), "g")
    plot(x, -Fails(1, :), "r")
    plot(x, -Fails(2, :), "g")
    hold off
    xlim([0 L])
    title("SFD vs. Material and Glue Shear Failures")
    xlabel("x (mm)")
    ylabel("V (N)")
    legend("SFD", "Material Shear", "Glue Shear")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 3, 3) % SFD w/ shear buckling failure
    plot(x, SFD, "k")
    hold on
    plot(x, Fails(3, :), "r")
    plot(x, -Fails(3, :), "r")
    hold off
    xlim([0 L])
    title("SFD vs. Shear Buckling Failure")
    xlabel("x (mm)")
    ylabel("V (N)")
    legend("SFD", "Shear Buck.")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 3, 4) % BMD
    plot(x, BMD, "k")
    xlim([0 L])
    title("Bending Moment over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("M (N mm)")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference
    
    subplot(2, 3, 5) % BMD w/ material moment failure
    plot(x, BMD, "k")
    hold on
    plot(x, Fails(4, :), "r")
    plot(x, Fails(5, :), "g")
    hold off
    xlim([0 L])
    title("BMD vs. Material Moment Failures")
    xlabel("x (mm)")
    ylabel("M (N mm)")
    legend("BMD", "Flange Moment", "Web Moment")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference

    subplot(2, 3, 6) % BMD w/ moment buckling failure (plate buckling)
    plot(x, BMD, "k")
    hold on
    plot(x, Fails(6, :), "r")
    plot(x, Fails(7, :), "g")
    plot(x, Fails(8, :), "b")
    plot(x, Fails(9, :), "y")
    plot(x, Fails(10, :), "m")
    hold off
    xlim([0 L])
    title("BMD vs. Plate Buckling Failures")
    xlabel("x (mm)")
    ylabel("M (N mm)")
    legend("SFD", "Mid Flange Buck.", "Side Flange Buck.", "Top Web Buck.", "Bottom Web Buck.", "Bottom Flange Buck.", 'Location', "southeast")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference

    set(gcf, 'Name', 'Failure Loads of Improved Design') % name of window
end

