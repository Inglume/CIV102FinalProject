%% 0. Initialize Parameters
L = 1310; % Length of bridge
n = L + 1; % Number of locations to evaluate bridge failure
x = linspace(0, L, n); % Define x coordinate
SFD = zeros(1, n); % Initialize SFD(x)
BMD = zeros(1, n); % Initialize SFD(x)

%% 1. Point Loading Analysis (SFD, BMD)

P = 1; % CHANGE THIS probapbly or do some loop
SFD2L = zeros(1, n);
BMD2L = zeros(1, n);
[SFD2L, BMD2L] = ApplyTwoLoads(1, x, SFD, BMD);
PlotDiagrams(x, L, SFD2L, BMD2L)


[SFDTrain, BMDTrain] = ApplyTrainLoad(x, SFD, BMD);
%PlotTrain(x, L, SFDTrain, BMDTrain)


%% 2. Define cross-sections

%xc = [0 550 L]; % Location, x, of cross-section change
%bft = [100 100 100]; % Top Flange Width
%tft = [2.54 2.54 2.54]; % Top Flange Thickness
%hw = [100 120 100]; % Web Height
%tw = [1.27 1.27 1.27]; % Web Thickness (Assuming 2 separate webs)
%bfb = [80 80 80]; % Bottom Flange Width
%tfb = [1.27 1.27 1.27]; % Bottom Flange Thickness
%a = [400 400 400]; % Diaphragm Spacing
% have to find a way to deal with diaphragms (or not)

%GeometricInputs(end + 1, :) = [0, 100, 2.54, 100, 1.27, 80, 1.27, 400]; % add something for length of glue tabs on top and bottom (thickness should be same as thickness of webs
%GeometricInputs(end + 1, :) = [550, 100, 2.54, 120, 1.27, 80, 1.27, 400];
%GeometricInputs(end + 1, :) = [L, 100, 2.54, 100, 1.27, 80, 1.27, 400];

% Design 0
GeometricInputs = [];

GeometricInputs(end + 1, :) = [0, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [30, 100, 1.27, 72.46, 1.27, 80, 1.27, 520];
GeometricInputs(end + 1, :) = [550, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [580, 100, 1.27, 72.46, 1.27, 80, 1.27, 480];
GeometricInputs(end + 1, :) = [1060, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [1090, 100, 1.27, 72.46, 1.27, 80, 1.27, 160];
GeometricInputs(end + 1, :) = [1280, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [L, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];

% Design 1.0 woowwowow

% Optional but you need to ensure that your geometric inputs are correctly implemented
% VisualizeBridge( {CrossSectionInputs} );
%% 3. Define Material Properties
SigT = 30; % tensile stress
SigC = 6; % compressive stress
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

CrossSectionProperties = SectionProperties(GeometricInputs, n);

%{
for i = 1 : size(GeometricInputs, 1) - 1
    cs = CrossSectionProperties(GeometricInputs(i, 1) + 1, :);
    sprintf("Cross Section @ %d mm - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Q: %.3g", GeometricInputs(i, 1), cs(8:11))
end
%}

%cs1 = CrossSectionProperties(1, :);
%cs2 = CrossSectionProperties(551, :);

%sprintf("Cross Section 1 - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Q: %d", cs1(8:11))
%sprintf("Cross Section 2 - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Q: %d", cs2(8:11))



%% 4. Calculate Failure Moments and Shear Forces
V_Mat = Vfail(CrossSectionProperties, TauU);
V_Mat(1)
V_Glue = VfailGlue(CrossSectionProperties, TauG);
V_Glue(1)



%{

V_Buck = VfailBuck(CrossSectionProperties, E, mu );
M_MatT = MfailMatT(CrossSectionProperties, SigT);
M_MatC = MfailMatC(CrossSectionProperties, SigC);
M_Buck1 = MfailBuck(CrossSectionProperties, E, mu, 1 );
M_Buck2 = MfailBuck(CrossSectionProperties, E, mu, 2 );
M_Buck3 = MfailBuck(CrossSectionProperties, E, mu, 3 );


%% 4.7 Calculate Failure Load
Pf = FailLoad(P, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3);

%% Visualization
VisualizePL(x, P, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, Pf);

%% 5. Curvature, Slope, Deflections
Defls = Deflections(x, BMD_PL, I, E);
%}
%% Functions

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
    
    [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xP, -P, x, SFD, BMD);
end

function [SFD, BMD] = ApplyTwoLoads(P, x, SFD, BMD) % P is force of each individual load
    xA = 15; % location of support A
    xB = xA + 1060; % location of support B
    xP1 = xA + 550;
    xP2 = xB + 190;
    By = (P * (xP1 - xA) + P * (xP2 - xA)) / (xB - xA); % support B located at 1050mm from support A, support A located at x = 0
    Ay = 2 * P - By; % force equilibrium on x
    
    [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xP1, -P, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xP2, -P, x, SFD, BMD);
end

function [trainSFD, trainBMD] = ApplyTrainLoad(x, SFD, BMD) % get rid of SFD, BMD here
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
    position = [15, 137]; % two cases of train loading
    
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
        BMD(i) = BMD(i - 1) + (SFD(i - 1) / 1000); % start changing BMD at (xP + 1) since moment at xP = 0 anyways relative to shear force
    end
end

function PlotDiagrams(x, L, SFD, BMD) % include curvature diagram
    subplot(2, 1, 1) % SFD
    plot(x, SFD)
    xlim([0 L])
    title("Shear Force over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("V (kN)")
    ax = gca;
    ax.XAxisLocation = 'origin';

    subplot(2, 1, 2) % BMD
    plot(x, BMD)
    xlim([0 L])
    title("Moment over Horizontal Distance")
    xlabel("x (mm)")
    ylabel("M (kN m)")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference
    
    set(gcf, 'Name', 'Force and Moment Diagrams') % name of window
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
    ylabel("V (kN)")
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
    ylabel("M (kN m)")
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(ax, 'YDir','reverse') % may not want it reversed, personal preference
    
    set(gcf, 'Name', 'Force and Moment Diagrams') % name of window
end

function VisualizeBridge(GeometricInputs)
% Optional. Provides a graphical interpretation of user geometric inputs
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
        ybot = 0; % sum of area times relative y divided by sum of area
        ytop = 0; % ybar from top
        I = 0; % sum of individual i plus area times distance squared
        Qmax = 0; % sum of area times distance
        yfla = 0; % distance from global centroid to centroid of flange above ybar
        yweb = 0; % distance from global centroid to centroid of web above ybar
        ytab = 0; % distance from global centroid to centroid of tab above ybar
        yfgl = 0; % distance from glue to flange
        Qglue = 0;% implement this stuf
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
    V_fail = I .* b ./ Qcent * TauU;
end

function [V_failGlue] = VfailGlue(CrossSectionProperties, TauG)
% Calculates shear forces at every value of x that would cause a glue shear failure
%   Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
%   Output: V_fail a 1-D array of length n
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 4) * 2 + 20;
    Qglue = CrossSectionProperties(:, 16);
    V_failGlue = I .* b ./ Qglue * TauG;
end

function [V_Buck] = VfailBuck(CrossSectionProperties, E, mu)  
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property)
    I = CrossSectionProperties(:, 10);
    b = CrossSectionProperties(:, 4) * 2;
    Qcent = CrossSectionProperties(:, 11);
    h = CrossSectionProperties(:, 2) + CrossSectionProperties(:, 3) + CrossSectionProperties(:, 6);
    a = CrossSectionProperties(:, 7);
    TauCrit = ((5 .* (pi^2) .* E) ./ (12 .* (1 - (mu^2)))) .* ((1.27 ./ h) .^ 2 + (1.27 ./ a) .^ 2);
    V_Buck = TauCrit .* I .* b ./ Qcent;
end

function [M_MatT] = MfailMatT(CrossSectionProperties, SigT, BMD)
% Calculates bending moments at every value of x that would cause a matboard tension failure
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array)
% Output: M_MatT a 1-D array of length n
%[I, ybot, ytop] = SectionalProperties;
    M_MatT = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom
            M_MatT(i) = SigT * I(i) / ybot(i);
        elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top
            M_MatT(i) = -SigT * I(i) / ytop(i);
        end
    end
end

function [M_MatC] = MfailMatC(CrossSectionProperties, SigC, BMD) % Similar to MfailMatT
    M_MatC = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        if BMD(i) > 0 % If the moment is positive, the compression failure will be at the top
            M_MatC(i) = SigC * I(i) / ytop(i);
        elseif BMD(i) < 0 % If the moment is negative, the compression failure will be at the bottom
            M_MatC(i) = -SigC * I(i) / ybot(i);
        end
    end
end

function [M_Buck] = MfailBuck1(CrossSectionProperties, E, mu, BMD) % case 1, middle part of flange
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    fail = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        fail(i) = ((4*(pi^2)*(E))/(12*(1-(mu^2))))*((t/b(i))^2);
        if yplate > 0
            if BMD(i) > 0
                M_Buck(i) = fail(i) * I(i) / yplate;
            else
                M_Buck(i) = 0;
            end
        elseif yplate < 0
            if BMD(i) < 0
                M_Buck(i) = fail(i) * I(i) / -yplate;
            else
                M_Buck(i) = 0;
            end
        end
    end
end

function [M_Buck] = MfailBuck2(CrossSectionProperties, E, mu, BMD) % case 2, ends of flange
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    fail = zeros(1, size(BMD, 2));    
    for i = 1 : size(BMD, 2)
        fail(i) = ((0.425*(pi^2)*(E))/(12*(1-(mu^2))))*((t/b(i))^2)
        if yplate > 0
            if BMD(i) > 0
                M_Buck(i) = fail(i) * I(i) / yplate
            else
                M_Buck(i) = 0
            end
        elseif yplate < 0
            if BMD(i) < 0
                M_Buck(i) = fail(i) * I(i) / -yplate
            else
                M_Buck(i) = 0
            end
        end
    end
end

function [M_Buck] = MfailBuck3(CrossSectionProperties, E, mu, BMD) % case 3, webs
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    fail = zeros(1, size(BMD, 2));
    for i = 1 : size(BMD, 2)
        fail(i) = ((6*(pi^2)*(E))/(12*(1-(mu^2))))*((t/b(i))^2);
        if yplate > 0
            if BMD(i) > 0
                M_Buck(i) = fail(i) * I(i) / yplate
            else
                M_Buck(i) = 0
            end
        else
            if BMD(i) < 0
                M_Buck(i) = fail(i) * I(i) / -yplate
            else
                M_Buck(i) = 0

            end
        end
    end
end
%}

%{
function VisualizePL(x, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, Pf)
% Plots all outputs of design process
end

function [Defls] = Deflections(x, BMD, I, E)
% Calculates deflections
% Input: I(1-D arrays), E (material property), BMD (1-D array)
% Output: Deflection for every value of x (1-D array) or for the midspan only
    Defls = 0;
end
%}