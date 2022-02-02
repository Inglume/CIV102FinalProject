%% 0. Initialize Parameters
L = 1310; % Length of bridge
n = L + 1; % Number of locations to evaluate bridge failure
x = linspace(0, L, n); % Define x coordinate
SFD = zeros(1, n); % Initialize SFD(x)
BMD = zeros(1, n); % Initialize SFD(x)

%% 1. Point Loading Analysis (SFD, BMD)

P = 1; % CHANGE THIS probapbly or do some loop
SFD = zeros(1, n);
BMD = zeros(1, n);
%[SFD2L, BMD2L] = ApplyTwoLoads(1, x, SFD, BMD);
%PlotDiagrams(x, L, SFD_2L, BMD_2L)

trainLength = 856;
t = zeros(0, trainLength);
%trainSFD = zeros(1, trainLength + 1);
%trainBMD = zeros(1, trainLength + 1);
[SFDTrain, BMDTrain] = ApplyTrainLoad(x, SFD, BMD);

%PlotDiagrams(x, L, SFDTrain, BMDTrain)
PlotTrain(x, L, SFDTrain, BMDTrain)

%% 2. Define cross-sections
% There are many (more elegant ways) to construct cross-section objects

GeometricInputs = [];

xc = [0 550 L]; % Location, x, of cross-section change
bft = [100 100 100]; % Top Flange Width
tft = [2.54 2.54 2.54]; % Top Flange Thickness
hw = [100 120 100]; % Web Height
tw = [1.27 1.27 1.27]; % Web Thickness (Assuming 2 separate webs)
bfb = [80 80 80]; % Bottom Flange Width
tfb = [1.27 1.27 1.27]; % Bottom Flange Thickness
a = [400 400 400]; % Diaphragm Spacing
% have to find a way to deal with diaphragms (or not)

% format geometric inputs in some sort of way
% each column corresponds to different cross section
% tried implementing in another way, wehre each row corresponds to
% different cross section

%GeometricInputs(end + 1, :) = [0, 100, 2.54, 100, 1.27, 80, 1.27, 400]; % add something for length of glue tabs on top and bottom (thickness should be same as thickness of webs
%GeometricInputs(end + 1, :) = [550, 100, 2.54, 120, 1.27, 80, 1.27, 400];
%GeometricInputs(end + 1, :) = [L, 100, 2.54, 100, 1.27, 80, 1.27, 400];

% Design 0
GeometricInputs(end + 1, :) = [0, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [30, 100, 1.27, 72.46, 1.27, 80, 1.27, 520];
GeometricInputs(end + 1, :) = [550, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [580, 100, 1.27, 72.46, 1.27, 80, 1.27, 480];
GeometricInputs(end + 1, :) = [1060, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [1090, 100, 1.27, 72.46, 1.27, 80, 1.27, 160];
GeometricInputs(end + 1, :) = [1280, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];
GeometricInputs(end + 1, :) = [L, 100, 1.27, 72.46, 1.27, 80, 1.27, 30];

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


% find a way to automate this for ease of use
%cs1 = CrossSectionProperties(1, :);
%cs2 = CrossSectionProperties(551, :);

%sprintf("Cross Section 1 - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Q: %d", cs1(8:11))
%sprintf("Cross Section 2 - ybot: %.3g mm ytop: %.3g mm I: %.3g mm^4 Q: %d", cs2(8:11))


%{
%% 4. Calculate Failure Moments and Shear Forces
V_Mat = Vfail(CrossSectionInputs, TauU);
V_Glue = VfailGlue(CrossSectionInputs, TauG);
V_Buck = VfailBuck(CrossSectionInputs, E, mu );
M_MatT = MfailMatT(CrossSectionInputs, SigT);
M_MatC = MfailMatC(CrossSectionInputs, SigC);
M_Buck1 = MfailBuck(CrossSectionInputs, E, mu, 1 );
M_Buck2 = MfailBuck(CrossSectionInputs, E, mu, 2 );
M_Buck3 = MfailBuck(CrossSectionInputs, E, mu, 3 );

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

function [trainSFD, trainBMD] = ApplyTrainLoad(x, SFD, BMD)
% Constructs SFD and BMD from application of train Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous SFD can be entered as input to
% construct SFD of multiple point loads
%   Output: SFD, BMD both 2-D arrays of length n and height i have no idea
% need to account for support reaction forces too
    P = 400;
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
    %trainSFD = SFD;
    %trainBMD = BMD;

    %PlotTrain(x, trainSFD, trainBMD)
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

function PlotTrain2 (x, L, SFD, BMD)
    trainLength = 856;
    t = tiledlayout(2,1);
    ax1 = nexttile;
    xlim([0 L])
    xlabel(ax1, "x (mm)")
    ylabel(ax1, "V (kN)")
    title(ax1, "Shear Force Diagram")
    ax1.XAxisLocation = 'origin';
    hold(ax1, 'on')
    ax2 = nexttile;
    xlim([0 L])
    xlabel(ax2, "x (mm)")
    ylabel(ax2, "M (kN m)")
    title(ax2, "Bending Moment Diagram")
    ax2.XAxisLocation = 'origin';
    set(ax2, 'YDir','reverse')
    hold(ax2, 'off')

    for i = 1 : L - size(SFD, 2)
        plot(ax1, x(i : trainLength + i + 1), SFD)
        plot(ax2, x(i : trainLength + i + 1), BMD)
    end

    title(t, "Diagrams for Train Loading Scenario")
end

function VisualizeBridge(GeometricInputs)
% Optional. Provides a graphical interpretation of user geometric inputs
end

% check this function
function CrossSectionProperties = SectionProperties(GeometricInputs, n) % include y plate 1, 2, 3, 4, 5, 6
% Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
%   Input: Geometric Inputs. Format will depend on user
%   Output: Sectional Properties at every value of x. Each property is a 1-D array of length n
    sect = zeros(n, 4); % columns: ybot, ytop, I, Q
    geom = zeros(n, size(GeometricInputs, 2) - 1);
    for i = 1 : (size(GeometricInputs, 1) - 1) % get row # (for each cross section)
        areas = zeros(3, 1); % length three (top, web, bottom)
        distances = zeros(3, 1); % length three (distance from bot)
        secondMomInert = zeros(3, 1); % second moment of inertias for each part
        areas(1) = GeometricInputs(i, 2) * GeometricInputs(i, 3);
        areas(2) = GeometricInputs(i, 5) * GeometricInputs(i, 4);
        areas(3) = GeometricInputs(i, 6) * GeometricInputs(i, 7);
        distances(3) = GeometricInputs(i, 7) / 2;
        distances(2) = GeometricInputs(i, 4) / 2 + distances(3);
        distances(1) = GeometricInputs(i, 3) / 2 + distances(2);
        secondMomInert(1) = (GeometricInputs(i, 2) * GeometricInputs(i, 3) ^ 2) / 12;
        secondMomInert(2) = (GeometricInputs(i, 4) * GeometricInputs(i, 5) ^ 2) / 12;
        secondMomInert(3) = (GeometricInputs(i, 6) * GeometricInputs(i, 7) ^ 2) / 12;
        ybot = 0; % sum of area times relative y divided by sum of area
        ytop = 0; % ybar from top
        I = 0; % sum of individual i plus area times distance squared
        yfla = 0; % distance from global centroid to centroid of flange above ybar
        yweb = 0; % distance from global centroid to centroid of web above ybar
        Q = 0; % sum of area times distance
        for k = 1 : 3
                ybot = ybot + areas(k) * distances(k);
        end
        ybot = ybot / sum(areas);
        ytop = GeometricInputs(i, 3) + GeometricInputs(i, 4) + GeometricInputs(i, 7) - ybot;
        for k = 1 : 3
                I = I + secondMomInert(k) + areas(k) * (distances(k) - ybot) ^ 2; % prbaobly somethign wrong with this calc
                % Q = Q + areas(k) * (distances(k) - ybot); % do Q just for areas above centroid
        end
        yfla = distances(1) - ybot;
        yweb = (GeometricInputs(i, 4) - ybot + GeometricInputs(i, 7)) / 2;
        Q = areas(1) * yfla;
        Q = Q + GeometricInputs(i, 5) * 2 * yweb ^ 2; % multiply thickness of web by height (yweb * 2), multiply by distance (yweb)
        for j = (GeometricInputs(i, 1) + 1) : GeometricInputs(i + 1, 1) + 1 % beginning of current cross section to beginning of next cross section
            geom(j, :) = GeometricInputs(i, 2 : end);
            sect(j, 1) = ybot;
            sect(j, 2) = ytop;
            sect(j, 3) = I;
            sect(j, 4) = Q;
        end
    end
    CrossSectionProperties = [geom sect];
end

function [V_fail] = Vfail(SectionalProperties, TauU)
% Calculates shear forces at every value of x that would cause a matboard shear failure
%   Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property)
%   Output: V_fail a 1-D array of length n
 I = {SectionalProperties};
 b = {SectionalProperties};
Qcent = {SectionalProperties};
 V_fail = TauU .* I .* b ./ Qcent;
end

function [V_Buck] = VfailBuck(SectionalProperties, E, mu)
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web
%   Input: Sectional Properties (list of 1-D arrays), E, mu (material property)
%   Output: V_Buck a 1-D array of length n
    V_Buck = 0;
end

function [M_MatT] = MfailMatT(SectionalProperties, SigT, BMD)
% Calculates bending moments at every value of x that would cause a matboard tension failure
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array)
% Output: M_MatT a 1-D array of length n
%[I, ybot, ytop] = SectionalProperties;
    for i = 1 : length(x)
        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom
            M_MatT(i) = SigT * I(i) / ybot(i);
        elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top
            M_MatT(i) = -SigT * I(i) / ytop(i);
        end
    end
end

function [M_MatT] = MfailMatC(SectionalProperties, SigC, BMD) % Similar to MfailMatT
    M_MatT = 0;
end

function [M_Buck] = MfailBuck(SectionalProperties, E, mu, BMD)
% Calculates bending moments at every value of x that would cause a buckling failure
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array)
% Output: M_MatBuck a 1-D array of length n
    M_Buck = 0;
end


function [Pf] = FailLoad(P, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3)
% Calculates the magnitude of the load P that will cause one of the failure mechanisms to occur
% Input: SFD, BMD under the currently applied points loads (P) (each 1-D array of length n)
% {V_Mat, V_Glue, … M_MatT, M_MatC, … } (each 1-D array of length n)
% Output: Failure Load value Pf
    Pf = 0;
end

function VisualizePL(x, SFD, BMD, V_Mat, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, Pf)
% Plots all outputs of design process
end

function [Defls] = Deflections(x, BMD, I, E)
% Calculates deflections
% Input: I(1-D arrays), E (material property), BMD (1-D array)
% Output: Deflection for every value of x (1-D array) or for the midspan only
    Defls = 0;
end