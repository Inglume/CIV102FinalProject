%% 0. Initialize Parameters
L = 1250; % Length of bridge
n = L + 1; % Number of locations to evaluate bridge failure
x = linspace(0, L, n); % Define x coordinate
SFD_PL = zeros(1, n); % Initialize SFD(x)
BMD_PL = zeros(1, n); % Initialize SFD(x)

%% 1. Point Loading Analysis (SFD, BMD)
P = 318;
[SFD_PL, BMD_PL] = ApplyPL(550, P, x, SFD_PL, BMD_PL); % Construct SFD, BMD
[SFD_PL, BMD_PL] = ApplyPL(L, P, x, SFD_PL, BMD_PL); % Construct SFD, BMD

%PlotDiagrams(x, L, SFD_PL, BMD_PL)

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

GeometricInputs(end + 1, :) = [0, 100, 2.54, 100, 1.27, 80, 1.27, 400];
GeometricInputs(end + 1, :) = [550, 100, 2.54, 120, 1.27, 80, 1.27, 400];
GeometricInputs(end + 1, :) = [L, 100, 2.54, 100, 1.27, 80, 1.27, 400];



% Optional but you need to ensure that your geometric inputs are correctly implemented
% VisualizeBridge( {CrossSectionInputs} );
%% 3. Define Material Properties
SigT = 30; % tensile stress
SigC = 6; % compressive stress
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

CrossSectionProperties = SectionProperties(GeometricInputs, n)

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
function [SFD, BMD] = ApplyPL(xP, P, x, SFD, BMD)
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
%   Input: location and magnitude of point load. The previous SFD can be entered as input to
% construct SFD of multiple point loads
%   Output: SFD, BMD both 1-D arrays of length n
% need to account for support reaction forces too
    xA = 0; % location of support A
    xB = 1050; % location of support B
    By = P * (xP - xA) / (xB - xA); % support B located at 1050mm from support A, support A located at x = 0
    Ay = P - By; % force equilibrium on x
    
    [SFD, BMD] = UpdateDiagrams(xA, Ay, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xB, By, x, SFD, BMD);
    [SFD, BMD] = UpdateDiagrams(xP, -P, x, SFD, BMD);
end

function [SFD, BMD] = UpdateDiagrams(xP, P, x, SFD, BMD)
    SFD(xP + 1) = SFD(xP + 1) + P;  % x location must be increased by one since matrices are 1-indexed
    
    for i = xP + 2 : length(x)
        SFD(i) = SFD(i) + P;
        BMD(i) = BMD(i - 1) + (SFD(i) / 1000); % start changing BMD at (xP + 1) since moment at xP = 0 anyways relative to shear force
    end
end


% xc: Location, x, of cross-section change
% bft: Top Flange Width
% tft: Top Flange Thickness
% hw: Web Height
% tw: Web Thickness (Assuming 2 separate webs)
% bfb: Bottom Flange Width
% tfb: Bottom Flange Thickness
% a: Diaphragm Spacing
function AddCrossSection(xc, bft, tft, hw, tw, bfb, tfb, a, GeometricInputs) % don't need htis shite
    GeometricInputs(end + 1, :) = [xc, bft, tft, hw, tw, bfb, tfb, a];
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

function TrainLoadDiagram(x, L, SFD, BMD)
    print("a");
end

function VisualizeBridge(GeometricInputs)
% Optional. Provides a graphical interpretation of user geometric inputs
end

% check this function
function CrossSectionProperties = SectionProperties(GeometricInputs, n) % include y plate 1, 2, 3, 4, 5, 6
% Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
%   Input: Geometric Inputs. Format will depend on user
%   Output: Sectional Properties at every value of x. Each property is a 1-D array of length n
% idk why it says each property shoudl be legnth n, because properties
% should be consistent across a length with uniform cross section (unless
% diaphragms have to be accounted for
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
        Q = 0; % sum of area times distance
        for k = 1 : 3
                ybot = ybot + areas(k) * distances(k);
        end
        ybot = ybot / sum(areas);
        ytop = GeometricInputs(3) + GeometricInputs(4) + GeometricInputs(7) - ybot;
        for k = 1 : 3
                I = I + secondMomInert(k) + areas(k) * (distances(k) - ybot) ^ 2; % prbaobly somethign wrong with this calc
                Q = Q + areas(k) * (distances(k) - ybot);
        end
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