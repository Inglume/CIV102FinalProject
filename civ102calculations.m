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

subplot(2, 1, 1) % SFD
plot(x, SFD_PL)
xlim([0 L])
subplot(2, 1, 2) % BMD
plot(x, BMD_PL)
xlim([0 L])

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
% have to find a way to deal with diaphragms

% format geometric inputs in some sort of way
% each column corresponds to different cross section
% tried implementing in another way, wehre each row corresponds to
% different cross section

AddCrossSection(0, 100, 2.54, 100, 1.27, 80, 1.27, 400, GeometricInputs)
AddCrossSection(550, 100, 2.54, 120, 1.27, 80, 1.27, 400, GeometricInputs)
% AddCrossSection(L, 100, 2.54, 100, 1.27, 80, 1.27, 400);

% Optional but you need to ensure that your geometric inputs are correctly implemented
% VisualizeBridge( {CrossSectionInputs} );
%% 3. Define Material Properties
SigT = 30; % tensile stress
SigC = 6; % compressive stress
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

%{

SectionalProperties = SectionProperties(GeometricInputs)
CrossSectionProperties = [GeometricInputs SectionalProperties] % geometric and sectional things combined (may want to do this differently)

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
    SFD(xP + 1) = SFD(xP + 1) + P;
    for i = xP + 2 : length(x)
        SFD(i) = SFD(i) + P;
        BMD(i) = BMD(i - 1) + (SFD(i) / 1000);
    end
end
function VisualizeBridge(GeometricInputs)
% Optional. Provides a graphical interpretation of user geometric inputs
end
function [SectionalProperties] = SectionProperties(GeometricInputs)
% Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
%   Input: Geometric Inputs. Format will depend on user
%   Output: Sectional Properties at every value of x. Each property is a 1-D array of length n
% idk why it says each property shoudl be legnth n, because properties
% should be consistent across a length with uniform cross section (unless
% diaphragms have to be accounted for
    for i = 1 : length(GeometricInputs) - 1 % get row # (for each cross section)
        for j = GeometricInputs(i, 0) : GeometricInputs(i, 1) % beginning of current cross section to beginning of next cross section
            ybar = 0;
            I = 0;
            Q = 0;
        end
    end
    SectionalProperties = [ybar I Q];
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

% xc: Location, x, of cross-section change
% bft: Top Flange Width
% tft: Top Flange Thickness
% hw: Web Height
% tw: Web Thickness (Assuming 2 separate webs)
% bfb: Bottom Flange Width
% tfb: Bottom Flange Thickness
% a: Diaphragm Spacing
function AddCrossSection(xc, bft, tft, hw, tw, bfb, tfb, a, GeometricInputs)
    GeometricInputs(end + 1, :) = [xc, bft, tft, hw, tw, bfb, tfb, a];
end

