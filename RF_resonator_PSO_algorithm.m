% PARTICLE SWARM OPTIMISATION COUPLED TO CST
%
% This script implements a Particle Swarm Optimisation (PSO) loop in MATLAB
% and uses CST Studio Suite through its COM interface to evaluate each design.
%
% The CST project should contain the design parameters listed below:
%   t, g, D, r
%
% The associated cost function updates these parameters in CST, runs the
% solver, extracts the S21 response, and evaluates the resonator performance.

clear; clc;

%% USER SETTINGS

projectFile = fullfile("cst_models", "resonator_template.cst");
touchstoneOutputFolder = fullfile("results", "touchstone_exports");

if ~exist(touchstoneOutputFolder, "dir")
    mkdir(touchstoneOutputFolder);
end

%% SET UP CONNECTION WITH CST

hCST = actxserver("CSTStudio.Application");
hMWS = invoke(hCST, "OpenFile", projectFile);

%% PSO PARAMETERS

targetMagdB = -30;      % target resonant magnitude in dB
maxIter = 10;           % number of PSO iterations
nParticles = 10;        % number of particles
nVar = 4;               % number of design variables
varSize = [1, nVar];

inertiaWeight = 1;      % particle inertia
damping = 0.9;          % damping coefficient
c1 = 1.4;               % cognitive acceleration coefficient
c2 = 1.4;               % social acceleration coefficient

%% DESIGN VARIABLE BOUNDS

% Design variables correspond to CST parameters:
%   pos(1) -> t
%   pos(2) -> g
%   pos(3) -> D
%   pos(4) -> r

[t_lb, t_ub, g_lb, g_ub, D_lb, D_ub, r_lb, r_ub] = ...
    deal(1, 2, 0.3, 0.8, 1, 4, 5, 8);

lowerBounds = [t_lb, g_lb, D_lb, r_lb];
upperBounds = [t_ub, g_ub, D_ub, r_ub];

velocityMax = 0.5 * (upperBounds - lowerBounds);
velocityMin = -velocityMax;

%% PSO INITIALISATION

emptyParticle.pos = [];
emptyParticle.vel = [];
emptyParticle.cost = [];
emptyParticle.best.pos = [];
emptyParticle.best.cost = [];

particles = repmat(emptyParticle, nParticles, 1);
globalBest.cost = inf;  % starting from worst value

for i = 1:nParticles

    filePrefix = fullfile( ...
        touchstoneOutputFolder, ...
        strcat("iter0_particle", num2str(i)) ...
    );

    particles(i).pos = lowerBounds + rand(varSize) .* (upperBounds - lowerBounds);
    particles(i).vel = zeros(1, nVar);

    particles(i).cost = cost_function_resonator( ...
        particles(i), targetMagdB, hMWS, filePrefix ...
    );

    particles(i).best.pos = particles(i).pos;
    particles(i).best.cost = particles(i).cost;

    if particles(i).best.cost < globalBest.cost
        globalBest = particles(i).best;
    end
end

tGlobalBest = zeros(1, maxIter);
gGlobalBest = zeros(1, maxIter);
DGlobalBest = zeros(1, maxIter);
rGlobalBest = zeros(1, maxIter);
costHistory = zeros(1, maxIter);

%% PSO MAIN LOOP
% Each PSO iteration updates the full swarm once. For every particle, the
% geometry parameters are updated, the CST model is simulated twice - once
% for each material under test, and the cost is evaluated from the
% resulting S-parameter response for each material and the differences between them.

for iter = 1:maxIter

    for particleIndex = 1:nParticles

        filePrefix = fullfile( ...
            touchstoneOutputFolder, ...
            strcat("iter", num2str(iter), "_particle", num2str(particleIndex)) ...
        );

        % Velocity update combines inertia, attraction to the particle's personal
        % best position, and attraction to the best position found by the swarm.
        particles(particleIndex).vel = ...
            inertiaWeight * particles(particleIndex).vel ...
            + c1 * rand(varSize) .* ...
              (particles(particleIndex).best.pos - particles(particleIndex).pos) ...
            + c2 * rand(varSize) .* ...
              (globalBest.pos - particles(particleIndex).pos);

        % Apply velocity limits
        particles(particleIndex).vel = max(particles(particleIndex).vel, velocityMin);
        particles(particleIndex).vel = min(particles(particleIndex).vel, velocityMax);

        % Update position
        particles(particleIndex).pos = ...
            particles(particleIndex).pos + particles(particleIndex).vel;

        % Clamp the updated particle position so that all CST geometry parameters
        % remain within physically meaningful design bounds.
        particles(particleIndex).pos = max(particles(particleIndex).pos, lowerBounds);
        particles(particleIndex).pos = min(particles(particleIndex).pos, upperBounds);

        % Evaluate updated particle
        particles(particleIndex).cost = cost_function_resonator( ...
            particles(particleIndex), targetMagdB, hMWS, filePrefix ...
        );

        % Update personal best
        if particles(particleIndex).cost < particles(particleIndex).best.cost
            particles(particleIndex).best.cost = particles(particleIndex).cost;
            particles(particleIndex).best.pos = particles(particleIndex).pos;

            % Update global best
            if particles(particleIndex).best.cost < globalBest.cost
                globalBest = particles(particleIndex).best;
            end
        end
    end

    tGlobalBest(iter) = globalBest.pos(1);
    gGlobalBest(iter) = globalBest.pos(2);
    DGlobalBest(iter) = globalBest.pos(3);
    rGlobalBest(iter) = globalBest.pos(4);
    costHistory(iter) = globalBest.cost;
    
    % Reduce inertia over time to shift the swarm gradually from exploration
    % toward exploitation.
    inertiaWeight = inertiaWeight * damping;
end

disp("Global best solution:");
disp(globalBest);

%% SAVE OPTIMISATION HISTORY

historyTable = table( ...
    (1:maxIter)', ...
    tGlobalBest(:), ...
    gGlobalBest(:), ...
    DGlobalBest(:), ...
    rGlobalBest(:), ...
    costHistory(:), ...
    'VariableNames', { ...
        'Iteration', 't', 'g', 'D', 'r', 'Cost' ...
    } ...
);

writetable(historyTable, fullfile("results", "pso_history.csv"));

%% PLOT RESULTS

figure;
plot(tGlobalBest);
title("t progression");
xlabel("Iteration");
ylabel("t");
grid on;

figure;
plot(gGlobalBest);
title("g progression");
xlabel("Iteration");
ylabel("g");
grid on;

figure;
plot(DGlobalBest);
title("D progression");
xlabel("Iteration");
ylabel("D");
grid on;

figure;
plot(rGlobalBest);
title("r progression");
xlabel("Iteration");
ylabel("r");
grid on;

figure;
plot(costHistory);
title("Cost progression");
xlabel("Iteration");
ylabel("Cost");
grid on;
