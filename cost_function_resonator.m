function cost = cost_function_resonator(particle, targetMagdB, hMWS, filePrefix)
% COST_FUNCTION_RESONATOR
%
% Updates CST model parameters, runs simulations for two material states,
% extracts the S21 resonance, exports Touchstone files, and computes the PSO
% cost.
%
% The cost penalises:
%   1. weak resonance depth,
%   2. small resonant-frequency shift between material states,
%   3. broad 3 dB bandwidth.

    disp("Current position:");
    disp(particle.pos);

    disp("Current velocity:");
    disp(particle.vel);

    % CST material states used to evaluate resonance shift.
    materialState1 = "material_state_1";
    materialState2 = "material_state_2";
    
    % Labels used when exporting Touchstone files.
    touchstoneLabel1 = "state_1";
    touchstoneLabel2 = "state_2";

    %% UPDATE CST PARAMETERS
    % Map the particle position vector to CST design parameters.
    % The CST project must contain parameters named t, g, D, and r.

    invoke(hMWS, "StoreParameter", "t", particle.pos(1));
    invoke(hMWS, "StoreParameter", "g", particle.pos(2));
    invoke(hMWS, "StoreParameter", "D", particle.pos(3));
    invoke(hMWS, "StoreParameter", "r", particle.pos(4));

    % Rebuild the CST model after changing the geometric parameters. If the
    % rebuild fails, the particle is rejected by assigning infinite cost.
    if ~invoke(hMWS, "RebuildOnParametricChange", "true", "true")
        disp("Parameter update failed.");
        cost = inf;
        return;
    end

    %% FIRST MATERIAL STATE
    % First material state: used to evaluate resonance depth and bandwidth.
    % The exported Touchstone file allows the simulated response to be reviewed
    % after the optimisation has finished. (MUT: Material under test)
    invoke(hMWS, "DeleteResults");

    solidObject = invoke(hMWS, "Solid");
    invoke(solidObject, "ChangeMaterial", "CSRR:MUT", materialState1);

    hSolver = invoke(hMWS, "Solver");

    if ~invoke(hSolver, "Start")
        disp("Simulation problem.");
        cost = inf;
        return;
    end

    [fRes1, magActualdB, bandwidth3dB] = extract_resonance_from_cst(hMWS);

    disp("First resonant frequency and 3 dB bandwidth:");
    disp(fRes1);
    disp(bandwidth3dB);

    export_touchstone(hMWS, strcat(filePrefix, "_", touchstoneLabel1));

    %% SECOND MATERIAL STATE
    % Second material state: used to calculate the resonance shift between two
    % dielectric loading conditions.
    invoke(hMWS, "DeleteResults");
    invoke(solidObject, "ChangeMaterial", "CSRR:MUT", materialState2);

    if ~invoke(hSolver, "Start")
        disp("Simulation problem.");
        cost = inf;
        return;
    end

    [fRes2, ~, ~] = extract_resonance_from_cst(hMWS);

    disp("Second resonant frequency:");
    disp(fRes2);

    resonanceShift = abs(fRes1 - fRes2);

    disp("Resonance shift:");
    disp(resonanceShift);

    export_touchstone(hMWS, strcat(filePrefix, "_", touchstoneLabel2));

    %% COST CALCULATION
    % Two minima that are too far apart are rejected because the
    % minima may correspond to different resonant modes rather than a shifted
    % version of the same mode.
    if resonanceShift < 0.4
        cost = (magActualdB - targetMagdB) ...
            + 2 / resonanceShift ...
            + (10 * bandwidth3dB)^1.8;

        disp("Cost:");
        disp(cost);
    else
        cost = inf;
        disp("Resonant minima are too far apart and may not correspond to the same mode.");
        disp("Cost:");
        disp(cost);
    end
end


function [fRes, magMindB, bandwidth3dB] = extract_resonance_from_cst(hMWS)
% Extract resonant frequency, minimum S21 magnitude, and 3 dB bandwidth.

    hResultTree = invoke(hMWS, "ResultTree");
    resultFile = invoke( ...
        hResultTree, ...
        "GetFileFromTreeItem", ...
        "1D Results\S-Parameters\S2,1" ...
    );

    hResultComplex = invoke(hMWS, "Result1DComplex", "");
    invoke(hResultComplex, "Load", resultFile);

    hResult = invoke(hResultComplex, "Magnitude");

    minIndex = invoke(hResult, "GetGlobalMinimum");
    fRes = invoke(hResult, "GetX", minIndex);
    magMindB = mag2db(invoke(hResult, "GetY", minIndex));

    numPoints = invoke(hResult, "GetN");
    frequencies = zeros(1, numPoints);
    magnitudes = zeros(1, numPoints);

    for i = 1:numPoints
        frequencies(i) = invoke(hResult, "GetX", i - 1);
        magnitudes(i) = invoke(hResult, "GetY", i - 1);
    end

    magnitudesdB = mag2db(magnitudes);
    mag3dB = magMindB + 3;
    bandwidthIndices = find(magnitudesdB <= mag3dB);

    bandwidth3dB = ...
        frequencies(max(bandwidthIndices)) - frequencies(min(bandwidthIndices));
end


function export_touchstone(hMWS, filePrefix)
% Export full-range S-parameter results as a Touchstone file.

    touchstone = invoke(hMWS, "TOUCHSTONE");

    invoke(touchstone, "Reset");
    invoke(touchstone, "FileName", filePrefix);
    invoke(touchstone, "Impedance", "50");
    invoke(touchstone, "Format", "db");
    invoke(touchstone, "FrequencyRange", "full");
    invoke(touchstone, "Renormalize", "True");
    invoke(touchstone, "UseARResults", "False");
    invoke(touchstone, "SetNSamples", "1001");
    invoke(touchstone, "Write");
end
