%% Resonator S-Parameter Averaging and Post-Processing
% Supporting analysis script for resonators designed/optimized using PSO + CST.
%
% The script reads repeated Touchstone (.s2p) files, extracts S21, averages
% repeated measurements for each group, identifies resonant frequency and
% minimum S21 magnitude, and plots selected experimental conditions.
%
% Validation measurements are processed in the same way as all other groups.
% The difference between "validation start" and "validation end" is then used
% as an additional uncertainty component.

clear; clc; close all;

%% USER SETTINGS

dataFolder = fullfile("data");
resultsFolder = fullfile("results");

% Include all groups that exist in the data folder, including validation checks.
groups = ["validation start", "baseline", "condition 1", "condition 2", "validation end"];

groupLabels = ["Validation start", "Baseline", "Condition 1", "Condition 2", "Validation end"];

% Only these groups are shown in the main plots.
plotGroupIndices = [2, 3, 4];

repeatLabels = [" m1", " m2", " m3"];

frequencyRangeGHz = [2.0, 2.8];

outputPort = 2;
inputPort = 1;

sgolayWindowSingleTrace = 1;
sgolayWindowAverage = 10;

plotTitle = "Microwave Resonator S21 Response";

% Manually specified baseline uncertainty components.
% Replace these with values from the relevant uncertainty budget.
baselineFrequencyUncertaintyGHz = 0.00372;
baselineMagnitudeUncertaintydB = 0.062;

if ~exist(resultsFolder, "dir")
    mkdir(resultsFolder);
end

%% READ FILES AND EXTRACT S21

numGroups = numel(groups);
numRepeats = numel(repeatLabels);
numFiles = numGroups * numRepeats;

s21Traces = [];
singleTraceResFreqGHz = zeros(1, numFiles);
singleTraceMinS21dB = zeros(1, numFiles);

fileIndex = 1;

for groupIndex = 1:numGroups
    for repeatIndex = 1:numRepeats

        fileName = groups(groupIndex) + repeatLabels(repeatIndex) + ".s2p";
        filePath = fullfile(dataFolder, fileName);

        if ~isfile(filePath)
            error("File not found: %s", filePath);
        end

        sParams = sparameters(filePath);

        if fileIndex == 1
            frequencyGHz = sParams.Frequencies / 1e9;
            frequencyMask = frequencyGHz > frequencyRangeGHz(1) & ...
                            frequencyGHz <= frequencyRangeGHz(2);
            selectedFrequencyGHz = frequencyGHz(frequencyMask);
            numFrequencyPoints = sum(frequencyMask);

            s21Traces = zeros(numFiles, numFrequencyPoints);
        end

        s21dB = 20 * log10(abs(rfparam(sParams, outputPort, inputPort)));
        s21dB = s21dB(frequencyMask);

        s21dBSmoothed = smoothdata(s21dB, "sgolay", sgolayWindowSingleTrace);

        s21Traces(fileIndex, :) = s21dB;

        [singleTraceMinS21dB(fileIndex), minIndex] = min(s21dBSmoothed);
        singleTraceResFreqGHz(fileIndex) = selectedFrequencyGHz(minIndex);

        fileIndex = fileIndex + 1;
    end
end

%% AVERAGE REPEATED MEASUREMENTS BY GROUP

averagedS21dB = zeros(numGroups, numFrequencyPoints);
averagedS21dBSmoothed = zeros(numGroups, numFrequencyPoints);

resonantFrequencyGHz = zeros(1, numGroups);
minimumS21dB = zeros(1, numGroups);

stdResonantFrequencyGHz = zeros(1, numGroups);
stdMinimumS21dB = zeros(1, numGroups);

for groupIndex = 1:numGroups

    rowIndices = (1:numRepeats) + (groupIndex - 1) * numRepeats;

    averagedS21dB(groupIndex, :) = mean(s21Traces(rowIndices, :), 1);

    averagedS21dBSmoothed(groupIndex, :) = smoothdata( ...
        averagedS21dB(groupIndex, :), ...
        "sgolay", ...
        sgolayWindowAverage ...
    );

    [minimumS21dB(groupIndex), minIndex] = min(averagedS21dBSmoothed(groupIndex, :));

    resonantFrequencyGHz(groupIndex) = selectedFrequencyGHz(minIndex);

    stdResonantFrequencyGHz(groupIndex) = std(singleTraceResFreqGHz(rowIndices));

    stdMinimumS21dB(groupIndex) = std(singleTraceMinS21dB(rowIndices));
end

%% VALIDATION DEVIATION AND TOTAL UNCERTAINTY

validationStartIndex = find(groups == "validation start", 1);
validationEndIndex = find(groups == "validation end", 1);

validationFrequencyDeviationGHz = abs( ...
    resonantFrequencyGHz(validationEndIndex) - ...
    resonantFrequencyGHz(validationStartIndex) ...
) / sqrt(3);

validationMagnitudeDeviationdB = abs( ...
    minimumS21dB(validationEndIndex) - ...
    minimumS21dB(validationStartIndex) ...
) / sqrt(3);

totalFrequencyUncertaintyGHz = sqrt( ...
    stdResonantFrequencyGHz.^2 + ...
    baselineFrequencyUncertaintyGHz.^2 + ...
    validationFrequencyDeviationGHz.^2 ...
);

totalMagnitudeUncertaintydB = sqrt( ...
    stdMinimumS21dB.^2 + ...
    baselineMagnitudeUncertaintydB.^2 + ...
    validationMagnitudeDeviationdB.^2 ...
);

%% EXPORT SUMMARY TABLE

analysisGroupIndices = setdiff( ...
    1:numGroups, ...
    [validationStartIndex, validationEndIndex] ...
);

summaryTable = table( ...
    groupLabels(analysisGroupIndices)', ...
    resonantFrequencyGHz(analysisGroupIndices)', ...
    stdResonantFrequencyGHz(analysisGroupIndices)', ...
    validationFrequencyDeviationGHz * ones(numel(analysisGroupIndices), 1), ...
    baselineFrequencyUncertaintyGHz * ones(numel(analysisGroupIndices), 1), ...
    totalFrequencyUncertaintyGHz(analysisGroupIndices)', ...
    minimumS21dB(analysisGroupIndices)', ...
    stdMinimumS21dB(analysisGroupIndices)', ...
    validationMagnitudeDeviationdB * ones(numel(analysisGroupIndices), 1), ...
    baselineMagnitudeUncertaintydB * ones(numel(analysisGroupIndices), 1), ...
    totalMagnitudeUncertaintydB(analysisGroupIndices)', ...
    'VariableNames', { ...
        'Group', ...
        'Res Freq GHz', ...
        'Repeat Freq Std', ...
        'Validation Freq Dev', ...
        'Baseline Freq Unc', ...
        'Total Freq Unc', ...
        'Min S21 dB', ...
        'Repeat Mag Std', ...
        'Validation Mag Dev', ...
        'Baseline Mag Unc', ...
        'Total Mag Unc' ...
    } ...
);

disp(summaryTable);
writetable(summaryTable, fullfile(resultsFolder, "resonator_s2p_summary.csv"));

%% PLOT AVERAGED S21 TRACES

figure;
hold on;

for groupIndex = plotGroupIndices
    plot( ...
        selectedFrequencyGHz, ...
        averagedS21dBSmoothed(groupIndex, :), ...
        "LineWidth", 2 ...
    );
end

hold off;
xlim(frequencyRangeGHz);
xlabel("Frequency (GHz)");
ylabel("S_{21} (dB)");
title(plotTitle);
legend(groupLabels(plotGroupIndices), "Location", "best");
fontsize(gca, 16, "points");
fontname(gca, "Times");
grid on;

saveas(gcf, fullfile(resultsFolder, "averaged_s21_traces.png"));
saveas(gcf, fullfile(resultsFolder, "averaged_s21_traces.fig"));

%% PLOT RESONANT FREQUENCY VS MINIMUM S21 MAGNITUDE

figure;
hold on;

for groupIndex = plotGroupIndices
    scatter( ...
        resonantFrequencyGHz(groupIndex), ...
        minimumS21dB(groupIndex), ...
        60, ...
        "LineWidth", 1.2 ...
    );
end

errorbar( ...
    resonantFrequencyGHz(plotGroupIndices), ...
    minimumS21dB(plotGroupIndices), ...
    totalMagnitudeUncertaintydB(plotGroupIndices), ...
    "vertical", ...
    ".", ...
    "Color", "k" ...
);

errorbar( ...
    resonantFrequencyGHz(plotGroupIndices), ...
    minimumS21dB(plotGroupIndices), ...
    totalFrequencyUncertaintyGHz(plotGroupIndices), ...
    "horizontal", ...
    ".", ...
    "Color", "k" ...
);

hold off;
xlabel("Resonant Frequency (GHz)");
ylabel("Minimum S_{21} (dB)");
title("Resonance Shift by Condition");
legend(groupLabels(plotGroupIndices), "Location", "best");
fontsize(gca, 16, "points");
fontname(gca, "Times");
grid on;

saveas(gcf, fullfile(resultsFolder, "resonance_shift_scatter.png"));
saveas(gcf, fullfile(resultsFolder, "resonance_shift_scatter.fig"));
