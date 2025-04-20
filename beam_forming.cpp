#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "matplotlib-cpp/matplotlibcpp.h"

#define PI 3.14159265358979323846

namespace plt = matplotlibcpp;

const int NUM_PHASE_CENTERS = 4;
const int NUM_POINTS = 1000;
const int NUM_ANGLES = 360;

using namespace std;

vector<double> convertToDb(const vector<double>& normalized_gain) {
    vector<double> gain_db(normalized_gain.size());
    for (size_t i = 0; i < normalized_gain.size(); i++) {
        gain_db[i] = 10 * log10(max(normalized_gain[i], numeric_limits<double>::epsilon()));
    }
    return gain_db;
}

int main() {
    vector<double> time(NUM_POINTS);
    vector<vector<double>> signals_real(NUM_PHASE_CENTERS, vector<double>(NUM_POINTS));
    vector<double> beamformed_uniform(NUM_POINTS, 0.0);
    vector<double> beamformed_hamming(NUM_POINTS, 0.0);
    vector<double> hamming_weights(NUM_PHASE_CENTERS);

    vector<double> angles(NUM_ANGLES);
    vector<double> beam_pattern_uniform(NUM_ANGLES, 0.0);
    vector<double> beam_pattern_hamming(NUM_ANGLES, 0.0);
    vector<double> array_response_uniform(NUM_ANGLES, 0.0);
    vector<double> array_response_hamming(NUM_ANGLES, 0.0);

    // Time vector
    for (int i = 0; i < NUM_POINTS; i++) {
        time[i] = i * 0.001;
    }

    // Hamming weights
    for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
        hamming_weights[j] = 0.54 - 0.46 * cos((2 * PI * j) / (NUM_PHASE_CENTERS - 1));
    }

    // Input signals
    for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
        double phase_shift = j * PI / 6.0;
        for (int i = 0; i < NUM_POINTS; i++) {
            signals_real[j][i] = sin(2 * PI * 5 * time[i] + phase_shift);
        }
    }

    // Beamformed signals
    for (int i = 0; i < NUM_POINTS; i++) {
        for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
            beamformed_uniform[i] += signals_real[j][i];
            beamformed_hamming[i] += hamming_weights[j] * signals_real[j][i];
        }
        beamformed_uniform[i] /= NUM_PHASE_CENTERS;
        beamformed_hamming[i] /= NUM_PHASE_CENTERS;
    }

    // Beam pattern and array response
    for (int k = 0; k < NUM_ANGLES; k++) {
        angles[k] = k - 180;
        double angle_rad = angles[k] * PI / 180.0;

        for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
            double phase_term = j * PI * sin(angle_rad);
            beam_pattern_uniform[k] += cos(phase_term);
            beam_pattern_hamming[k] += hamming_weights[j] * cos(phase_term);
        }

        beam_pattern_uniform[k] = pow(beam_pattern_uniform[k], 2);
        beam_pattern_hamming[k] = pow(beam_pattern_hamming[k], 2);

        array_response_uniform[k] = beam_pattern_uniform[k];
        array_response_hamming[k] = beam_pattern_hamming[k];
    }

    // Normalize array responses
    double max_uniform = *max_element(array_response_uniform.begin(), array_response_uniform.end());
    double max_hamming = *max_element(array_response_hamming.begin(), array_response_hamming.end());

    for (int k = 0; k < NUM_ANGLES; k++) {
        array_response_uniform[k] /= max_uniform;
        array_response_hamming[k] /= max_hamming;
    }

    // Convert to dB
    vector<double> array_response_uniform_db = convertToDb(array_response_uniform);
    vector<double> array_response_hamming_db = convertToDb(array_response_hamming);

    // Plot individual signals
    plt::figure_size(1200, 600);
    for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
        plt::plot(time, signals_real[j], {{"label", "Signal " + to_string(j + 1)}});
    }
    plt::title("Individual Input Signals (Real Part)");
    plt::xlabel("Time (s)");
    plt::ylabel("Amplitude");
    plt::legend();
    plt::grid(true);
    plt::show();

    // Plot beamformed signals
    plt::figure_size(1200, 600);
    plt::plot(time, beamformed_uniform, {{"label", "Beamformed Signal (Uniform)"}, {"color", "blue"}, {"linewidth", "2"}});
    plt::plot(time, beamformed_hamming, {{"label", "Beamformed Signal (Hamming)"}, {"color", "red"}, {"linewidth", "2"}});
    plt::title("Beamformed Output Signal (With and Without Hamming Weighting)");
    plt::xlabel("Time (s)");
    plt::ylabel("Amplitude");
    plt::legend();
    plt::grid(true);
    plt::show();

    // Plot array response
    plt::figure_size(1200, 600);
    plt::plot(angles, array_response_uniform_db, {{"label", "Without Hamming"}, {"linestyle", "--"}, {"color", "blue"}, {"linewidth", "2"}});
    plt::plot(angles, array_response_hamming_db, {{"label", "With Hamming"}, {"color", "orange"}, {"linewidth", "2"}});
    plt::title("Antenna Pattern with and without Hamming Weighting");
    plt::xlabel("Angle of Arrival (degrees)");
    plt::ylabel("Normalized Array Response (dB)");
    plt::legend();
    plt::grid(true);
    plt::show();

    return 0;
}
