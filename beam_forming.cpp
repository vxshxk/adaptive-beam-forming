#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "matplotlib-cpp/matplotlibcpp.h"

#define PI 3.14159265358979323846

namespace plt = matplotlibcpp;

const int NUM_PHASE_CENTERS = 4; // Number of antenna elements
const int NUM_POINTS = 1000;     // Number of time samples
const int NUM_ANGLES = 360;      // Number of angles for beam pattern

using namespace std;

// Convert normalized gain to dB, ensuring we don't take log(0)
vector<double> convertToDb(const vector<double>& normalized_gain) {
    vector<double> gain_db(normalized_gain.size());
    for (size_t i = 0; i < normalized_gain.size(); i++) {
        gain_db[i] = 10 * log10(max(normalized_gain[i], std::numeric_limits<double>::epsilon())); 
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

    // Generate time vector
    for (int i = 0; i < NUM_POINTS; i++) {
        time[i] = i * 0.001; // Assume time step of 1ms
    }

    // Compute Hamming window weights
    for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
        hamming_weights[j] = 0.54 - 0.46 * cos((2 * PI * j) / (NUM_PHASE_CENTERS - 1));
    }

    // Generate individual input signals (sinusoids with different phase shifts)
    for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
        double phase_shift = j * PI / 6.0; // Different phase shifts
        for (int i = 0; i < NUM_POINTS; i++) {
            signals_real[j][i] = sin(2 * PI * 5 * time[i] + phase_shift);
        }
    }

    // Compute beamformed signals (with and without Hamming weighting)
    for (int i = 0; i < NUM_POINTS; i++) {
        for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
            beamformed_uniform[i] += signals_real[j][i]; // Uniform weighting
            beamformed_hamming[i] += hamming_weights[j] * signals_real[j][i]; // Hamming weighted
        }
        beamformed_uniform[i] /= NUM_PHASE_CENTERS; // Normalize
        beamformed_hamming[i] /= NUM_PHASE_CENTERS; // Normalize
    }

    // Compute Antenna Beam Pattern & Array Response vs. AoA
    for (int k = 0; k < NUM_ANGLES; k++) {
        angles[k] = k - 180; // Angle from -180 to 180 degrees
        double angle_rad = angles[k] * PI / 180.0;
        
        for (int j = 0; j < NUM_PHASE_CENTERS; j++) {
            double phase_term = j * PI * sin(angle_rad);
            beam_pattern_uniform[k] += cos(phase_term); // Without Hamming
            beam_pattern_hamming[k] += hamming_weights[j] * cos(phase_term); // With Hamming
        }

        // Normalize beam patterns
        beam_pattern_uniform[k] = pow(beam_pattern_uniform[k], 2);
        beam_pattern_hamming[k] = pow(beam_pattern_hamming[k], 2);

        // Compute normalized array response
        array_response_uniform[k] = beam_pattern_uniform[k] / beam_pattern_uniform[NUM_ANGLES / 2]; // Normalize to max value
        array_response_hamming[k] = beam_pattern_hamming[k] / beam_pattern_hamming[NUM_ANGLES / 2]; // Normalize to max value
    }

    // Plot individual input signals
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

    // Plot beamformed output signals with and without Hamming weighting
    plt::figure_size(1200, 600);
    plt::plot(time, beamformed_uniform, {{"label", "Beamformed Signal (Uniform)"}, {"color", "blue"}, {"linewidth", "2"}});
    plt::plot(time, beamformed_hamming, {{"label", "Beamformed Signal (Hamming)"}, {"color", "red"}, {"linewidth", "2"}});
    plt::title("Beamformed Output Signal (With and Without Hamming Weighting)");
    plt::xlabel("Time (s)");
    plt::ylabel("Amplitude");
    plt::legend();
    plt::grid(true);
    plt::show();

    vector<double> angles;
    for (double theta = -90; theta <= 90; theta += 1.0) {
        angles.push_back(theta);
    }

    // Simulated normalized array responses (dummy example, replace with actual computations)
    vector<double> array_response_uniform = {
        -25, -23, -21, -18, -15, -12, -9, -7, -5, -3, -2, -1, 0, -1, -2, -3, -5, -7, -9, -12, -15, -18, -21, -23, -25, 
        -26, -28, -30, -32, -33, -34, -35, -36, -37, -38, -39, -40, -41, -41, -42, -43, -43, -44, -45, -45, -46, -47, -47
    };
    vector<double> array_response_hamming = {
        -40, -38, -35, -30, -25, -20, -15, -12, -9, -6, -4, -2, 0, -2, -4, -6, -9, -12, -15, -20, -25, -30, -35, -38, -40, 
        -42, -43, -44, -45, -46, -47, -48, -48, -49, -50, -51, -51, -52, -53, -53, -54, -55, -55, -56, -57, -57, -58, -59, -59
    };
    // Convert responses to dB
    vector<double> array_response_uniform_db = convertToDb(array_response_uniform);
    vector<double> array_response_hamming_db = convertToDb(array_response_hamming);

    // Plot Normalized Array Response vs. Angle of Arrival
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
