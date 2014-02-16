class Molby::Molecule

#
# UFF force field parameters
# Taken from RDKit source code
# http://rdkit.svn.sourceforge.net/viewvc/rdkit/trunk/Code/ForceField/UFF/Params.cpp?revision=2044&view=markup
# Accessed on 2012.5.8.
#
# Atom    an   r1     theta0  x1     D1     zeta    Z1     Vi     Uj    Xi      Hard    Radius Description
UFFParams = [
["H_",    1,   0.354, 180,    2.886, 0.044, 12,     0.712, 0,     0,    4.528,  6.9452, 0.371, "H generic"],
["H_b",   1,   0.46,  83.5,   2.886, 0.044, 12,     0.712, 0,     0,    4.528,  6.9452, 0.371, "H bridging B-H-B"],
["He4+4", 2,   0.849, 90,     2.362, 0.056, 15.24,  0.098, 0,     0,    9.66,   14.92,  1.3,   "He"],
["Li",    3,   1.336, 180,    2.451, 0.025, 12,     1.026, 0,     2,    3.006,  2.386,  1.557, "Li"],
["Be3+2", 4,   1.074, 109.47, 2.745, 0.085, 12,     1.565, 0,     2,    4.877,  4.443,  1.24,  "Be2+ tetrahedral"],
["B_3",   5,   0.838, 109.47, 4.083, 0.18,  12.052, 1.755, 0,     2,    5.11,   4.75,   0.822, "B tetrahedral"],
["B_2",   5,   0.828, 120,    4.083, 0.18,  12.052, 1.755, 0,     2,    5.11,   4.75,   0.822, "B trigonal"],
["C_3",   6,   0.757, 109.47, 3.851, 0.105, 12.73,  1.912, 2.119, 2,    5.343,  5.063,  0.759, "C sp3"],
["C_R",   6,   0.729, 120,    3.851, 0.105, 12.73,  1.912, 0,     2,    5.343,  5.063,  0.759, "C aromatic"],
["C_2",   6,   0.732, 120,    3.851, 0.105, 12.73,  1.912, 0,     2,    5.343,  5.063,  0.759, "C sp2"],
["C_1",   6,   0.706, 180,    3.851, 0.105, 12.73,  1.912, 0,     2,    5.343,  5.063,  0.759, "C sp"],
["N_3",   7,   0.7,   106.7,  3.66,  0.069, 13.407, 2.544, 0.45,  2,    6.899,  5.88,   0.715, "N sp3"],
["N_R",   7,   0.699, 120,    3.66,  0.069, 13.407, 2.544, 0,     2,    6.899,  5.88,   0.715, "N aromatic"],
["N_2",   7,   0.685, 111.2,  3.66,  0.069, 13.407, 2.544, 0,     2,    6.899,  5.88,   0.715, "N sp2"],
["N_1",   7,   0.656, 180,    3.66,  0.069, 13.407, 2.544, 0,     2,    6.899,  5.88,   0.715, "N sp"],
["O_3",   8,   0.658, 104.51, 3.5,   0.06,  14.085, 2.3,   0.018, 2,    8.741,  6.682,  0.669, "O tetrahedral"],
["O_3_z", 8,   0.528, 146,    3.5,   0.06,  14.085, 2.3,   0.018, 2,    8.741,  6.682,  0.669, "O tetrahedral (zeolite)"],
["O_R",   8,   0.68,  110,    3.5,   0.06,  14.085, 2.3,   0,     2,    8.741,  6.682,  0.669, "O aromatic"],
["O_2",   8,   0.634, 120,    3.5,   0.06,  14.085, 2.3,   0,     2,    8.741,  6.682,  0.669, "O trigonal"],
["O_1",   8,   0.639, 180,    3.5,   0.06,  14.085, 2.3,   0,     2,    8.741,  6.682,  0.669, "O linear"],
["F_",    9,   0.668, 180,    3.364, 0.05,  14.762, 1.735, 0,     2,    10.874, 7.474,  0.706, "F"],
["Ne4+4", 10,  0.92,  90,     3.243, 0.042, 15.44,  0.194, 0,     2,    11.04,  10.55,  1.768, "Ne"],
["Na",    11,  1.539, 180,    2.983, 0.03,  12,     1.081, 0,     1.25, 2.843,  2.296,  2.085, "Na"],
["Mg3+2", 12,  1.421, 109.47, 3.021, 0.111, 12,     1.787, 0,     1.25, 3.951,  3.693,  1.5,   "Mg2+ tetrahedral"],
["Al3",   13,  1.244, 109.47, 4.499, 0.505, 11.278, 1.792, 0,     1.25, 4.06,   3.59,   1.201, "Al tetrahedral"],
["Si3",   14,  1.117, 109.47, 4.295, 0.402, 12.175, 2.323, 1.225, 1.25, 4.168,  3.487,  1.176, "Si tetrahedral"],
["P_3+3", 15,  1.101, 93.8,   4.147, 0.305, 13.072, 2.863, 2.4,   1.25, 5.463,  4,      1.102, "P3+ tetrahedral"],
["P_3+5", 15,  1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.4,   1.25, 5.463,  4,      1.102, "P5+ tetrahedral"],
["P_3+q", 15,  1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.4,   1.25, 5.463,  4,      1.102, "P organometallic"],
["S_3+2", 16,  1.064, 92.1,   4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928,  4.486,  1.047, "S2+ tetrahedral"],
["S_3+4", 16,  1.049, 103.2,  4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928,  4.486,  1.047, "S4+ tetrahedral"],
["S_3+6", 16,  1.027, 109.47, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928,  4.486,  1.047, "S6+ tetrahedral"],
["S_R",   16,  1.077, 92.2,   4.035, 0.274, 13.969, 2.703, 0,     1.25, 6.928,  4.486,  1.047, "S aromatic"],
["S_2",   16,  0.854, 120,    4.035, 0.274, 13.969, 2.703, 0,     1.25, 6.928,  4.486,  1.047, "S trigonal"],
["Cl",    17,  1.044, 180,    3.947, 0.227, 14.866, 2.348, 0,     1.25, 8.564,  4.946,  0.994, "Cl"],
["Ar4+4", 18,  1.032, 90,     3.868, 0.185, 15.763, 0.3,   0,     1.25, 9.465,  6.355,  2.108, "Ar"],
["K_",    19,  1.953, 180,    3.812, 0.035, 12,     1.165, 0,     0.7,  2.421,  1.92,   2.586, "K"],
["Ca6+2", 20,  1.761, 90,     3.399, 0.238, 12,     2.141, 0,     0.7,  3.231,  2.88,   2,     "Ca2+ octahedral"],
["Sc3+3", 21,  1.513, 109.47, 3.295, 0.019, 12,     2.592, 0,     0.7,  3.395,  3.08,   1.75,  "Sc3+ tetrahedral"],
["Ti3+4", 22,  1.412, 109.47, 3.175, 0.017, 12,     2.659, 0,     0.7,  3.47,   3.38,   1.607, "Ti4+ tetrahedral"],
["Ti6+4", 22,  1.412, 90,     3.175, 0.017, 12,     2.659, 0,     0.7,  3.47,   3.38,   1.607, "Ti6+ octahedral"],
["V_3+5", 23,  1.402, 109.47, 3.144, 0.016, 12,     2.679, 0,     0.7,  3.65,   3.41,   1.47,  "V5+ tetrahedral"],
["Cr6+3", 24,  1.345, 90,     3.023, 0.015, 12,     2.463, 0,     0.7,  3.415,  3.865,  1.402, "Cr3+ octahedral"],
["Mn6+2", 25,  1.382, 90,     2.961, 0.013, 12,     2.43,  0,     0.7,  3.325,  4.105,  1.533, "Mn2+ octahedral"],
["Fe3+2", 26,  1.27,  109.47, 2.912, 0.013, 12,     2.43,  0,     0.7,  3.76,   4.14,   1.393, "Fe2+ tetrahedral"],
["Fe6+2", 26,  1.335, 90,     2.912, 0.013, 12,     2.43,  0,     0.7,  3.76,   4.14,   1.393, "Fe2+ octahedral"],
["Co6+3", 27,  1.241, 90,     2.872, 0.014, 12,     2.43,  0,     0.7,  4.105,  4.175,  1.406, "Co3+ octahedral"],
["Ni4+2", 28,  1.164, 90,     2.834, 0.015, 12,     2.43,  0,     0.7,  4.465,  4.205,  1.398, "Ni2+ square planar"],
["Cu3+1", 29,  1.302, 109.47, 3.495, 0.005, 12,     1.756, 0,     0.7,  4.2,    4.22,   1.434, "Cu+ tetrahedral"],
["Zn3+2", 30,  1.193, 109.47, 2.763, 0.124, 12,     1.308, 0,     0.7,  5.106,  4.285,  1.4,   "Zn2+ tetrahedral"],
["Ga3+3", 31,  1.26,  109.47, 4.383, 0.415, 11,     1.821, 0,     0.7,  3.641,  3.16,   1.211, "Ga3+ tetrahedral"],
["Ge3",   32,  1.197, 109.47, 4.28,  0.379, 12,     2.789, 0.701, 0.7,  4.051,  3.438,  1.189, "Ge tetrahedral"],
["As3+3", 33,  1.211, 92.1,   4.23,  0.309, 13,     2.864, 1.5,   0.7,  5.188,  3.809,  1.204, "As3+ tetrahedral"],
["Se3+2", 34,  1.19,  90.6,   4.205, 0.291, 14,     2.764, 0.335, 0.7,  6.428,  4.131,  1.224, "Se2+ tetrahedral"],
["Br",    35,  1.192, 180,    4.189, 0.251, 15,     2.519, 0,     0.7,  7.79,   4.425,  1.141, "Br"],
["Kr4+4", 36,  1.147, 90,     4.141, 0.22,  16,     0.452, 0,     0.7,  8.505,  5.715,  2.27,  "Kr"],
["Rb",    37,  2.26,  180,    4.114, 0.04,  12,     1.592, 0,     0.2,  2.331,  1.846,  2.77,  "Rb"],
["Sr6+2", 38,  2.052, 90,     3.641, 0.235, 12,     2.449, 0,     0.2,  3.024,  2.44,   2.415, "Sr2+ octahedral"],
["Y_3+3", 39,  1.698, 109.47, 3.345, 0.072, 12,     3.257, 0,     0.2,  3.83,   2.81,   1.998, "Y3+ tetrahedral"],
["Zr3+4", 40,  1.564, 109.47, 3.124, 0.069, 12,     3.667, 0,     0.2,  3.4,    3.55,   1.758, "Zr4+ tetrahedral"],
["Nb3+5", 41,  1.473, 109.47, 3.165, 0.059, 12,     3.618, 0,     0.2,  3.55,   3.38,   1.603, "Nb5+ tetrahedral"],
["Mo6+6", 42,  1.467, 90,     3.052, 0.056, 12,     3.4,   0,     0.2,  3.465,  3.755,  1.53,  "Mo6+ octahedral"],
["Mo3+6", 42,  1.484, 109.47, 3.052, 0.056, 12,     3.4,   0,     0.2,  3.465,  3.755,  1.53,  "Mo6+ tetrahedral"],
["Tc6+5", 43,  1.322, 90,     2.998, 0.048, 12,     3.4,   0,     0.2,  3.29,   3.99,   1.5,   "Tc5+ octahedral"],
["Ru6+2", 44,  1.478, 90,     2.963, 0.056, 12,     3.4,   0,     0.2,  3.575,  4.015,  1.5,   "Ru2+ octahedral"],
["Rh6+3", 45,  1.332, 90,     2.929, 0.053, 12,     3.5,   0,     0.2,  3.975,  4.005,  1.509, "Rh3+ octahedral"],
["Pd4+2", 46,  1.338, 90,     2.899, 0.048, 12,     3.21,  0,     0.2,  4.32,   4,      1.544, "Pd2+ square planar"],
["Ag1+1", 47,  1.386, 180,    3.148, 0.036, 12,     1.956, 0,     0.2,  4.436,  3.134,  1.622, "Ag+ linear"],
["Cd3+2", 48,  1.403, 109.47, 2.848, 0.228, 12,     1.65,  0,     0.2,  5.034,  3.957,  1.6,   "Cd2+ tetrahedral"],
["In3+3", 49,  1.459, 109.47, 4.463, 0.599, 11,     2.07,  0,     0.2,  3.506,  2.896,  1.404, "In3+ tetrahedral"],
["Sn3",   50,  1.398, 109.47, 4.392, 0.567, 12,     2.961, 0.199, 0.2,  3.987,  3.124,  1.354, "Sn tetrahedral"],
["Sb3+3", 51,  1.407, 91.6,   4.42,  0.449, 13,     2.704, 1.1,   0.2,  4.899,  3.342,  1.404, "Sb3+ tetrahedral"],
["Te3+2", 52,  1.386, 90.25,  4.47,  0.398, 14,     2.882, 0.3,   0.2,  5.816,  3.526,  1.38,  "Te2+ tetrahedral"],
["I_",    53,  1.382, 180,    4.5,   0.339, 15,     2.65,  0,     0.2,  6.822,  3.762,  1.333, "I"],
["Xe4+4", 54,  1.267, 90,     4.404, 0.332, 12,     0.556, 0,     0.2,  7.595,  4.975,  2.459, "Xe"],
["Cs",    55,  2.57,  180,    4.517, 0.045, 12,     1.573, 0,     0.1,  2.183,  1.711,  2.984, "Cs"],
["Ba6+2", 56,  2.277, 90,     3.703, 0.364, 12,     2.727, 0,     0.1,  2.814,  2.396,  2.442, "Ba2+ octahedral"],
["La3+3", 57,  1.943, 109.47, 3.522, 0.017, 12,     3.3,   0,     0.1,  2.8355, 2.7415, 2.071, "La3+ tetrahedral"],
["Ce6+3", 58,  1.841, 90,     3.556, 0.013, 12,     3.3,   0,     0.1,  2.774,  2.692,  1.925, "Ce3+ octahedral"],
["Pr6+3", 59,  1.823, 90,     3.606, 0.01,  12,     3.3,   0,     0.1,  2.858,  2.564,  2.007, "Pr3+ octahedral"],
["Nd6+3", 60,  1.816, 90,     3.575, 0.01,  12,     3.3,   0,     0.1,  2.8685, 2.6205, 2.007, "Nd3+ octahedral"],
["Pm6+3", 61,  1.801, 90,     3.547, 0.009, 12,     3.3,   0,     0.1,  2.881,  2.673,  2,     "Pm3+ octahedral"],
["Sm6+3", 62,  1.78,  90,     3.52,  0.008, 12,     3.3,   0,     0.1,  2.9115, 2.7195, 1.978, "Sm3+ octahedral"],
["Eu6+3", 63,  1.771, 90,     3.493, 0.008, 12,     3.3,   0,     0.1,  2.8785, 2.7875, 2.227, "Eu3+ octahedral"],
["Gd6+3", 64,  1.735, 90,     3.368, 0.009, 12,     3.3,   0,     0.1,  3.1665, 2.9745, 1.968, "Gd3+ octahedral"],
["Tb6+3", 65,  1.732, 90,     3.451, 0.007, 12,     3.3,   0,     0.1,  3.018,  2.834,  1.954, "Tb3+ octahedral"],
["Dy6+3", 66,  1.71,  90,     3.428, 0.007, 12,     3.3,   0,     0.1,  3.0555, 2.8715, 1.934, "Dy3+ octahedral"],
["Ho6+3", 67,  1.696, 90,     3.409, 0.007, 12,     3.416, 0,     0.1,  3.127,  2.891,  1.925, "Ho3+ octahedral"],
["Er6+3", 68,  1.673, 90,     3.391, 0.007, 12,     3.3,   0,     0.1,  3.1865, 2.9145, 1.915, "Er3+ octahedral"],
["Tm6+3", 69,  1.66,  90,     3.374, 0.006, 12,     3.3,   0,     0.1,  3.2514, 2.9329, 2,     "Tm3+ octahedral"],
["Yb6+3", 70,  1.637, 90,     3.355, 0.228, 12,     2.618, 0,     0.1,  3.2889, 2.965,  2.158, "Yb3+ octahedral"],
["Lu6+3", 71,  1.671, 90,     3.64,  0.041, 12,     3.271, 0,     0.1,  2.9629, 2.4629, 1.896, "Lu3+ octahedral"],
["Hf3+4", 72,  1.611, 109.47, 3.141, 0.072, 12,     3.921, 0,     0.1,  3.7,    3.4,    1.759, "Hf4+ tetrahedral"],
["Ta3+5", 73,  1.511, 109.47, 3.17,  0.081, 12,     4.075, 0,     0.1,  5.1,    2.85,   1.605, "Ta5+ tetrahedral"],
["W_6+6", 74,  1.392, 90,     3.069, 0.067, 12,     3.7,   0,     0.1,  4.63,   3.31,   1.538, "W6+ octahedral"],
["W_3+4", 74,  1.526, 109.47, 3.069, 0.067, 12,     3.7,   0,     0.1,  4.63,   3.31,   1.538, "W4+ tetrahedral"],
["W_3+6", 74,  1.38,  109.47, 3.069, 0.067, 12,     3.7,   0,     0.1,  4.63,   3.31,   1.538, "W6+ tetrahedral"],
["Re6+5", 75,  1.372, 90,     2.954, 0.066, 12,     3.7,   0,     0.1,  3.96,   3.92,   1.6,   "Re5+ octahedral"],
["Re3+7", 75,  1.314, 109.47, 2.954, 0.066, 12,     3.7,   0,     0.1,  3.96,   3.92,   1.6,   "Re7+ tetrahedral"],
["Os6+6", 76,  1.372, 90,     3.12,  0.037, 12,     3.7,   0,     0.1,  5.14,   3.63,   1.7,   "Os6+ octahedral"],
["Ir6+3", 77,  1.371, 90,     2.84,  0.073, 12,     3.731, 0,     0.1,  5,      4,      1.866, "Ir3+ octahedral"],
["Pt4+2", 78,  1.364, 90,     2.754, 0.08,  12,     3.382, 0,     0.1,  4.79,   4.43,   1.557, "Pt2+ square planar"],
["Au4+3", 79,  1.262, 90,     3.293, 0.039, 12,     2.625, 0,     0.1,  4.894,  2.586,  1.618, "Au3+ tetrahedral"],
["Hg1+2", 80,  1.34,  180,    2.705, 0.385, 12,     1.75,  0,     0.1,  6.27,   4.16,   1.6,   "Hg2+ linear"],
["Tl3+3", 81,  1.518, 120,    4.347, 0.68,  11,     2.068, 0,     0.1,  3.2,    2.9,    1.53,  "Tl3+ tetrahedral"],
["Pb3",   82,  1.459, 109.47, 4.297, 0.663, 12,     2.846, 0.1,   0.1,  3.9,    3.53,   1.444, "Pb tetrahedral"],
["Bi3+3", 83,  1.512, 90,     4.37,  0.518, 13,     2.47,  1,     0.1,  4.69,   3.74,   1.514, "Bi3+ tetrahedral"],
["Po3+2", 84,  1.5,   90,     4.709, 0.325, 14,     2.33,  0.3,   0.1,  4.21,   4.21,   1.48,  "Po2+ tetrahedral"],
["At",    85,  1.545, 180,    4.75,  0.284, 15,     2.24,  0,     0.1,  4.75,   4.75,   1.47,  "At"],
["Rn4+4", 86,  1.42,  90,     4.765, 0.248, 16,     0.583, 0,     0.1,  5.37,   5.37,   2.2,   "Rn"],
["Fr",    87,  2.88,  180,    4.9,   0.05,  12,     1.847, 0,     0,    2,      2,      2.3,   "Fr"],
["Ra6+2", 88,  2.512, 90,     3.677, 0.404, 12,     2.92,  0,     0,    2.843,  2.434,  2.2,   "Ra2+ octahedral"],
["Ac6+3", 89,  1.983, 90,     3.478, 0.033, 12,     3.9,   0,     0,    2.835,  2.835,  2.108, "Ac3+ octahedral"],
["Th6+4", 90,  1.721, 90,     3.396, 0.026, 12,     4.202, 0,     0,    3.175,  2.905,  2.018, "Th4+ octahedral"],
["Pa6+4", 91,  1.711, 90,     3.424, 0.022, 12,     3.9,   0,     0,    2.985,  2.905,  1.8,   "Pa4+ octahedral"],
["U_6+4", 92,  1.684, 90,     3.395, 0.022, 12,     3.9,   0,     0,    3.341,  2.853,  1.713, "U4+ octahedral"],
["Np6+4", 93,  1.666, 90,     3.424, 0.019, 12,     3.9,   0,     0,    3.549,  2.717,  1.8,   "Np4+ octahedral"],
["Pu6+4", 94,  1.657, 90,     3.424, 0.016, 12,     3.9,   0,     0,    3.243,  2.819,  1.84,  "Pu4+ octahedral"],
["Am6+4", 95,  1.66,  90,     3.381, 0.014, 12,     3.9,   0,     0,    2.9895, 3.0035, 1.942, "Am4+ octahedral"],
["Cm6+3", 96,  1.801, 90,     3.326, 0.013, 12,     3.9,   0,     0,    2.8315, 3.1895, 1.9,   "Cm3+ octahedral"],
["Bk6+3", 97,  1.761, 90,     3.339, 0.013, 12,     3.9,   0,     0,    3.1935, 3.0355, 1.9,   "Bk3+ octahedral"],
["Cf6+3", 98,  1.75,  90,     3.313, 0.013, 12,     3.9,   0,     0,    3.197,  3.101,  1.9,   "Cf3+ octahedral"],
["Es6+3", 99,  1.724, 90,     3.299, 0.012, 12,     3.9,   0,     0,    3.333,  3.089,  1.9,   "Es3+ octahedral"],
["Fm6+3", 100, 1.712, 90,     3.286, 0.012, 12,     3.9,   0,     0,    3.4,    3.1,    1.9,   "Fm3+ octahedral"],
["Md6+3", 101, 1.689, 90,     3.274, 0.011, 12,     3.9,   0,     0,    3.47,   3.11,   1.9,   "Md3+ octahedral"],
["No6+3", 102, 1.679, 90,     3.248, 0.011, 12,     3.9,   0,     0,    3.475,  3.175,  1.9,   "No3+ octahedral"],
["Lw6+3", 103, 1.698, 90,     3.236, 0.011, 12,     3.9,   0,     0,    3.5,    3.2,    1.9,   "Lr3+ octahedral"]
]

#  Calculate UFF bond length
def uff_bond_length(r1, r2, x1, x2, bond_order)
  bond_order_correction = -0.1332 * (r1 + r2) * log(bond_order)
  sq = sqrt(x1) - sqrt(x2)
  electronegativity_correction = r1 * r2 * (sqrt(x1) - sqrt(x2)) ** 2 / (x1 * r1 + x2 * r2)
  return r1 + r2 + bond_order_correction + electronegativity_correction
end

#  UFF bond force constant
def uff_bond_force(idx1, idx2, bond_order)
  r1 = UFFParams[idx1][2]
  r2 = UFFParams[idx2][2]
  r12 = uff_bond_length(r1, r2, UFFParams[idx1][10], UFFParams[idx2][10], bond_order)
  return 332.06 * UFFParams[idx1][7] * UFFParams[idx2][7] / (r12 ** 3)
end

#  UFF angle force constant (equilibrium angle should be given in degree)
def uff_angle_force(idx1, idx2, idx3, bond_order_12, bond_order_23, angle)
  r1 = UFFParams[idx1][2]
  r2 = UFFParams[idx2][2]
  r3 = UFFParams[idx3][2]
  cost = cos(angle * 3.1415927 / 180.0)
  r12 = uff_bond_length(r1, r2, UFFParams[idx1][10], UFFParams[idx2][10], bond_order_12)
  r23 = uff_bond_length(r2, r3, UFFParams[idx2][10], UFFParams[idx3][10], bond_order_23)
  r13 = sqrt(r12 * r12 + r23 * r23 - 2 * r12 * r23 * cost)
  return 332.06 * UFFParams[idx1][7] * UFFParams[idx3][7] / (r13 ** 5) * (r12 * r23 * (1.0 - cost * cost) - r13 * r13 * cost)
end

def guess_uff_parameter_dialog(current_value, indices)

  #  TODO: this dialog is soon to be made obsolete
  
  indices = indices.split(/-/)

  if indices.length == 2
    par_type = "bond"
  elsif indices.length == 3
    par_type = "angle"
  else
    message_box("UFF parameter guess of this type is not implemented.", "Guess UFF Force", "OK")
	return nil
  end
  
  #  Look up the UFFParams entry
  uff_types = indices.map { |i|
    ap = atoms[i]
    sel = []
    UFFParams.each_with_index { |p, j|
      if p[1] == ap.atomic_number
        sel.push(j)
      end
    }
    if sel.length == 0
      message_box("No UFF parameter is available for atom #{ap.name}", "Guess UFF Force", "OK")
      return nil
    end
    sel
  }
  names = indices.map { |i| sprintf("%d:%s", atoms[i].res_seq, atoms[i].name) }
  types = indices.map { |i| atoms[i].atom_type }
  utypes = indices.map { |i| atoms[i].uff_type }
  names_str = names.join("-")
  types_str = types.join("-")
  recalc = lambda { |i1, i2, i3, b1, b2|
    i1 = uff_types[0][i1] rescue 0
	i2 = uff_types[1][i2] rescue 0
	i3 = uff_types[2][i3] rescue 0
    if par_type == "bond"
      k = uff_bond_force(i1, i2, b1)
    else
      k = uff_angle_force(i1, i2, i3, b1, b2, current_value.to_f)
    end
    sprintf("%.5f", k)
  }
  hash = Dialog.run("Guess UFF Force", "Accept", "Cancel") {
    action_proc = lambda { |it|
	  t1 = value("uff_type_0").to_i rescue 0
	  t2 = value("uff_type_1").to_i rescue 0
	  t3 = value("uff_type_2").to_i rescue 0
	  b1 = value("bond_order_0").to_f rescue 1.0
	  b2 = value("bond_order_1").to_f rescue 1.0
	  current_value = value("current_value") rescue nil
      set_value("guessed", recalc.call(t1, t2, t3, b1, b2)) rescue nil
    }
    type_selects = []
    uff_types.each_with_index { |p, i|
      type_selects.push(item(:text, :title=>names[i] + "->"))
      if p.length == 1
        type_selects.push(item(:text, :title=>UFFParams[p[0]][13]))
      else
	    subitems = p.map { |pp| UFFParams[pp][0] }
		uff_idx = subitems.index(utypes[i]) || 0
        type_selects.push(item(:popup, :subitems=>p.map { |pp| UFFParams[pp][13] }, :tag=>"uff_type_#{i}", :action=>action_proc, :value=>uff_idx))
      end
    }
    bond_orders = []
    (0..indices.length - 2).each { |i|
      bond_orders.push(item(:text, :title=>names[i] + "-" + names[i + 1]))
      bond_orders.push(item(:textfield, :width=>60, :value=>"1.0", :tag=>"bond_order_#{i}", :action=>action_proc))
    }
    layout(1,
      item(:text, :title=>"Guess UFF force parameter from atom types and bond order"),
      layout(2,
        item(:text, :title=>"UFF Atom Types:"),
        item(:text, :title=>"Bond order:"),
        layout(2, *type_selects),
        layout(2, *bond_orders)
      ),
	  (par_type == "bond" ? nil :
		layout(2,
		  item(:text, :title=>"Equilibrium angle = "),
		  item(:textfield, :width=>100, :value=>current_value, :tag=>"current_value", :action=>action_proc))
	  ),
      layout(2,
        item(:text, :title=>"Guessed force = "),
        item(:textfield, :editable=>false, :width=>100, :value=>recalc.call(0, 0, 0, 1.0, 1.0), :tag=>"guessed")
      )
    )
  }
  if hash[:status] == 0
    3.times { |i|
	  idx = indices[i]
	  next unless idx
	  ii = uff_types[i][hash["uff_type_#{i}"].to_i]
	  if ii
        atoms[idx].uff_type = UFFParams[ii][0]
      end
	}
    return hash["guessed"], current_value
  else
    return nil
  end
end

def guess_uff_parameters
  mol = self
  #  Look up the non-conventional atoms
  exclude = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]
  arena = mol.md_arena
  arena.prepare(true)
  xatoms = IntGroup[]
  xbonds = xangles = xdihedrals = xfragments = []
  h = Dialog.new("Guess UFF Parameters: #{mol.name}", nil, nil, :resizable=>true) {
    update_xatoms = lambda {
      xatoms = mol.atom_group { |ap| !exclude.member?(ap.atomic_number) }
      xfragments = mol.fragments(xatoms)
      xbonds = (0...mol.nbonds).select { |i| b = mol.bonds[i]; xatoms.include?(b[0]) || xatoms.include?(b[1]) }
      xangles = (0...mol.nangles).select { |i| a = mol.angles[i]; xatoms.include?(a[0]) || xatoms.include?(a[1]) || xatoms.include?(a[2]) }
      xdihedrals = (0...mol.ndihedrals).select { |i| d = mol.dihedrals[i]; xatoms.include?(d[0]) || xatoms.include?(d[1]) || xatoms.include?(d[2]) || xatoms.include?(d[3]) }
      xbonds.each { |i| xatoms.add(mol.bonds[i]) }
      xangles.each { |i|
	    ang = mol.angles[i]
		xatoms.add(ang)
		2.times { |j|
		  b0 = ang[j]
		  b1 = ang[j + 1]
		  k = 0
		  mol.bonds.each { |b|
		    break if (b[0] == b0 && b[1] && b1) || (b[0] == b1 && b[1] && b0)
			k += 1
		  }
		  if k < mol.nbonds && !xbonds.include?(k)
		    xbonds.push(k)
		  end
		}
	  }
	  xbonds.sort!
      item_with_tag("table")[:refresh] = true
    }
    columns = {
      "atoms"=>[["no", 40], ["name", 60], ["type", 40], ["element", 40], ["int_charge", 40], ["uff_type", 120], ["weight", 70], ["eps", 70], ["r", 70], ["eps14", 70], ["r14", 70]],
      "bonds"=>[["no", 40], ["bond", 110], ["nums", 80], ["types", 90], ["par_types", 90], ["order", 80], ["k", 80], ["r0", 80], ["real_r", 80]],
      "angles"=>[["no", 40], ["angle", 140], ["nums", 100], ["types", 100], ["par_types", 100], ["k", 80], ["a0", 80], ["real_a", 80]],
      "dihedrals"=>[["no", 40], ["dihedral", 180], ["nums", 100], ["types", 110], ["par_types", 110], ["k", 50], ["period", 20], ["phi0", 50], ["real_phi", 60]],
      "fragments"=>[["no", 40], ["fragment", 360]]
    }
    tab = "atoms"
    update_selection = lambda {
      g = mol.selection
      sel = IntGroup[]
      case tab
      when "atoms"
        xatoms.each_with_index { |n, i|
          sel.add(i) if g.include?(n)
        }
      when "bonds"
        xbonds.each_with_index { |n, i|
          b = mol.bonds[n]
          sel.add(i) if g.include?(b[0]) && g.include?(b[1])
        }
      when "angles"
        xangles.each_with_index { |n, i|
          an = mol.angles[n]
          sel.add(i) if g.include?(an[0]) && g.include?(an[1]) && g.include?(an[2])
        }
      when "dihedrals"
        xdihedrals.each_with_index { |n, i|
          di = mol.dihedrals[n]
          sel.add(i) if g.include?(di[0]) && g.include?(di[1]) && g.include?(di[2]) && g.include?(di[3])
        }
      when "fragments"
        xfragments.each_with_index { |f, i|
          sel.add(i) if (f - g).count == 0
        }
      end
      it = item_with_tag("table")
      if sel != it[:selection]
        @dont_change_mol_selection = true
        it[:selection] = sel
      end
    }
    selection_changed = lambda { |it|
      if @dont_change_mol_selection
        @dont_change_mol_selection = false
        return
      end
      sel = it[:selection]
      g = IntGroup[]
      case tab
      when "atoms"
        sel.each { |idx|
          g.add(xatoms[idx])
        }
      when "bonds"
        sel.each { |idx|
          g.add(mol.bonds[xbonds[idx]])
        }
      when "angles"
        sel.each { |idx|
          g.add(mol.angles[xangles[idx]])
        }
      when "dihedrals"
        sel.each { |idx|
          g.add(mol.dihedrals[xdihedrals[idx]])
        }
      when "fragments"
        sel.each { |idx|
          g.add(xfragments[idx])
        }
      end
      mol.selection = g
    }
    select_tab = lambda { |tag|
      table = item_with_tag("table")
      table[:columns] = columns[tag]
      tab = tag
      table[:refresh] = true
      update_selection.call
    }
    tab_button_pressed = lambda { |it|
      ["atoms", "bonds", "angles", "dihedrals", "fragments"].each { |tag|
        next if tag == it[:tag]
        item_with_tag(tag)[:value] = 0
      }
      select_tab.call(it[:tag])
    }
    uff_popup_titles = []
    uff_type_for_title = Hash.new
    uff_title_for_type = Hash.new
    uff_title_for_type[""] = "-- select --"
    uff_popup = lambda { |an|
      if uff_popup_titles[an] == nil
	    if an == 0
		  titles = ["(no type)"]
		else
          titles = []
          Molby::Molecule::UFFParams.each { |u|
            if u[1] == an
              titles.push(u[13])
              uff_type_for_title[u[13]] = u[0]
              uff_title_for_type[u[0]] = u[13]
            end
          }
		end
        uff_popup_titles[an] = titles
      end
      uff_popup_titles[an]
    }
    get_count = lambda { |it|
      case tab
      when "atoms"
        return xatoms.count
      when "bonds"
        return xbonds.count
      when "angles"
        return xangles.count
      when "dihedrals"
        return xdihedrals.count
      when "fragments"
        return xfragments.count
      end
      return 0
    }
    get_value = lambda { |it, row, col|
      case tab
      when "atoms"
        idx = xatoms[row]
        ap0 = mol.atoms[idx]
        case col
        when 0
          return idx.to_s
        when 1
          return "#{ap0.res_seq}:#{ap0.name}"
        when 2
          return ap0.atom_type
        when 3
          return ap0.element
        when 4
          return ap0.int_charge.to_s
        when 5
          return (ap0.atomic_number == 0 ? "(no type)" : (uff_title_for_type[ap0.uff_type] || "..."))
        when 6
          return sprintf("%.3f", ap0.weight)
        when 7
          return sprintf("%.3f", arena.vdw_par(idx).eps)
        when 8
          return sprintf("%.3f", arena.vdw_par(idx).r_eq)
        when 9
          return sprintf("%.3f", arena.vdw_par(idx).eps14)
        when 10
          return sprintf("%.3f", arena.vdw_par(idx).r_eq14)
        end
      when "bonds"
        idx = xbonds[row]
        b = mol.bonds[idx]
        ap0 = mol.atoms[b[0]]
        ap1 = mol.atoms[b[1]]
        case col
        when 0
          return idx.to_s
        when 1
          return "#{ap0.res_seq}:#{ap0.name}-#{ap1.res_seq}:#{ap1.name}"
        when 2
          return "#{b[0]}-#{b[1]}"
        when 3
          return "#{ap0.atom_type}-#{ap1.atom_type}"
        when 4
          atom_types = arena.bond_par(idx).atom_types
          return "#{atom_types[0]}-#{atom_types[1]}"
        when 5
          return (o = mol.get_bond_order(idx)) && sprintf("%.3f", o)
        when 6
          return sprintf("%.3f", arena.bond_par(idx).k)
        when 7
          return sprintf("%.3f", arena.bond_par(idx).r0)
        when 8
          return sprintf("%.3f", mol.calc_bond(b[0], b[1]))
        end
      when "angles"
        idx = xangles[row]
        an = mol.angles[idx]
        ap0 = mol.atoms[an[0]]
        ap1 = mol.atoms[an[1]]
        ap2 = mol.atoms[an[2]]
        case col
        when 0
          return idx.to_s
        when 1
          return "#{ap0.res_seq}:#{ap0.name}-#{ap1.res_seq}:#{ap1.name}-#{ap2.res_seq}:#{ap2.name}"
        when 2
          return "#{an[0]}-#{an[1]}-#{an[2]}"
        when 3
          return "#{ap0.atom_type}-#{ap1.atom_type}-#{ap2.atom_type}"
        when 4
          atom_types = arena.angle_par(idx).atom_types
          return "#{atom_types[0]}-#{atom_types[1]}-#{atom_types[2]}"
        when 5
          return sprintf("%.3f", arena.angle_par(idx).k)
        when 6
          return sprintf("%.2f", arena.angle_par(idx).a0)
        when 7
          return sprintf("%.2f", mol.calc_angle(an[0], an[1], an[2]))
        end
      when "dihedrals"
        idx = xdihedrals[row]
        di = mol.dihedrals[idx]
        ap0 = mol.atoms[di[0]]
        ap1 = mol.atoms[di[1]]
        ap2 = mol.atoms[di[2]]
        ap3 = mol.atoms[di[3]]
        case col
        when 0
          return idx.to_s
        when 1
          return "#{ap0.res_seq}:#{ap0.name}-#{ap1.res_seq}:#{ap1.name}-#{ap2.res_seq}:#{ap2.name}-#{ap2.res_seq}:#{ap2.name}"
        when 2
          return "#{di[0]}-#{di[1]}-#{di[2]}-#{di[3]}"
        when 3
          return "#{ap0.atom_type}-#{ap1.atom_type}-#{ap2.atom_type}-#{ap3.atom_type}"
        when 4
          atom_types = arena.dihedral_par(idx).atom_types
          return "#{atom_types[0]}-#{atom_types[1]}-#{atom_types[2]}-#{atom_types[3]}"
        when 5
          return sprintf("%.3f", arena.dihedral_par(idx).k)
        when 6
          return arena.dihedral_par(idx).period.to_s
        when 7
          return sprintf("%.2f", arena.dihedral_par(idx).phi0)
        when 8
          return sprintf("%.2f", mol.calc_dihedral(di[0], di[1], di[2], di[3]))
        end
      when "fragments"
        case col
        when 0
          return row.to_s
        when 1
          return xfragments[row].to_s[9..-2]  #  Remove "IntGroup[" and "]"
        end
      end
      return "..."
    }
    is_item_editable = lambda { |it, row, col|
      case tab
      when "atoms"
        case col
        when 1..4, 6..10
          return true
        end
      when "bonds"
        case col
        when 5..7
          return true
        end
      when "angles"
        case col
        when 5..6
          return true
        end
      when "dihedrals"
        case col
        when 5..7
          return true
        end
      end
      return false
    }
    modify_parameter = lambda { |mol, partype, idx, attr, val|
      arena = mol.md_arena
      case partype
      when "vdw"
        pref = arena.vdw_par(idx)
        pen = mol.parameter.vdws
      when "bond"
        pref = arena.bond_par(idx)
        pen = mol.parameter.bonds
      when "angle"
        pref = arena.angle_par(idx)
        pen = mol.parameter.angles
      when "dihedral"
        pref = arena.dihedral_par(idx)
        pen = mol.parameter.dihedrals
      when "improper"
        pref = arena.improper_par(idx)
        pen = mol.parameter.impropers
      end
      pref_new = pen.lookup(pref.atom_types, :create, :local, :missing, :nobasetype, :nowildcard)
      pref.keys.each { |k|
        next if k == :source || k == :index || k == :par_type
        pref_new.set_attr(k, pref.get_attr(k))
      }
      if attr != nil
        pref_new.set_attr(attr, val)
        arena.prepare(true)
      end
      return pref_new
    }
    set_value = lambda { |it, row, col, val|
      case tab
      when "atoms"
        idx = xatoms[row]
        ap0 = mol.atoms[idx]
        case col
        when 1
          val = val.to_s
          if val =~ /\d:/
            val = Regexp.last_match.post_match
          end
          ap0.name = val
        when 2
          ap0.atom_type = val
        when 3
          ap0.element = val
        when 4
          ap0.int_charge = val.to_i
        when 6
          ap0.weight = val.to_f
        when 7
          pp = modify_parameter.call(mol, "vdw", idx, :eps, val)
          if pp.eps14 == 0.0
            pp.eps14 = val.to_f
          end
        when 8
          pp = modify_parameter.call(mol, "vdw", idx, :r_eq, val)
          if pp.r_eq14 == 0.0
            pp.r_eq14 = val.to_f
          end
        when 9
          modify_parameter.call(mol, "vdw", idx, :eps14, val)
        when 10
          modify_parameter.call(mol, "vdw", idx, :r_eq14, val)
        end
      when "bonds"
        idx = xbonds[row]
        case col
        when 5
          mol.assign_bond_order(idx, val.to_f)
        when 6
          modify_parameter.call(mol, "bond", idx, :k, val)
        when 7
          modify_parameter.call(mol, "bond", idx, :r0, val)
        end
      when "angles"
        idx = xangles[row]
        case col
        when 5
          modify_parameter.call(mol, "angle", idx, :k, val)
        when 6
          modify_parameter.call(mol, "angle", idx, :a0, val)
        end
      when "dihedrals"
        idx = xdihedrals[row]
        case col
        when 5
          modify_parameter.call(mol, "dihedral", idx, :k, val)
        when 6
          modify_parameter.call(mol, "dihedral", idx, :period, val)
        when 7
          modify_parameter.call(mol, "dihedral", idx, :phi0, val)
        end
      end
    }
    has_popup_menu = lambda { |it, row, col|
      if tab == "atoms" && col == 5
        #  UFF type popup
        ap = mol.atoms[xatoms[row]]
        val = uff_popup.call(ap.atomic_number)
        return val
      else
        return nil
      end
    }
    popup_menu_selected = lambda { |it, row, col, sel|
      return if tab != "atoms" || col != 5
      ap = mol.atoms[xatoms[row]]
      title = uff_popup.call(ap.atomic_number)[sel]
      ap.uff_type = (uff_type_for_title[title] || "")
    }
    guess_uff_types = lambda { |recalc_all|
      xatoms.each { |idx|
        ap = mol.atoms[idx]
        next if !recalc_all && ap.uff_type != ""
        u = uff_popup.call(ap.atomic_number)
        if u.length == 1
          ap.uff_type = (uff_type_for_title[u[0]] || "")
          next
        end
        case ap.atom_type
        when "c3", "cx", "cy"
          ap.uff_type = "C_3"
        when "ca", "cc", "cd", "cp", "cq"
          ap.uff_type = "C_R"
        when "c2", "ce", "cf", "cu", "cv", "c"
          ap.uff_type = "C_2"
        when "c1", "cg", "ch"
          ap.uff_type = "C_1"
        when "n3", "n4", "nh"
          ap.uff_type = "N_3"
        when "nb", "nc", "nd"
          ap.uff_type = "N_R"
        when "n", "n2", "na", "ne", "nf"
          ap.uff_type = "N_2"
        when "n1"
          ap.uff_type = "N_1"
        when "oh", "os", "ow"
          ap.uff_type = "O_2"
        else
          ap.uff_type = ""
        end
      }
    }
    set_color = lambda { |it, row, col|
      @red_color ||= [1.0, 0.2, 0.2]
      @yellow_color ||= [1.0, 1.0, 0.6]
      arena = mol.md_arena
      case tab
      when "atoms"
        pp = arena.vdw_par(xatoms[row])
      when "bonds"
        pp = arena.bond_par(xbonds[row])
      when "angles"
        pp = arena.angle_par(xangles[row])
      when "dihedrals"
        pp = arena.dihedral_par(xdihedrals[row])
      end
      src = pp.source
      if src == nil
        return [nil, @yellow_color]
      elsif src == false
        return [nil, @red_color]
      else
        return nil
      end
    }
    guess_parameters_for_fragments = lambda { |*d|
      name = mol.name
      xfragments.each_with_index { |frag, i|
        fmol = mol.extract(frag)
		mol.selection = frag
		frag_str = frag.to_s[9..-2]  #  Remove "IntGroup[" and "]"
		mes = "Guess MM/MD Parameters for #{mol.name}.fragment #{i}"
		n = fmol.ambertools_dialog("antechamber", mes, frag_str)
		break if n == 0
		next if n == -1
        n = fmol.invoke_antechamber(false, mes)
        break if n == 1
        calc_charge = get_global_settings("antechamber.calc_charge").to_i
        guess_atom_types = get_global_settings("antechamber.guess_atom_types").to_i
        if calc_charge
          #  Copy partial charges
          frag.each_with_index { |n, i|
            mol.atoms[n].charge = fmol.atoms[i].charge
          }
        end
        if guess_atom_types
          #  Copy atom types and local parameters
          frag.each_with_index { |n, i|
            mol.atoms[n].atom_type = fmol.atoms[i].atom_type
          }
          [:bond, :angle, :dihedral, :improper, :vdw].each { |ptype|
            case ptype
            when :bond
              pen = fmol.parameter.bonds
            when :angle
              pen = fmol.parameter.angles
            when :dihedral
              pen = fmol.parameter.dihedrals
            when :improper
              pen = fmol.parameter.impropers
            when :vdw
              pen = fmol.parameter.vdws
            end
            pen.each { |pref|
              next if pref.source != nil
              pref_new = mol.parameter.lookup(ptype, pref.atom_types, :local, :missing, :create, :nowildcard, :nobasetype)
              if pref_new != nil
                pref.keys.each { |k|
                  next if k == :index || k == :par_type || k == :source
                  pref_new.set_attr(k, pref.get_attr(k))
                }
              end
            }
          }
		  guess_uff_types.call(false)
        end
      }
    }
    guess_parameters_for_metals = lambda { |*d|
      catch(:exit) {
        #  Atoms
        xatoms.each { |idx|
		  pref = arena.vdw_par(idx) 
		  next if pref.source != false   #  Already defined
          ap0 = mol.atoms[idx]
          next if exclude.member?(ap0.atomic_number)
		  next if ap0.anchor_list != nil
          uff_type = ap0.uff_type
          u = UFFParams.find { |u| u[0] == uff_type }
          if u == nil
            error_message_box("The UFF type for atom #{idx} (#{ap0.name}) is not defined.")
            throw(:exit)
          end
          pref = mol.parameter.lookup(:vdw, ap0.index, :local, :missing, :create, :nowildcard, :nobasetype)
          pref.atom_type = idx
          pref.eps = pref.eps14 = u[5]
          pref.r_eq = pref.r_eq14 = u[4] * 0.5
		  pref.atomic_number = ap0.atomic_number
		  pref.weight = ap0.weight
        }
		#  Bonds
		pars = []
		xbonds.each { |idx|
		  pref = arena.bond_par(idx)
		  next if pref.source != false && pref.k > 0.0   #  Already defined
		  b = mol.bonds[idx]
		  is = []
		  aps = [mol.atoms[b[0]], mol.atoms[b[1]]]
		  2.times { |i|
		    if aps[i].anchor_list != nil
			  is[i] = -1
			  next
			end
		    uff_type = aps[i].uff_type
		    UFFParams.each_with_index { |u, j|
			  if u[0] == uff_type
			    is[i] = j
				break
			  end
			}
		    if is[i] == nil
			  error_message_box("The UFF type for atom #{b[i]} (#{aps[i].name}) is not defined.")
			  throw(:exit)
			end
		  }
		  bo = mol.get_bond_order(idx)
		  if bo == nil || bo == 0.0
		    bo = 1.0
		  end
		  if pref.source != false && pref.r0 > 0.0
		    len = pref.r0
	      else
		    len = mol.calc_bond(b[0], b[1])
	      end
		  if is[0] == -1 && is[1] == -1
		    #  Bond between anchors: no force
			force = 0.0
	      elsif is[0] == -1 || is[1] == -1
		    i = (is[0] == -1 ? 1 : 0)
			case aps[i].atomic_number
			when 0..23
			  force = 135.0
			when 24
			  force = 150.0
			when 25..27
			  force = 205.0
			when 28..36
			  force = 140.0
			when 37..54
			  force = 205.0
			else
			  force = 260.0
			end
		  else
		    force = mol.uff_bond_force(is[0], is[1], bo)
		  end
		  pars[idx] = [b, force, len]
		}
		pars.each { |pp|
		  next if pp == nil
		  pref = mol.parameter.lookup(:bond, pp[0], :local, :missing, :create, :nowildcard, :nobasetype)
		  if pref.source == false
 		    pref.atom_types = pp[0]
		  end
		  pref.k = pp[1]
		  pref.r0 = pp[2]
	    }
		#  Angles
		pars.clear
		xangles.each { |idx|
		  pref = arena.angle_par(idx)
		  next if pref.source != false && pref.k > 0.0   #  Already defined
		  a = mol.angles[idx]
		  is = []
		  aps = [mol.atoms[a[0]], mol.atoms[a[1]], mol.atoms[a[2]]]
		  3.times { |i|
		    if aps[i].anchor_list != nil
			  is[i] = -1
			  next
			end
		    uff_type = aps[i].uff_type
		    UFFParams.each_with_index { |u, j|
			  if u[0] == uff_type
			    is[i] = j
				break
			  end
			}
		    if is[i] == nil
			  error_message_box("The UFF type for atom #{a[i]} (#{aps[i].name}) is not defined.")
			  throw(:exit)
			end
		  }
		  bos = [1.0, 1.0]
		  2.times { |i|
		    #  Look up the bond and get the bond order
			mol.bonds.each_with_index { |b, j|
			  if (a[i] == b[0] && a[i + 1] == b[1]) || (a[i] == b[1] && a[i + 1] == b[0])
			    bos[i] = mol.get_bond_order(j)
				if bos[i] == nil || bos[i] == 0.0
				  bos[i] = 1.0
				end
			  end
			}
		  }
		  if pref.source != false && pref.a0 > 0.0
		    ang = pref.a0
		  else
		    ang = mol.calc_angle(a[0], a[1], a[2])
		  end
		  met = [aps[0].atomic_number, aps[1].atomic_number, aps[2].atomic_number].max
		  if is[1] == -1 && is[0] != -1 && is[2] != 0
		    #  Metal-Dummy-C angle
			case met
			when 0..23
			  force = 85.0
			when 24
			  force = 93.0
			when 25..27
			  force = 100.0
			when 28..36
			  force = 28.0
			when 37..54
			  force = 125.0
			else
			  force = 130.0
			end
		  elsif is[1] != -1 && is[0] == -1 && is[2] == -1
		    #  Dummy-Metal-Dummy angle
			case met
			when 0..23
			  force = 50.0
			when 24
			  force = 40.0
			when 25..27
			  force = 40.0
			when 28..36
			  force = 40.0
			when 37..54
			  force = 52.0
			else
			  force = 56.0
			end
		  elsif is[0] >= 0 && is[1] >= 0 && is[2] >= 0
		    force = mol.uff_angle_force(is[0], is[1], is[2], bos[0], bos[1], ang)
		  else
		    force = 0.0
		  end
		  pars[idx] = [a, force, ang]
		}
		pars.each { |pp|
		  next if pp == nil
		  pref = mol.parameter.lookup(:angle, pp[0], :local, :missing, :create, :nowildcard, :nobasetype)
		  if pref.source == false
		    pref.atom_types = pp[0]
		  end
		  pref.k = pp[1]
		  pref.a0 = pp[2]
	    }
		xdihedrals.each { |idx|
		  pref = arena.dihedral_par(idx)
		  next if pref.source != false   #  Already defined
		  d = mol.dihedrals[idx]
		  is = []
		  aps = [mol.atoms[d[0]], mol.atoms[d[1]], mol.atoms[d[2]], mol.atoms[d[3]]]
		  4.times { |i|
		    if aps[i].anchor_list != nil
			  is[i] = -1
			  next
			end
		    uff_type = aps[i].uff_type
		    UFFParams.each_with_index { |u, j|
			  if u[0] == uff_type
			    is[i] = j
				break
			  end
			}
		  }
		  if is[1] == -1 && is[2] == -1 && is[0] != -1 && is[3] != -1
		    #  X-##-##-X dihedral
			met = (aps[1].connects & aps[2].connects)[0]
			if met != nil
			  case mol.atoms[met].atomic_number
			  when 0..36
			    force = 0.36
		      when 37..54
			    force = 3.4
			  else
			    force = 3.4
			  end
			  pref = mol.parameter.lookup(:dihedral, ["X", d[1], d[2], "X"], :local, :missing, :create)
			  pref.atom_types = ["X", d[1], d[2], "X"]
			  pref.period = 5
			  pref.phi0 = 180.0
			  pref.k = force
			  arena.prepare(true)
			end
		  else
		    pref = mol.parameter.lookup(:dihedral, d, :local, :missing, :create, :nowildcard, :nobasetype)
			pref.atom_types = d
			pref.period = 2
			pref.phi0 = 180.0
			pref.k = 0.0
		  end
	    }
      }
	  arena.prepare(true)
	  item_with_tag("table")[:refresh] = true
    }
    layout(1,
#      layout(2,
#        item(:text, :title=>"Total charge: "),
#        item(:textfield, :width=>"80", :tag=>"total_charge")),
      layout(5,
        item(:togglebutton, :width=>80, :height=>24, :title=>"Atoms", :tag=>"atoms",
          :value=>1,
          :action=>tab_button_pressed),
        item(:togglebutton, :width=>80, :height=>24, :title=>"Bonds", :tag=>"bonds",
          :action=>tab_button_pressed),
        item(:togglebutton, :width=>80, :height=>24, :title=>"Angles", :tag=>"angles", 
          :action=>tab_button_pressed),
        item(:togglebutton, :width=>80, :height=>24, :title=>"Dihedrals", :tag=>"dihedrals", 
          :action=>tab_button_pressed),
        item(:togglebutton, :width=>80, :height=>24, :title=>"Fragments", :tag=>"fragments", 
          :action=>tab_button_pressed),
        :padding=>0, :margin=>0),
      item(:table, :width=>640, :height=>240, :flex=>[0,0,0,0,1,1], :tag=>"table", 
        :columns=>columns["atoms"],
        :on_count=>get_count,
        :on_get_value=>get_value,
        :on_set_value=>set_value,
        :is_item_editable=>is_item_editable,
        :on_selection_changed=>selection_changed,
        :has_popup_menu=>has_popup_menu,
        :on_popup_menu_selected=>popup_menu_selected,
        :on_set_color=>set_color),
      item(:button, :title=>"Run Antechamber for Non-Metal Fragments",
        :action=>guess_parameters_for_fragments, :flex=>[1,1,1,0,0,0], :align=>:center),
      item(:button, :title=>"Guess UFF Parameters for Metal Atoms",
        :action=>guess_parameters_for_metals, :flex=>[1,1,1,0,0,0], :align=>:center),
      item(:button, :title=>"Close", :action=>lambda { |it| hide }, :flex=>[1,1,1,0,0,0], :align=>:center),
      :flex=>[0,0,0,0,1,1]
    )
    size = self.size
    set_min_size(size)
    set_size(size[0] + 100, size[1] + 50);
    listen(mol, "documentModified", lambda { |d| update_xatoms.call; update_selection.call })
    listen(mol, "documentWillClose", lambda { |d| hide })
    update_xatoms.call
    guess_uff_types.call(true)
    show
  }
end

end

