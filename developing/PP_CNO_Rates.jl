net = NuclearNetwork([:H1, :D2, :He3, :He4,
:Be7, :Li7, :B8
# :C12,   :C13,   
# :N13,   :N14,   :N15,
# :O14,   :O15,   :O16,   :O17,   :O18,   
# :F17,   :F18,   :F19,   
],
[
# PP I
(:jina_rates, :H1_H1_to_D2_betplus_w_x_0),
(:jina_rates, :H1_H1_to_D2_xxec_w_x_0),
(:jina_rates, :H1_D2_to_He3_de04_n_x_0),
# (:jina_rates, :H1_D2_to_He3_de04_x_x_0),
(:jina_rates, :He3_He3_to_H1_H1_He4_nacr_n_x_0),
# PP II
(:jina_rates, :He4_He3_to_Be7_cd08_n_x_0),
(:jina_rates, :He4_He3_to_Be7_cd08_n_x_1),
(:jina_rates, :Be7_to_Li7_xxec_w_x_0),
(:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_0),
(:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_0),
# (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_1),
# (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_1),
# PP III
(:jina_rates, :H1_Be7_to_B8_nacr_r_x_0),
(:jina_rates, :H1_Be7_to_B8_nacr_n_x_0),
(:jina_rates, :B8_to_He4_He4_wc12_w_x_0),
# PP IV
(:jina_rates, :H1_He3_to_He4_betplus_w_x_0),


# CNO Cycle 1

# (:jina_rates, :H1_C12_to_N13_ls09_r_x_0),
# (:jina_rates, :H1_C12_to_N13_ls09_n_x_0),

# (:jina_rates, :N13_to_C13_wc12_w_x_0),

# (:jina_rates, :H1_C13_to_N14_nacr_r_x_0),
# (:jina_rates, :H1_C13_to_N14_nacr_r_x_1),
# (:jina_rates, :H1_C13_to_N14_nacr_n_x_0),

# (:jina_rates, :H1_N14_to_O15_im05_r_x_0),
# (:jina_rates, :H1_N14_to_O15_im05_n_x_0),
# (:jina_rates, :H1_N14_to_O15_im05_n_x_1),
# (:jina_rates, :H1_N14_to_O15_im05_r_x_1),

# (:jina_rates, :O15_to_N15_wc12_w_x_0),

# (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_0),
# (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_1),
# (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_2),
# (:jina_rates, :H1_N15_to_He4_C12_nacr_n_x_0),

# CNO Cycle 2

# (:jina_rates, :H1_N14_to_O15_im05_r_x_0),
# (:jina_rates, :H1_N14_to_O15_im05_n_x_0),
# (:jina_rates, :H1_N14_to_O15_im05_n_x_1),
# (:jina_rates, :H1_N14_to_O15_im05_r_x_1),

# (:jina_rates, :O15_to_N15_wc12_w_x_0),

# (:jina_rates, :H1_N15_to_O16_li10_r_x_0),
# (:jina_rates, :H1_N15_to_O16_li10_r_x_1),
# (:jina_rates, :H1_N15_to_O16_li10_n_x_0),

# (:jina_rates, :H1_O16_to_F17_ia08_n_x_0),

# (:jina_rates, :F17_to_O17_wc12_w_x_0),

# (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_0),
# (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_1),
# (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_2),
# (:jina_rates, :H1_O17_to_He4_N14_il10_n_x_0),

# CNO Cycle 3

# (:jina_rates, :H1_N15_to_O16_li10_r_x_0),
# (:jina_rates, :H1_N15_to_O16_li10_r_x_1),
# (:jina_rates, :H1_N15_to_O16_li10_n_x_0),

# (:jina_rates, :H1_O16_to_F17_ia08_n_x_0),

# (:jina_rates, :F17_to_O17_wc12_w_x_0),

# (:jina_rates, :H1_O17_to_F18_il10_r_x_0),
# (:jina_rates, :H1_O17_to_F18_il10_r_x_1),
# (:jina_rates, :H1_O17_to_F18_il10_n_x_0),

# (:jina_rates, :F18_to_O18_wc12_w_x_0),

# (:jina_rates, :H1_O18_to_He4_N15_il10_n_x_0),
# (:jina_rates, :H1_O18_to_He4_N15_il10_r_x_0),
# (:jina_rates, :H1_O18_to_He4_N15_il10_r_x_1),
# (:jina_rates, :H1_O18_to_He4_N15_il10_r_x_2),

# CNO Cycle 4

# (:jina_rates, :H1_O16_to_F17_ia08_n_x_0),

# (:jina_rates, :F17_to_O17_wc12_w_x_0),

# (:jina_rates, :H1_O17_to_F18_il10_r_x_0),
# (:jina_rates, :H1_O17_to_F18_il10_r_x_1),
# (:jina_rates, :H1_O17_to_F18_il10_n_x_0),

# (:jina_rates, :F18_to_O18_wc12_w_x_0),

# (:jina_rates, :H1_O18_to_F19_il10_r_x_0),
# (:jina_rates, :H1_O18_to_F19_il10_r_x_1),
# (:jina_rates, :H1_O18_to_F19_il10_r_x_2),
# (:jina_rates, :H1_O18_to_F19_il10_n_x_0),

# (:jina_rates, :H1_F19_to_He4_O16_nacr_r_x_0),
# (:jina_rates, :H1_F19_to_He4_O16_nacr_x_x_0),
# (:jina_rates, :H1_F19_to_He4_O16_nacr_x_x_1),
# (:jina_rates, :H1_F19_to_He4_O16_nacr_x_x_2),
# (:jina_rates, :H1_F19_to_He4_O16_nacr_x_x_3),])
])