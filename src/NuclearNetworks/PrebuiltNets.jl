networks::Dict{Symbol,NuclearNetwork} = Dict()

# Basic nets based on reactions from the Kippenhahn textbook
networks[:basic_H_burn] = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
networks[:basic_He_burn] = NuclearNetwork([:He4, :C12, :O16, :Ne20], [(:kipp_rates, :kipp_3alphaA99),
                                                         (:kipp_rates, :kipp_C12alpha),
                                                         (:kipp_rates, :kipp_O16alpha)])
networks[:basic_CO_burn] = NuclearNetwork([:He4, :C12, :O16, :Mg24], [(:kipp_rates, :kipp_CC),
                                                                      (:kipp_rates, :kipp_OO)])

# Various nets based on reactions from JINA
networks[:JINA_PPI] = NuclearNetwork([:H1, :D2, :He3, :He4], [(:jina_rates, :H1_H1_to_D2_betplus_w_x_0),
                                                              (:jina_rates, :H1_D2_to_He3_de04_n_x_0),
                                                              (:jina_rates, :H1_D2_to_He3_de04_x_x_0),
                                                              (:jina_rates, :He3_He3_to_H1_H1_He4_nacr_n_x_0),
                                                            ])
networks[:JINA_PPII] = NuclearNetwork([:H1, :He3, :He4, :Li7, :Be7], [(:jina_rates, :He4_He3_to_Be7_cd08_n_x_0),
                                                              (:jina_rates, :He4_He3_to_Be7_cd08_n_x_1),
                                                              (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_0),
                                                              (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_0),
                                                              (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_1),
                                                              (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_1),
                                                            ])
networks[:JINA_PPIII] = NuclearNetwork([:H1, :He4, :Be7, :Be8, :B8], [(:jina_rates, :H1_Be7_to_B8_nacr_r_x_0),
                                                              (:jina_rates, :H1_Be7_to_B8_nacr_n_x_0),
                                                              (:jina_rates, :B8_to_Be8_wc17_w_x_0),
                                                              (:jina_rates, :B8_to_He4_He4_wc12_w_x_0),
                                                            ])
networks[:JINA_PPIV] = NuclearNetwork([:H1, :He3, :He4], [(:jina_rates, :H1_He3_to_He4_betplus_w_x_0)])

networks[:JINA_CNOI] = NuclearNetwork([:H1, :He4, :C12, :C13, :N13, :N14, :N15, :O15],
                                            [
                                                (:jina_rates, :H1_C12_to_N13_ls09_r_x_0),
                                                (:jina_rates, :H1_C12_to_N13_ls09_n_x_0),
                                                (:jina_rates, :N13_to_C13_wc12_w_x_0),
                                                (:jina_rates, :H1_C13_to_N14_nacr_r_x_0),
                                                (:jina_rates, :H1_C13_to_N14_nacr_r_x_1),
                                                (:jina_rates, :H1_C13_to_N14_nacr_n_x_0),
                                                (:jina_rates, :H1_N14_to_O15_im05_r_x_0),
                                                (:jina_rates, :H1_N14_to_O15_im05_n_x_0),
                                                (:jina_rates, :H1_N14_to_O15_im05_n_x_1),
                                                (:jina_rates, :H1_N14_to_O15_im05_r_x_1),
                                                (:jina_rates, :O15_to_N15_wc12_w_x_0),
                                                (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_0),
                                                (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_1),
                                                (:jina_rates, :H1_N15_to_He4_C12_nacr_r_x_2),
                                                (:jina_rates, :H1_N15_to_He4_C12_nacr_n_x_0),
                                            ])

networks[:JINA_CNOII] = NuclearNetwork([:H1, :He4, :N14, :N15, :O15, :O16, :O17, :F17],
                                            [
                                                (:jina_rates, :H1_N14_to_O15_im05_r_x_0),
                                                (:jina_rates, :H1_N14_to_O15_im05_n_x_0),
                                                (:jina_rates, :H1_N14_to_O15_im05_n_x_1),
                                                (:jina_rates, :H1_N14_to_O15_im05_r_x_1),
                                                (:jina_rates, :O15_to_N15_wc12_w_x_0),
                                                (:jina_rates, :H1_N15_to_O16_li10_r_x_0),
                                                (:jina_rates, :H1_N15_to_O16_li10_r_x_1),
                                                (:jina_rates, :H1_N15_to_O16_li10_n_x_0),
                                                (:jina_rates, :H1_O16_to_F17_ia08_n_x_0),
                                                (:jina_rates, :F17_to_O17_wc12_w_x_0),
                                                (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_0),
                                                (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_1),
                                                (:jina_rates, :H1_O17_to_He4_N14_il10_r_x_2),
                                                (:jina_rates, :H1_O17_to_He4_N14_il10_n_x_0),
                                            ])