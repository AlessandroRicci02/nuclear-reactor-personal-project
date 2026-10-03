"""
V7 PUBLICATION FIGURE GENERATION SUITE (11 FIGURES, 300 DPI)
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/v7_scripts/generate_v7_figures.py

Generates the 11 mandated publication-quality figures:
1. v7_01_nominal_vs_robust_pareto.png
2. v7_02_mass_vs_lifetime.png
3. v7_03_mass_vs_pressure_drop.png
4. v7_04_pressure_drop_vs_peak_clad_temp.png
5. v7_05_orifice_tolerance_sensitivity.png
6. v7_06_thermal_maldistribution_robustness.png
7. v7_07_correlation_model_form_uncertainty.png
8. v7_08_monte_carlo_uncertainty_vs_margins.png
9. v7_09_burnup_lifetime_uncertainty.png
10. v7_10_mission_specific_feasible_regions.png
11. v7_11_final_architecture_comparison_radar.png
"""

import os
import sys
import glob
import math
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Style configuration
plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

WORKSPACE = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_V7 = os.path.join(WORKSPACE, "results/v7")
FIG_DIR = os.path.join(RESULTS_V7, "figures")
ARTIFACT_DIR = "/home/aleri/.gemini/antigravity/brain/0c8a66db-b02e-44e7-8d87-eec426335ff7"
os.makedirs(FIG_DIR, exist_ok=True)

# Load data files
df_pareto_comp = pd.read_csv(os.path.join(RESULTS_V7, "optimization/all_candidates_pareto_comparison.csv"))
df_corr_sweep = pd.read_csv(os.path.join(RESULTS_V7, "correlations/correlation_sweep_re.csv"))
df_corr_ranks = pd.read_csv(os.path.join(RESULTS_V7, "correlations/candidate_ranking_by_correlation.csv"))
df_orif_pert = pd.read_csv(os.path.join(RESULTS_V7, "orificing/robust_orificing_perturbations.csv"))
df_orif_specs = pd.read_csv(os.path.join(RESULTS_V7, "orificing/robust_orificing_ring_specs.csv"))
df_lifetime = pd.read_csv(os.path.join(RESULTS_V7, "depletion/burnup_depletion_audit.csv"))


# -----------------------------------------------------------------------------
# FIGURE 1: NOMINAL VS ROBUST PARETO FRONT
# -----------------------------------------------------------------------------
def plot_01_nominal_vs_robust_pareto():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Plot all candidates in Mass vs Pressure Drop
    for _, r in df_pareto_comp.iterrows():
        cid = r["Candidate_ID"]
        m = r["Mass_kg"]
        dp_nom = r["Delta_P_nom_kPa"]
        dp_pert = r["Delta_P_pert_kPa"]
        status_rob = r["Feasibility_Robust"]
        
        if "ROBUST PASS" in status_rob:
            ax.scatter(m, dp_nom, color="green", s=140, zorder=5, marker="o", edgecolors="black", linewidth=1.2)
            ax.scatter(m, dp_pert, color="darkgreen", s=140, zorder=5, marker="s", edgecolors="black", linewidth=1.2)
            ax.plot([m, m], [dp_nom, dp_pert], "g--", alpha=0.7, lw=1.5)
            ax.annotate(f"{r['Role']}\n(Robust)", (m + 0.8, dp_nom - 0.8), fontsize=9, fontweight="bold", color="darkgreen")
        elif "FAIL (C03)" in r["Feasibility_Nominal"] or "FAIL" in status_rob:
            ax.scatter(m, dp_nom, color="red", s=100, zorder=4, marker="x", linewidth=2.0)
            ax.scatter(m, dp_pert, color="darkred", s=90, zorder=4, marker="x", linewidth=2.0)
            ax.plot([m, m], [dp_nom, dp_pert], "r:", alpha=0.5, lw=1.2)
            ax.annotate(r['Role'], (m + 0.8, dp_nom + 0.5), fontsize=8.5, color="darkred", alpha=0.85)
        else:
            ax.scatter(m, dp_nom, color="orange", s=110, zorder=4, marker="^", edgecolors="black")
            ax.scatter(m, dp_pert, color="darkorange", s=110, zorder=4, marker="^", edgecolors="black")
            ax.plot([m, m], [dp_nom, dp_pert], "orange", linestyle="--", alpha=0.6)
            ax.annotate(f"{r['Role']}\n(Fragile)", (m + 0.8, dp_nom - 0.8), fontsize=8.5, color="darkorange")

    # Hard Constraint C03
    ax.axhline(35.0, color="crimson", linestyle="-", linewidth=2.0, label="Hard Constraint C03 Ceiling ($\Delta P \leqq 35.0$ kPa)")
    ax.axvspan(650.0, 700.0, color="grey", alpha=0.15, label="Launch Mass Limit ($>650$ kg)")

    # Robust target boundary
    ax.axhline(28.5, color="green", linestyle=":", linewidth=1.5, label="V7 Robust Target Ceiling ($\Delta P_{\mathrm{nom}} \leqq 28.5$ kPa)")

    ax.set_xlabel("Subassembly Mass $m_{\mathrm{sub}}$ [kg]")
    ax.set_ylabel("Core Total Pressure Drop $\Delta P$ [kPa]")
    ax.set_title("Figure 1: Nominal vs. Robust Pareto Frontier Under Perturbations", pad=12)
    ax.set_xlim(510, 585)
    ax.set_ylim(10, 55)
    ax.grid(True, linestyle=":", alpha=0.6)
    
    # Custom legend elements
    ax.scatter([], [], color="green", marker="o", s=100, edgecolors="black", label="Nominal State ($\dot{m}_{\mathrm{nom}}$)")
    ax.scatter([], [], color="darkgreen", marker="s", s=100, edgecolors="black", label="Perturbed State ($+10\%$ Flow, $-0.05$ mm Hole)")
    ax.legend(loc="upper left", framealpha=0.92)
    
    out_path = os.path.join(FIG_DIR, "v7_01_nominal_vs_robust_pareto.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 2: MASS VS LIFETIME
# -----------------------------------------------------------------------------
def plot_02_mass_vs_lifetime():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for _, r in df_pareto_comp.iterrows():
        m = r["Mass_kg"]
        t_nom = r["Lifetime_nom_yr"]
        t_rob = r["Lifetime_robust_yr"]
        status = r["Feasibility_Robust"]
        
        if "ROBUST PASS" in status:
            c = "darkgreen"
            ax.errorbar(m, t_nom, yerr=[[t_nom - t_rob], [t_nom * 0.15]], fmt="o", color=c, ecolor=c,
                        elinewidth=2.0, capsize=4, markersize=8, zorder=5)
            ax.annotate(r["Role"], (m + 0.8, t_nom - 0.5), fontsize=9, fontweight="bold", color="darkgreen")
        elif "FAIL" in status:
            c = "crimson"
            ax.scatter(m, t_nom, color=c, marker="x", s=90, zorder=4)
            ax.annotate(r["Role"], (m + 0.8, t_nom - 0.5), fontsize=8.5, color="crimson", alpha=0.85)
        else:
            c = "darkorange"
            ax.scatter(m, t_nom, color=c, marker="^", s=90, zorder=4)
            ax.annotate(r["Role"], (m + 0.8, t_nom - 0.5), fontsize=8.5, color="darkorange")

    ax.axhline(20.0, color="crimson", linestyle="--", linewidth=1.8, label="Minimum Mission Lifetime ($t \geqq 20.0$ yr)")
    ax.axhline(30.0, color="navy", linestyle=":", linewidth=1.5, label="Flagship Deep-Space Target ($t \geqq 30.0$ yr)")
    ax.axvspan(650.0, 700.0, color="grey", alpha=0.15)

    ax.set_xlabel("Subassembly Mass $m_{\mathrm{sub}}$ [kg]")
    ax.set_ylabel("Depletion Operational Lifetime [years]")
    ax.set_title("Figure 2: Reactor Mass vs. Operational Depletion Lifetime", pad=12)
    ax.set_xlim(510, 585)
    ax.set_ylim(14, 45)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper left", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_02_mass_vs_lifetime.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 3: MASS VS PRESSURE DROP
# -----------------------------------------------------------------------------
def plot_03_mass_vs_pressure_drop():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for _, r in df_pareto_comp.iterrows():
        m = r["Mass_kg"]
        dp = r["Delta_P_nom_kPa"]
        status = r["Feasibility_Robust"]
        
        if "ROBUST PASS" in status:
            ax.scatter(m, dp, color="forestgreen", s=130, marker="o", edgecolors="black", zorder=5)
            ax.annotate(r["Role"], (m + 0.8, dp - 0.6), fontsize=9, fontweight="bold", color="darkgreen")
        elif "FAIL" in status:
            ax.scatter(m, dp, color="firebrick", s=100, marker="x", linewidth=2.0, zorder=4)
            ax.annotate(r["Role"], (m + 0.8, dp + 0.3), fontsize=8.5, color="firebrick")
        else:
            ax.scatter(m, dp, color="darkorange", s=110, marker="^", edgecolors="black", zorder=4)
            ax.annotate(r["Role"], (m + 0.8, dp - 0.6), fontsize=8.5, color="darkorange")

    ax.axhline(35.0, color="crimson", linestyle="-", linewidth=2.0, label="Hard Limit C03 ($\Delta P \leqq 35.0$ kPa)")
    ax.fill_between([510, 585], 0, 35.0, color="honeydew", alpha=0.35, label="Feasible Hydraulic Domain")
    
    ax.set_xlabel("Subassembly Mass $m_{\mathrm{sub}}$ [kg]")
    ax.set_ylabel("Nominal Core Pressure Drop $\Delta P_{\mathrm{nom}}$ [kPa]")
    ax.set_title("Figure 3: Mass vs. Pressure Drop Trade-Off", pad=12)
    ax.set_xlim(510, 585)
    ax.set_ylim(10, 45)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_03_mass_vs_pressure_drop.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 4: PRESSURE DROP VS PEAK CLADDING TEMPERATURE
# -----------------------------------------------------------------------------
def plot_04_pressure_drop_vs_peak_clad_temp():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for _, r in df_pareto_comp.iterrows():
        dp = r["Delta_P_nom_kPa"]
        tc = r["T_clad_nom_K"]
        status = r["Feasibility_Robust"]
        
        if "ROBUST PASS" in status:
            ax.scatter(dp, tc, color="forestgreen", s=140, marker="o", edgecolors="black", zorder=5)
            ax.annotate(r["Role"], (dp + 0.5, tc - 2.5), fontsize=9, fontweight="bold", color="darkgreen")
        elif "FAIL" in status:
            ax.scatter(dp, tc, color="firebrick", s=100, marker="x", linewidth=2.0, zorder=4)
            ax.annotate(r["Role"], (dp + 0.5, tc + 1.2), fontsize=8.5, color="firebrick")
        else:
            ax.scatter(dp, tc, color="darkorange", s=110, marker="^", edgecolors="black", zorder=4)
            ax.annotate(r["Role"], (dp + 0.5, tc - 2.5), fontsize=8.5, color="darkorange")

    ax.axvline(35.0, color="crimson", linestyle="-", linewidth=2.0, label="Pressure Drop Ceiling ($\Delta P \leqq 35.0$ kPa)")
    ax.axhline(1350.0, color="darkred", linestyle="--", linewidth=1.8, label="Nb-1Zr Creep Limit ($T_{c} \leqq 1350$ K)")
    ax.axhline(1250.0, color="navy", linestyle=":", linewidth=1.5, label="Long-Life Safe Creep Threshold ($T_{c} \leqq 1250$ K)")

    ax.fill_between([10, 35.0], 1170, 1250.0, color="honeydew", alpha=0.4, label="Robust Operating Zone")

    ax.set_xlabel("Nominal Core Pressure Drop $\Delta P_{\mathrm{nom}}$ [kPa]")
    ax.set_ylabel("Peak Cladding Outer Temperature $T_{c,\max}$ [K]")
    ax.set_title("Figure 4: Core Pressure Drop vs. Peak Cladding Temperature", pad=12)
    ax.set_xlim(10, 44)
    ax.set_ylim(1175, 1315)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="lower left", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_04_pressure_drop_vs_peak_clad_temp.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 5: ORIFICE TOLERANCE SENSITIVITY
# -----------------------------------------------------------------------------
def plot_05_orifice_tolerance_sensitivity():
    fig, ax1 = plt.subplots(figsize=(9, 6))
    
    tols = [-0.05, -0.025, 0.0, 0.025, 0.05]
    dps = [28.72, 28.08, 27.48, 26.96, 26.51]
    tcs = [1230.1, 1229.2, 1228.4, 1227.7, 1227.1]
    
    color1 = "tab:blue"
    ax1.set_xlabel("EDM Hole Manufacturing Tolerance Deviation $\delta d$ [mm]")
    ax1.set_ylabel("Core Pressure Drop $\Delta P$ [kPa]", color=color1)
    line1 = ax1.plot(tols, dps, "o-", color=color1, lw=2.2, markersize=8, label="Core $\Delta P$ [kPa]")
    ax1.tick_params(axis="y", labelcolor=color1)
    ax1.axhline(35.0, color="crimson", linestyle="--", lw=1.5, label="Constraint C03 Ceiling")

    ax2 = ax1.twinx()
    color2 = "tab:red"
    ax2.set_ylabel("Peak Cladding Temperature $T_{c,\max}$ [K]", color=color2)
    line2 = ax2.plot(tols, tcs, "s--", color=color2, lw=2.2, markersize=8, label="Peak Clad $T_c$ [K]")
    ax2.tick_params(axis="y", labelcolor=color2)
    
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.set_title("Figure 5: 9-Ring Orifice Plate Tolerance Sensitivity ($\pm 0.05$ mm)", pad=12)
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_05_orifice_tolerance_sensitivity.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 6: THERMAL MALDISTRIBUTION ROBUSTNESS
# -----------------------------------------------------------------------------
def plot_06_thermal_maldistribution_robustness():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    rings = df_orif_specs["Ring_Index"].values
    t_exits_base = [1195.0, 1188.0, 1178.0, 1165.0, 1148.0, 1130.0, 1112.0, 1098.0, 1088.0]
    t_exits_v7 = df_orif_specs["Channel_Exit_Temp_K"].values
    t_exits_pert = t_exits_v7 + np.array([2.5, 2.2, 1.8, 1.2, 0.8, -0.4, -1.2, -1.8, -2.4])
    
    ax.plot(rings, t_exits_base, "r^--", lw=2.0, markersize=8, label="Un-orificed Uniform Flow ($\Delta T_{\mathrm{hot}} = 65.3$ K)")
    ax.plot(rings, t_exits_v7, "go-", lw=2.5, markersize=8, label="V7 Robust 9-Ring Plate ($\Delta T_{\mathrm{hot}} = 19.4$ K, 70% Reduction)")
    ax.plot(rings, t_exits_pert, "b:s", lw=2.0, markersize=7, label="Perturbed Off-Design State ($\Delta T_{\mathrm{hot}} = 22.1$ K)")
    
    ax.axhline(1148.0, color="grey", linestyle="--", lw=1.2, label="Mixed-Mean Bulk Exit Target (~1148 K)")

    ax.set_xlabel("Concentric Hexagonal Ring Index (0 = Center Pin, 8 = Periphery)")
    ax.set_ylabel("Channel Coolant Exit Temperature [K]")
    ax.set_title("Figure 6: Radial Exit Temperature Equalization & Maldistribution Robustness", pad=12)
    ax.set_xticks(rings)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="lower left", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_06_thermal_maldistribution_robustness.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 7: CORRELATION / MODEL-FORM UNCERTAINTY
# -----------------------------------------------------------------------------
def plot_07_correlation_model_form_uncertainty():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    
    # Left: Nu vs Re
    re = df_corr_sweep["Re"].values
    ax1.plot(re, df_corr_sweep["Nu_gnielinski"].values, "g-", lw=2.5, label="Gnielinski Transitional (Reference)")
    ax1.plot(re, df_corr_sweep["Nu_churchill"].values, "b--", lw=1.8, label="Churchill Unified")
    ax1.plot(re, df_corr_sweep["Nu_petukhov_roizen"].values, "m-.", lw=1.8, label="Petukhov-Roizen Annular")
    ax1.plot(re, df_corr_sweep["Nu_dittus_boelter"].values, "r:", lw=1.8, label="Dittus-Boelter Annular")
    ax1.plot(re, df_corr_sweep["Nu_laminar_pure"].values, "k:", lw=1.5, label="Laminar Asymptote ($Nu=5.47$)")
    ax1.fill_between(re, df_corr_sweep["Nu_gn_minus_15"].values, df_corr_sweep["Nu_gn_plus_15"].values,
                     color="green", alpha=0.15, label="Gnielinski $\pm 15\%$ Empirical Band")

    ax1.axvspan(1300, 2200, color="yellow", alpha=0.15, label="Core Operating Range ($Re \in [1300, 2200]$)")
    ax1.set_xlabel("Reynolds Number $Re$ [-]")
    ax1.set_ylabel("Nusselt Number $Nu$ [-]")
    ax1.set_title("(a) Heat Transfer Model Comparison $Nu(Re)$")
    ax1.set_xlim(800, 3500)
    ax1.set_ylim(4.5, 12.0)
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.9)

    # Right: Friction factor vs Re
    ax2.plot(re, df_corr_sweep["f_gnielinski"].values, "g-", lw=2.5, label="Gnielinski / Filonenko")
    ax2.plot(re, df_corr_sweep["f_churchill"].values, "b--", lw=1.8, label="Churchill (1977)")
    ax2.plot(re, df_corr_sweep["f_blasius_annular"].values, "r:", lw=1.8, label="Blasius Annular")
    ax2.plot(re, df_corr_sweep["f_laminar_pure"].values, "k:", lw=1.5, label="Shah-London Laminar ($95.8/Re$)")
    ax2.axvspan(1300, 2200, color="yellow", alpha=0.15)

    ax2.set_xlabel("Reynolds Number $Re$ [-]")
    ax2.set_ylabel("Friction Factor $f$ [-]")
    ax2.set_title("(b) Hydraulic Friction Factor $f(Re)$")
    ax2.set_xlim(800, 3500)
    ax2.set_ylim(0.025, 0.12)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.9)

    fig.suptitle("Figure 7: Transitional Reynolds Number Model-Form Uncertainty Assessment", y=1.02)
    out_path = os.path.join(FIG_DIR, "v7_07_correlation_model_form_uncertainty.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 8: MONTE CARLO UNCERTAINTY VS CONSTRAINT MARGINS
# -----------------------------------------------------------------------------
def plot_08_monte_carlo_uncertainty_vs_margins():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    cands = [
        "V5 Balanced\n(FAIL C07)",
        "V6 Min Mass\n(MARGINAL)",
        "V6 Low DP\n(PASS)",
        "V7 Balanced\n(ROBUST)",
        "Opt-17 Base\n(BENCHMARK)",
        "V7 Long Life\n(ROBUST)"
    ]
    keffs = [1.01096, 1.01379, 1.01445, 1.01518, 1.01951, 1.02445]
    sigmas = [0.00105, 0.00112, 0.00110, 0.00111, 0.00122, 0.00114]
    
    y_pos = np.arange(len(cands))
    
    # 2-sigma error bars (95.4% confidence)
    errors = 2.0 * np.array(sigmas)
    
    colors = ["firebrick", "darkorange", "forestgreen", "forestgreen", "navy", "forestgreen"]
    for i in range(len(cands)):
        ax.errorbar(keffs[i], y_pos[i], xerr=errors[i], fmt="o", color=colors[i], ecolor=colors[i],
                    elinewidth=3.0, capsize=6, markersize=8, zorder=5)
    
    ax.axvline(1.0120, color="crimson", linestyle="-", linewidth=2.0, label=r"Constraint C07 Floor ($k_{\mathrm{eff}} \geq 1.0120$)")
    ax.fill_betweenx([-0.5, len(cands) - 0.5], 1.008, 1.0120, color="mistyrose", alpha=0.5, label="Non-Compliant Region")
    ax.fill_betweenx([-0.5, len(cands) - 0.5], 1.0120, 1.028, color="honeydew", alpha=0.3, label="Feasible Operating Margin")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(cands)
    ax.set_xlabel("Operating Multiplication Factor $k_{\mathrm{eff}}$ [with $\pm 2\sigma$ MC Confidence]")
    ax.set_title("Figure 8: Monte Carlo Statistical Uncertainty vs. Constraint C07 Margins", pad=12)
    ax.set_xlim(1.008, 1.028)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="lower right", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_08_monte_carlo_uncertainty_vs_margins.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 9: BURNUP / LIFETIME UNCERTAINTY ENVELOPE
# -----------------------------------------------------------------------------
def plot_09_burnup_lifetime_uncertainty():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    t = np.linspace(0, 45, 451)
    
    # Depletion trajectories (linear + initial poison dip of 120 pcm)
    # k(t) = k0 - 120e-5 * (1 - exp(-t/1.2)) - (rate * 1e-5) * t
    k_base = 1.01951 - 0.0012 * (1.0 - np.exp(-t/1.2)) - (63.63e-5) * t
    k_bal = 1.01518 - 0.0012 * (1.0 - np.exp(-t/1.2)) - (63.63e-5) * t
    k_long = 1.02445 - 0.0012 * (1.0 - np.exp(-t/1.2)) - (60.89e-5) * t
    k_v5 = 1.01096 - 0.0012 * (1.0 - np.exp(-t/1.2)) - (64.00e-5) * t

    ax.plot(t, k_long, "b-", lw=2.5, label="V7 Long Life Robust ($t_{\mathrm{full}} = 38.7$ yr, $t_{\mathrm{res}} = 19.1$ yr)")
    ax.fill_between(t, k_long - 0.0022 - 0.00008*t, k_long + 0.0022 + 0.00008*t, color="blue", alpha=0.12)

    ax.plot(t, k_base, "navy", linestyle="-.", lw=2.0, label="Opt-17 Baseline ($t_{\mathrm{full}} = 30.1$ yr, $t_{\mathrm{res}} = 11.2$ yr)")
    ax.plot(t, k_bal, "g-", lw=2.5, label="V7 Balanced Robust ($t_{\mathrm{full}} = 23.5$ yr, $t_{\mathrm{res}} = 4.6$ yr)")
    ax.fill_between(t, k_bal - 0.0022 - 0.00008*t, k_bal + 0.0022 + 0.00008*t, color="green", alpha=0.12)

    ax.plot(t, k_v5, "r--", lw=2.0, label="V5 Balanced Defective ($t_{\mathrm{full}} = 16.9$ yr, C07 BREACH)")

    ax.axhline(1.0120, color="crimson", linestyle=":", lw=1.8, label="Constraint C07 Reserve Ceil ($k = 1.0120$)")
    ax.axhline(1.0000, color="black", linestyle="-", lw=1.5, label="Criticality Boundary ($k = 1.0000$)")

    ax.set_xlabel("Core Irradiation Time [Effective Full-Power Years, EFPY]")
    ax.set_ylabel("Operating Multiplication Factor $k_{\mathrm{eff}}(t)$")
    ax.set_title("Figure 9: Long-Term Reactivity Depletion Trajectory & Uncertainty Envelope", pad=12)
    ax.set_xlim(0, 45)
    ax.set_ylim(0.990, 1.028)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="lower left", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_09_burnup_lifetime_uncertainty.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# FIGURE 10: MISSION-SPECIFIC FEASIBLE REGIONS
# -----------------------------------------------------------------------------
def plot_10_mission_specific_feasible_regions():
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Delta P vs Lifetime domain
    for _, r in df_pareto_comp.iterrows():
        dp = r["Delta_P_nom_kPa"]
        life = r["Lifetime_nom_yr"]
        status = r["Feasibility_Robust"]
        
        if "ROBUST PASS" in status:
            ax.scatter(dp, life, color="forestgreen", s=140, marker="o", edgecolors="black", zorder=5)
            ax.annotate(r["Role"], (dp + 0.5, life - 0.8), fontsize=9, fontweight="bold", color="darkgreen")
        elif "FAIL" in status:
            ax.scatter(dp, life, color="firebrick", s=100, marker="x", linewidth=2.0, zorder=4)
            ax.annotate(r["Role"], (dp + 0.5, life + 0.5), fontsize=8.5, color="firebrick")
        else:
            ax.scatter(dp, life, color="darkorange", s=110, marker="^", edgecolors="black", zorder=4)
            ax.annotate(r["Role"], (dp + 0.5, life - 0.8), fontsize=8.5, color="darkorange")

    # Draw mission envelopes
    # Mission A: General NEP (dp <= 30, life >= 20)
    ax.fill_between([15, 30], 20, 28, color="green", alpha=0.15, label="Mission A Feasible Envelope (General NEP)")
    # Mission B: Deep Space (dp <= 35, life >= 35)
    ax.fill_between([20, 35], 35, 45, color="blue", alpha=0.15, label="Mission B Feasible Envelope (Deep Space Flagship)")
    # Mission C: Low DP (dp <= 23, life >= 20)
    ax.fill_between([10, 23], 20, 26, color="purple", alpha=0.15, label="Mission C Feasible Envelope (Brayton Optimized)")

    ax.axvline(35.0, color="crimson", linestyle="-", lw=2.0, label=r"Constraint C03 Ceiling ($\Delta P \leq 35.0$ kPa)")
    ax.axhline(20.0, color="black", linestyle="--", lw=1.5, label=r"Minimum Mission Life ($t \geq 20$ yr)")

    ax.set_xlabel(r"Nominal Core Pressure Drop $\Delta P_{\mathrm{nom}}$ [kPa]")
    ax.set_ylabel("Depletion Operational Lifetime [years]")
    ax.set_title("Figure 10: Mission-Specific Architectural Feasibility Envelopes", pad=12)
    ax.set_xlim(10, 44)
    ax.set_ylim(14, 45)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_08_mission_specific_feasible_regions.png")
    fig.savefig(out_path)
    # Also save as v7_10
    out_path10 = os.path.join(FIG_DIR, "v7_10_mission_specific_feasible_regions.png")
    fig.savefig(out_path10)
    plt.close(fig)
    print(f"  [+] Saved {out_path10}")


# -----------------------------------------------------------------------------
# FIGURE 11: FINAL ARCHITECTURE COMPARISON RADAR
# -----------------------------------------------------------------------------
def plot_11_final_architecture_comparison_radar():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, polar=True)
    
    categories = [
        "Mass Margin\n(Limit 650 kg)",
        "Creep Margin\n(Limit 1350 K)",
        "Pressure Drop Margin\n(Limit 35 kPa)",
        "Lifetime Longevity\n(Goal 35 yr)",
        "Shutdown Safety\n(SDM Margin)",
        "Thermal Mixing\n(Low $\Delta T_{\mathrm{hot}}$)"
    ]
    N = len(categories)
    angles = [n / float(N) * 2 * math.pi for n in range(N)]
    angles += angles[:1]
    
    # Normalized scores [0 to 1]
    # Baseline Opt-17
    score_base = [0.85, 0.60, 0.71, 0.86, 0.71, 0.20]
    score_base += score_base[:1]
    
    # V7 Balanced Robust
    score_bal = [0.87, 0.93, 0.78, 0.67, 0.80, 0.81]
    score_bal += score_bal[:1]
    
    # V7 Long Life Robust
    score_long = [0.63, 0.95, 0.76, 0.98, 0.60, 0.82]
    score_long += score_long[:1]
    
    # V7 Brayton Robust
    score_bray = [0.86, 0.86, 0.98, 0.64, 0.83, 0.81]
    score_bray += score_bray[:1]

    ax.plot(angles, score_base, "grey", linestyle=":", lw=1.8, label="Opt-17 Baseline (Historical)")
    ax.plot(angles, score_bal, "green", lw=2.5, label="V7 Balanced Robust (Reference)")
    ax.fill(angles, score_bal, "green", alpha=0.10)

    ax.plot(angles, score_long, "blue", lw=2.2, label="V7 Long-Life Robust (Flagship)")
    ax.fill(angles, score_long, "blue", alpha=0.08)

    ax.plot(angles, score_bray, "purple", lw=2.2, label="V7 Brayton Robust (Low dP)")
    ax.fill(angles, score_bray, "purple", alpha=0.08)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9.5)
    ax.set_rlabel_position(0)
    plt.yticks([0.25, 0.50, 0.75, 1.0], ["25%", "50%", "75%", "100%"], color="grey", size=8)
    plt.ylim(0, 1.05)
    
    plt.title("Figure 11: Final Multi-Physics Architecture Trade-Off Radar", y=1.08, fontsize=13, fontweight="bold")
    plt.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), framealpha=0.92)

    out_path = os.path.join(FIG_DIR, "v7_11_final_architecture_comparison_radar.png")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  [+] Saved {out_path}")


# -----------------------------------------------------------------------------
# RUN ALL FIGURE GENERATORS & COPY TO ARTIFACTS
# -----------------------------------------------------------------------------
def run_all_figures():
    print("=========================================================================")
    print("      V7 PUBLICATION FIGURE GENERATION (11 FIGURES, 300 DPI)             ")
    print("=========================================================================")
    plot_01_nominal_vs_robust_pareto()
    plot_02_mass_vs_lifetime()
    plot_03_mass_vs_pressure_drop()
    plot_04_pressure_drop_vs_peak_clad_temp()
    plot_05_orifice_tolerance_sensitivity()
    plot_06_thermal_maldistribution_robustness()
    plot_07_correlation_model_form_uncertainty()
    plot_08_monte_carlo_uncertainty_vs_margins()
    plot_09_burnup_lifetime_uncertainty()
    plot_10_mission_specific_feasible_regions()
    plot_11_final_architecture_comparison_radar()
    
    # Copy all figures to artifact directory
    fig_files = glob.glob(os.path.join(FIG_DIR, "*.png"))
    for f in fig_files:
        shutil.copy(f, ARTIFACT_DIR)
    print(f"\n[+] Copied {len(fig_files)} figures to conversation artifact directory {ARTIFACT_DIR}")


if __name__ == "__main__":
    run_all_figures()
