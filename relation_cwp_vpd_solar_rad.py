# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 09:37:18 2025

@author: lorena
"""
    
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, spearmanr #,pearsonr


fn = r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\CWP_yield_wheat_2023.xlsx" #_no_barley
df = pd.read_excel(fn, engine="openpyxl")

 
### separate per unique crop:
    
crop_types = df["crop_type"].unique()
df_rice = df[df["crop_descrip"] == "rice"]
variables = ["CWP_kg_per_m3", "vpd_kPa_ndvi", "net_solar_rad_J_m2"]

for crop in crop_types:
    
    globals()[f"df_cwp_{crop}"] = df[df["crop_type"] == crop][["CWP_kg_per_m3"]]
    globals()[f"df_vpd_{crop}"] = df[df["crop_type"] == crop][["vpd_kPa_ndvi"]]
    globals()[f"df_sol_{crop}"] = df[df["crop_type"] == crop][["net_solar_rad_J_m2"]]


df_cwp_rice = df_rice[["CWP_kg_per_m3"]]
df_vpd_rice = df_rice[["vpd_kPa_ndvi"]]
df_sol_rice = df_rice[["net_solar_rad_J_m2"]]



###### CHECK IF THE DATA HAS A NORMAL DISTRIBUTION ##############

# source: https://www.statology.org/normality-test-python/
## shapiro wilk test: If the p-value of the test is greater than α = .05, then the data is assumed to be normally distributed.



results = []

# Shapiro-Wilk test for each crop type and variable
for crop in crop_types:
    crop_df = df[df["crop_type"] == crop]
    for var in variables:
        data = crop_df[var].dropna()
        if len(data) >= 3:  # Shapiro test requires at least 3 data points
            stat, pval = shapiro(data)
            normality = "normal" if pval > 0.05 else "not normal"
            results.append({
                "crop_type": crop,
                "variable": var,
                "p_value": pval,
                "normality": normality
            })

summary_df = pd.DataFrame(results)

summary_df

results_rice = []



if not df_rice.empty:
    for var in variables:
        data = df_rice[var].dropna()
        if len(data) >= 3:  # Shapiro requires at least 3 values
            stat, pval = shapiro(data)
            normality = "normal" if pval > 0.05 else "not normal"
            results_rice.append({
                "crop_type": "rice",
                "variable": var,
                "p_value": pval,
                "normality": normality
            })

rice_summary_df = pd.DataFrame(results_rice)

crops_cwp_vpd_sol_distrib = pd.concat([summary_df, rice_summary_df], ignore_index=True)
# crops_cwp_vpd_sol_distrib_transposed = crops_cwp_vpd_sol_distrib.transpose()
crops_cwp_vpd_sol_distrib.to_excel(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\crops_cwp_vpd_sol_normallity_distrib.xlsx", index = False) #_no_barley



##################################################### IF DATA IS NORMAL AND LINEAR USE PEARSON, IF DATA IS NOT NORMALLY DISTRIBUTED AND NOT LINEAR USE SPEARMAN! ###############

# ############## PEARSON ###################
# # source: https://www.geeksforgeeks.org/python/python-pearson-correlation-test-between-two-variables/

################ SPEARMAN 

# Source: https://www.geeksforgeeks.org/data-science/spearmans-rank-correlation/

# It is designed to capture monotonic relationships between variables. 
# Monotonic relation measures the effect of change in one variable on another variable


"""
What is Spearman's Correlation
Spearman's Rank Correlation is a statistical measure of the strength and direction of the 
monotonic relationship between two continuous variables. Therefore, these attributes are ranked 
or put in the order of their preference. It is denoted by the symbol "rho" (ρ) and can take values 
between -1 to +1. A positive value of rho indicates that there exists a positive relationship between 
the two variables, while a negative value of rho indicates a negative relationship. A rho value of 0 indicates 
no association between the two variables.

""" 


#### USING SPEARMAN CORRELATION FOR ALL CROPS: 

results_correlation = []
crop_types = crops_cwp_vpd_sol_distrib["crop_type"].unique()


for crop in crop_types:
    
    crop_df = df[df["crop_type"] == crop]
  
    cwp = crop_df["CWP_kg_per_m3"].dropna()
    vpd = crop_df["vpd_kPa_ndvi"].dropna()
    sol = crop_df["net_solar_rad_J_m2"].dropna()
   
    cwp_vpd = crop_df[["CWP_kg_per_m3", "vpd_kPa_ndvi"]].dropna()
    cwp_sol = crop_df[["CWP_kg_per_m3", "net_solar_rad_J_m2"]].dropna()

    # Get normality status (actually not used because I am using spearman for all the crops)
    cwp_norm = crops_cwp_vpd_sol_distrib[(crops_cwp_vpd_sol_distrib["crop_type"] == crop) & (crops_cwp_vpd_sol_distrib["variable"] == "CWP_kg_per_m3")]["normality"].values[0]
    vpd_norm = crops_cwp_vpd_sol_distrib[(crops_cwp_vpd_sol_distrib["crop_type"] == crop) & (crops_cwp_vpd_sol_distrib["variable"] == "vpd_kPa_ndvi")]["normality"].values[0]
    sol_norm = crops_cwp_vpd_sol_distrib[(crops_cwp_vpd_sol_distrib["crop_type"] == crop) & (crops_cwp_vpd_sol_distrib["variable"] == "net_solar_rad_J_m2")]["normality"].values[0]

    # CWP vs VPD
    if len(cwp_vpd) >= 2:
        corr, pval = spearmanr(cwp_vpd["CWP_kg_per_m3"], cwp_vpd["vpd_kPa_ndvi"])
        method = "spearman"
        results_correlation.append({
            "crop_type": crop,
            "variable_pair": "CWP vs VPD",
            "correlation": corr,
            "p_value": pval,
            "method": method
        })
    
    # CWP vs Solar Radiation
    if len(cwp_sol) >= 2:
        corr, pval = spearmanr(cwp_sol["CWP_kg_per_m3"], cwp_sol["net_solar_rad_J_m2"])
        method = "spearman"
        results_correlation.append({
            "crop_type": crop,
            "variable_pair": "CWP vs Solar Radiation",
            "correlation": corr,
            "p_value": pval,
            "method": method
        })
    
    correlation_summary_df = pd.DataFrame(results_correlation)
    correlation_summary_df
        

#### SPEARMAN CORRELATION FOR RICE ONLY:

rice_kharif = df[(df["crop_descrip"] == "rice") & (df["season"] == "kharif")]
rice_rabi = df[(df["crop_descrip"] == "rice") & (df["season"] == "rabi")]


results_correlation_rice = []


def compute_correlation(data, crop_type, var1, var2, label):
    subset = data[[var1, var2]].dropna()
    if len(subset) >= 2:
        corr, pval = spearmanr(subset[var1], subset[var2])
        results_correlation.append({
            "crop_type": crop_type,
            "variable_pair": label,
            "correlation": corr,
            "p_value": pval,
            "method": "spearman"
        })

# correlations for rice_kharif
compute_correlation(rice_kharif, "rice_kharif", "CWP_kg_per_m3", "vpd_kPa_ndvi", "CWP vs VPD")
compute_correlation(rice_kharif, "rice_kharif", "CWP_kg_per_m3", "net_solar_rad_J_m2", "CWP vs Solar Radiation")

# correlations for rice_rabi
compute_correlation(rice_rabi, "rice_rabi", "CWP_kg_per_m3", "vpd_kPa_ndvi", "CWP vs VPD")
compute_correlation(rice_rabi, "rice_rabi", "CWP_kg_per_m3", "net_solar_rad_J_m2", "CWP vs Solar Radiation")

correlation_summary_df_rice = pd.DataFrame(results_correlation_rice)        
correlation_crops_df = pd.concat([correlation_summary_df,correlation_summary_df_rice], ignore_index=True)     
correlation_crops_df.to_excel(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\correlation_p_value_crops_only_spearman.xlsx", index = False) # _no_barley  



####### PLOTS ##########
    
        
##### CWP X VPD PLOT ALL CROPS AND YEARS    


colors = ["#FF5733",  # Vibrant Red-Orange
          "black",  # Bright Lime Green "#33FF57"
          "#3357FF",  # Vivid Blue
          "#FF33A8",  # Hot Pink
          "#FFD433",  # Bright Yellow
          "#24b2ac"]  # Electric Cyan


unique_crops = df["crop_type"].unique()
palette = dict(zip(unique_crops, colors[:len(unique_crops)]))
sns.set(style="whitegrid")
fig, axes = plt.subplots(1, 1, figsize=(16, 8))

# Scatter plot: CWP vs VPD
sns.scatterplot(data=df, x="vpd_kPa_ndvi", y="CWP_kg_per_m3", hue="crop_type", palette=palette, ax=axes, s=60)
axes.set_title("CWP vs VPD by Crop Type")
axes.set_xlabel("VPD (kPa)")
axes.set_ylabel("CWP (kg/m³)")

plt.tight_layout()
plt.savefig(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\cwp_vpd.png") # _no_barley
plt.close()
# plt.show()


#### CWP X SNSSR PLOT ALL CROPS AND YEARS    


fig, axes = plt.subplots(1, 1, figsize=(16, 8))

# Scatter plot: CWP vs Solar Radiation
sns.scatterplot(data=df, x="net_solar_rad_J_m2", y="CWP_kg_per_m3", hue="crop_type", palette=palette, ax=axes, s=60)
axes.set_title("CWP vs SNSSR by Crop Type")
axes.set_xlabel("SNSSR (J/m²)")
axes.set_ylabel("CWP (kg/m³)")

plt.tight_layout()
plt.savefig(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\cwp_sol_rad.png")
plt.close()
# plt.show()


####### PLOT CORRELATION CWP X VPD AND CWP X SNSSR PER CROP TYPE 

import numpy as np    
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set font
plt.rcParams["font.family"] = "Space Grotesk"
plt.rcParams["font.sans-serif"] = ["Space Grotesk", "sans-serif"]


# Create a DataFrame (replace None with actual correlation values if available)
correlation_df = correlation_summary_df


####### PLOT CORRELATION CWP X VPD  PER CROP TYPE 

# Create a single scatter plot for the available variable pair
# Set the style for seaborn
sns.set(style="whitegrid")
fig, ax = plt.subplots(figsize=(10, 5))
sns.scatterplot(data=correlation_df[correlation_df["variable_pair"] == "CWP vs VPD"],
                x="crop_type", y="correlation", ax=ax, color="#8b9ada")
ax.set_title("Correlation: CWP vs VPD")
ax.tick_params(axis='x', rotation=45, labelsize=10, pad=-5)
ax.set_ylim(-1, 1)
ax.set_yticks(np.arange(-1, 1.1, 0.2))  # include 1.0

plt.tight_layout()
plt.savefig(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\corr_cwp_vpd_per_crop.png") # _no_barley
plt.close()
# plt.show()

####### PLOT CORRELATION CWP X SNSSR PER CROP TYPE 

# Create a single scatter plot for the available variable pair
# Set the style for seaborn
sns.set(style="whitegrid")
fig, ax = plt.subplots(figsize=(10, 5))
sns.scatterplot(data=correlation_df[correlation_df["variable_pair"] == "CWP vs Solar Radiation"],
                x="crop_type", y="correlation", ax=ax, color="#8b9ada")
ax.set_title("Correlation: CWP vs SNSSR")
ax.tick_params(axis='x', rotation=45, labelsize=10, pad=-5)
ax.set_ylim(-1, 1)
ax.set_yticks(np.arange(-1, 1.1, 0.2))  # include 1.0

plt.tight_layout()
plt.savefig(r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\corr_cwp_snssr_per_crop.png") # _no_barley
plt.close()
# plt.show()

# scatter plots for each crop and variable pair
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))
sns.scatterplot(data=correlation_df[correlation_df["variable_pair"] == "CWP vs VPD"],
                x="crop_type", y="correlation", ax=axes[0], color="#8b9ada")
axes[0].set_title("Correlation: CWP vs VPD")
axes[0].tick_params(axis='x', rotation=45, labelsize = 10, pad = -5)

axes[0].set_ylim(-1, 1)
axes[0].set_yticks(np.arange(-1, 1, 0.2))

axes[1].set_title("Correlation: CWP vs SNSSR") # col "net_solar_rad_J_m2"
axes[1].tick_params(axis='x', rotation=45, labelsize = 10, pad = -5)

axes[1].set_ylim(-1, 1)
axes[1].set_yticks(np.arange(-1, 1, 0.2))









