# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 14:06:33 2025

@author: lorena
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

##### add info below ####

file_path = r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\CWP_yield_wheat_2023.xlsx" # path to your file
output_dir = r"C:\Users\marta\OneDrive - Hydrosat\lorena_one_drive\cwp_project\output\corr_plots\test\rice_first_code_sec_total_final_test" #"#CWP\all_countries_per_year_rice_test_other_crops_together"

crops = ["cotton", "barley", "maize", "cotton", "wheat", "rice_rabi", "rice_kharif"]  # Add the crops you want
start_year = 2019 # first year of the plot
end_year = 2023 # last year of the plot
plt.rcParams["font.family"] = "Space Grotesk" # font Space Grotesk

##################

os.makedirs(output_dir, exist_ok=True)
df = pd.read_excel(file_path, engine="openpyxl")
years = list(range(start_year, end_year + 1))  

for crop in crops:
    df_crop = df[(df["crop_type"] == crop)] #& (df["aez_used"] == "country")] 

    x_min = 0
    x_max = df_crop["yield"].max()
    if crop.lower() in ["barley", "cotton", "wheat"]:
        x_interval = 500 # in the HyWater Global Analysis graphs, the x-axis interval is more spaced (3000), but I prefer a smaller one to make it easier to read the values.
    elif crop.lower() in ["maize"]: 
        x_interval = 1000 # adding a different interval in the plot because maize yield is much higher than the other crops
    else:
        x_interval = 500

    y_min = 0
    y_max = df_crop["CWP_kg_per_m3"].max()
    if crop.lower() in ["barley", "cotton"]:
        y_interval = 0.1
    elif crop.lower() in ["wheat", "maize"]:
        y_interval = 0.2
    else:
        y_interval = 0.1

    for year in years:
        df2 = df_crop[df_crop["year"] == year]

        if df2.empty:
            continue

        plt.figure(figsize=(8, 6))
        plt.scatter(df2["yield"], df2["CWP_kg_per_m3"], s=80, color="#8b9ada") # tried to use the same licac color as the in Hywater portal - Global Analysis

        for _, row in df2.iterrows():
            plt.text(
                row["yield"], row["CWP_kg_per_m3"],
                row["country"],
                fontsize=9, ha="right", va="bottom"
            )

        plt.xticks(np.arange(x_min, x_max + x_interval, x_interval))
        plt.yticks(np.arange(y_min, y_max + y_interval, y_interval))
        plt.xlim(x_min, x_max + x_interval)
        plt.ylim(y_min, y_max + y_interval)

        plt.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
        plt.xlabel("Yield (kg/ha)")
        plt.ylabel("Water Productivity (kg/m3)")
        plt.title(f"Water productivity for {crop} in {year}")
        plt.tight_layout()

        file_name = f"Water productivity for {crop} in {year}.png"
        save_path = os.path.join(output_dir, file_name)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()




