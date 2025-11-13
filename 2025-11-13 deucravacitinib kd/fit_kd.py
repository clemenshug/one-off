from collections import defaultdict
import multiprocessing as mp
from multiprocessing import Queue, Process
import os
import json
import contextlib
import gzip
import pickle
import time
import sqlite3

import pytensor
import numpy as np
import pandas as pd
import pymc as pm
import pymc.math
import nutpie
import pystow
import tqdm
import matplotlib
import matplotlib.pyplot as plt
import synapseclient

matplotlib.use("Agg")

np.random.seed(42)

# pytensor.config.cxx = "/usr/bin/g++"
pytensor.config.cxx = "/usr/bin/clang++"


def model_predict(kd, hill_slope, doses):
    return 100 - 100 / (1 + (kd / doses) ** hill_slope)


# (Intercept) -0.0429160  0.0188197   -2.28   0.0458 *
# V2           0.1499235  0.0002977  503.52   <2e-16 ***
# From KINOMEscan General Introduction.pptx
def error_model(mu):
    # constrain to minimum of 1
    return pymc.math.maximum(0.15 * mu - 0.043, 1)


def simulated_data(kd, hill_slope, doses, error_std):
    # Model for data generation
    mu_true = model_predict(kd, hill_slope, doses)
    responses = np.random.normal(mu_true, scale=error_std)
    data_errors = np.full_like(responses, error_std)
    return doses, responses, data_errors


def create_model_template():
    """Create a PyMC model template that can handle varying numbers of data points."""
    with pm.Model() as model:
        # Priors for kd and hill_slope using lognormal (via normal on log)
        kd_log = pm.Normal("kd_log", mu=np.log(1e-6), sigma=3)
        kd = pm.Deterministic("kd", pm.math.exp(kd_log))

        hill_slope_log = pm.Normal("hill_slope_log", mu=0, sigma=0.5)
        hill_slope = pm.Deterministic("hill_slope", pm.math.exp(hill_slope_log))

        # Use mutable data containers that can be updated with different shapes
        doses_data = pm.Data("doses_data", np.array([1e-9]))
        responses_data = pm.Data("responses_data", np.array([0.0]))

        mu = model_predict(kd, hill_slope, doses_data)

        # Likelihood of observed data:
        sigma = error_model(mu)
        obs = pm.Normal("obs", mu=mu, sigma=sigma, observed=responses_data)

    return model


def train_model(doses, responses, compiled_sampler=None):
    """Train the model using a pre-compiled sampler to avoid recompilation.

    Args:
        doses: Array of dose values
        responses: Array of response values
        compiled_sampler: Pre-compiled nutpie sampler, if None will create new one

    Returns:
        trace: PyMC trace object
        compiled_sampler: The compiled sampler (for reuse)
    """
    if compiled_sampler is None:
        # Create the model template and compile it
        model = create_model_template()
        compiled_sampler = nutpie.compile_pymc_model(model)

    # Sample using the compiled sampler
    trace = nutpie.sample(
        compiled_sampler.with_data(doses_data=doses, responses_data=responses),
        draws=2000,
        tune=1000,
        chains=4,
        cores=1,
        progress_bar=False,
    )

    return trace, compiled_sampler


def process_data(df, okl_only=True):
    # Organize data by DiscoveRx Gene Symbol and then dose
    # into a dict of two lists
    data = {}
    for _, row in df.iterrows():
        compound_id = row["Compound Name"]
        gene = row["DiscoveRx Gene Symbol"]
        dose = row["Compound Concentration (nM)"]
        response = row["Percent Control"]
        if compound_id not in data:
            data[compound_id] = {}
        if gene not in data[compound_id]:
            data[compound_id][gene] = {"doses": [], "responses": []}

        data[compound_id][gene]["doses"].append(dose)
        data[compound_id][gene]["responses"].append(response)
    return data


deucravacitinib_kd_data = pd.read_csv(
    "2025-11-13 deucravacitinib kd/Deucravacitinib_KINOMEscan_100nM_10uM_raw.csv"
)

deucravacitinib_kd_dict = process_data(deucravacitinib_kd_data)

model_res_raw = {}
for compound_id, compound_data in deucravacitinib_kd_dict.items():
    sampler = None
    for target_id, target_data in compound_data.items():
        res, sampler = train_model(
            # convert to mols
            np.array(target_data["doses"]) * 1e-9,
            target_data["responses"],
            compiled_sampler=sampler,
        )
        model_res_raw[(compound_id, target_id)] = res

model_res_df = []
for (compound_id, target_id), trace in model_res_raw.items():
    summary = pm.summary(
        trace, var_names=["kd_log", "hill_slope_log"]
    ).rename_axis("metric").reset_index()

    summary["Compound Name"] = compound_id
    summary["DiscoveRx Gene Symbol"] = target_id
    model_res_df.append(summary)
model_res_df = pd.concat(model_res_df, axis=0)
model_res_df_wide = model_res_df[["Compound Name", "DiscoveRx Gene Symbol", "metric", "mean", "hdi_3%", "hdi_97%"]].pivot(
    index=["Compound Name", "DiscoveRx Gene Symbol"],
    columns="metric",
    values=["mean", "hdi_3%", "hdi_97%"],
)
model_res_df_wide.columns = [
    "_".join(col).strip() for col in model_res_df_wide.columns.to_flat_index()
]
model_res_df_wide.rename(
    columns={
        "mean_kd_log": "kd_log",
        "hdi_3%_kd_log": "kd_log_lower",
        "hdi_97%_kd_log": "kd_log_upper",
        "mean_hill_slope_log": "hill_slope_log",
        "hdi_3%_hill_slope_log": "hill_slope_log_lower",
        "hdi_97%_hill_slope_log": "hill_slope_log_upper",
    },
    inplace=True
)

model_res_df_wide["kd"] = np.exp(model_res_df_wide["kd_log"]) * 1e9  # Convert to nM
model_res_df_wide["kd_lower"] = np.exp(model_res_df_wide["kd_log_lower"]) * 1e9
model_res_df_wide["kd_upper"] = np.exp(model_res_df_wide["kd_log_upper"]) * 1e9

model_res_df_wide.reset_index().to_csv(
    "2025-11-13 deucravacitinib kd/deucravacitinib_kd_estimates.csv",
    index=False,
)

model_res_df_wide.sort_values("kd")[["kd", "kd_lower", "kd_upper"]]



def plot_results(trace, dose_range, doses, responses, fname):
    # Extract posterior samples of Kd and HillSlope
    # These have shape (chains, draws), so flatten them into a single dimension.
    kd_samples = trace.posterior["kd"].values.reshape(-1)
    hill_slope_samples = trace.posterior["hill_slope"].values.reshape(-1)

    # For each posterior sample, compute the predicted response at each dose in dose_range
    # pred_matrix will have shape (n_samples, len(dose_range))
    pred_matrix = []
    for kd, hill_slope in zip(kd_samples, hill_slope_samples):
        mu = model_predict(kd, hill_slope, dose_range)
        pred_matrix.append(mu)
    pred_matrix = np.array(pred_matrix)  # shape: (n_samples, n_dose_points)

    # Compute mean and 95% credible intervals
    pred_mean = np.mean(pred_matrix, axis=0)
    pred_lower = np.percentile(pred_matrix, 2.5, axis=0)
    pred_upper = np.percentile(pred_matrix, 97.5, axis=0)

    kd_mean = np.mean(kd_samples)
    kd_std = np.std(kd_samples)

    # Use explicit fig, ax API and ensure figures are closed to avoid memory leaks
    fig, ax = plt.subplots(figsize=(6, 4))

    # Plot the measurements
    ax.scatter(np.array(doses)*1e-9, responses, color="black", label="Observed Data")

    # Plot the posterior predictive mean and credible interval over the specified dose range
    ax.plot(dose_range, pred_mean, color="red", label="Posterior Predictive Mean")
    ax.fill_between(
        dose_range,
        pred_lower,
        pred_upper,
        color="red",
        alpha=0.3,
        label="95% Credible Interval",
    )

    # Put the Kd estimate in the plot as a solid vertical line with dashed lines for +/- 1 SD
    ax.axvline(kd_mean, color="green", linestyle="-", label="Estimated Kd")
    ax.axvline(
        kd_mean - kd_std, color="green", linestyle="--", label="Estimated Kd +/- 1 SD"
    )
    ax.axvline(kd_mean + kd_std, color="green", linestyle="--")

    ax.set_xscale("log")
    ax.set_ylim(-5, 105)
    ax.set_xlabel("Dose")
    ax.set_ylabel("Response")
    ax.set_title("Dose-Response with Posterior Predictive")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fname)
    plt.close(fig)

plot_results(
    model_res_raw[("Deucravacitinib", "TYK2(JH2domain-pseudokinase)")],
    dose_range=np.logspace(-10, -5, 100),
    doses=deucravacitinib_kd_dict["Deucravacitinib"]["TYK2(JH2domain-pseudokinase)"]["doses"],
    responses=deucravacitinib_kd_dict["Deucravacitinib"]["TYK2(JH2domain-pseudokinase)"]["responses"],
    fname="2025-11-13 deucravacitinib kd/deucravacitinib_tyk2_dose_response.png",
)
