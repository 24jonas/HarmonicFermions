import json

notebook_path = "/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/JonasPolyConv-opt.ipynb"

with open(notebook_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

new_source_cell3 = [
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "from joblib import Parallel, delayed\n",
    "\n",
    "def make_tau_grid(tau_start, tau_end, tau_step):\n",
    "    if tau_step <= 0:\n",
    "        raise ValueError(\"tau_step must be positive.\")\n",
    "    if tau_end < tau_start:\n",
    "        raise ValueError(\"tau_end must be >= tau_start.\")\n",
    "\n",
    "    taus = []\n",
    "    tau = tau_start\n",
    "    while tau <= tau_end + 1e-15:\n",
    "        taus.append(round(tau, 12))\n",
    "        tau += tau_step\n",
    "    return taus\n",
    "\n",
    "def finite_difference_user_sign(x, y):\n",
    "    if len(x) != len(y):\n",
    "        raise ValueError(\"x and y must have the same length.\")\n",
    "    if len(x) < 2:\n",
    "        raise ValueError(\"Need at least 2 points.\")\n",
    "\n",
    "    x_fd = x[:-1]\n",
    "    y_fd = []\n",
    "    for i in range(len(x) - 1):\n",
    "        dx = x[i + 1] - x[i]\n",
    "        y_fd.append((y[i] - y[i + 1]) / dx)\n",
    "    return np.array(x_fd), np.array(y_fd)\n",
    "\n",
    "def compute_logZ_worker(tau, N, n, d, logZ_kwargs):\n",
    "    return fermion_logZ_numeric(tau=tau, N=N, n=n, d=d, **logZ_kwargs)\n",
    "\n",
    "def compute_logZ_worker_all(tau, N, n, d, logZ_kwargs):\n",
    "    logZ_kwargs['return_all'] = True\n",
    "    _, logZ_list, _ = fermion_logZ_numeric(tau=tau, N=N, n=n, d=d, **logZ_kwargs)\n",
    "    return logZ_list\n",
    "\n",
    "def plot_logZ_and_fd_multiN(tau_start, tau_end, tau_step, n, d, N_list,\n",
    "                            show_logZ=True, show_fd=True, store_all_n=False, save_filename=\"plot_data\", **logZ_kwargs):\n",
    "    taus = np.array(make_tau_grid(tau_start, tau_end, tau_step))\n",
    "    results = {}\n",
    "\n",
    "    if store_all_n:\n",
    "        results_fd_all_n = []\n",
    "        results_logZ_all_n = []\n",
    "\n",
    "    for N in N_list:\n",
    "        if store_all_n:\n",
    "            logZ_kwargs['return_all'] = True\n",
    "            logZ_vals = Parallel(n_jobs=-1)(\n",
    "                delayed(compute_logZ_worker_all)(tau, N, n, d, logZ_kwargs)\n",
    "                for tau in taus\n",
    "            )\n",
    "            logZ_vals_arr = np.array(logZ_vals).T \n",
    "            \n",
    "            fd_vals_arr = []\n",
    "            for i in range(n + 1):\n",
    "                tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals_arr[i])\n",
    "                fd_vals_arr.append(fd_vals)\n",
    "            fd_vals_arr = np.array(fd_vals_arr)\n",
    "            \n",
    "            results_logZ_all_n.append(logZ_vals_arr)\n",
    "            results_fd_all_n.append(fd_vals_arr)\n",
    "            \n",
    "            results[N] = {\n",
    "                \"taus\": taus,\n",
    "                \"logZ\": logZ_vals_arr[n],\n",
    "                \"tau_fd\": tau_fd,\n",
    "                \"fd\": fd_vals_arr[n],\n",
    "            }\n",
    "        else:\n",
    "            logZ_kwargs['return_all'] = False\n",
    "            logZ_vals = Parallel(n_jobs=-1)(\n",
    "                delayed(compute_logZ_worker)(tau, N, n, d, logZ_kwargs)\n",
    "                for tau in taus\n",
    "            )\n",
    "            tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals)\n",
    "    \n",
    "            results[N] = {\n",
    "                \"taus\": taus,\n",
    "                \"logZ\": logZ_vals,\n",
    "                \"tau_fd\": tau_fd,\n",
    "                \"fd\": fd_vals,\n",
    "            }\n",
    "\n",
    "    if show_logZ:\n",
    "        plt.figure(figsize=(8, 5))\n",
    "        for N in N_list:\n",
    "            plt.plot(results[N][\"taus\"], results[N][\"logZ\"], marker=\".\", label=f\"N={N}\")\n",
    "        plt.xlabel(\"tau\")\n",
    "        plt.ylabel(\"log Z_n\")\n",
    "        plt.title(f\"log Z_n vs tau (n={n}, d={d})\")\n",
    "        plt.grid(True)\n",
    "        plt.legend()\n",
    "        plt.show()\n",
    "\n",
    "    if show_fd:\n",
    "        plt.figure(figsize=(8, 5))\n",
    "        for N in N_list:\n",
    "            plt.plot(results[N][\"tau_fd\"], results[N][\"fd\"], marker=\".\", label=f\"N={N}\")\n",
    "        plt.xlabel(\"tau\")\n",
    "        plt.ylabel(\"(log Z(tau) - log Z(tau + dt)) / dt\")\n",
    "        plt.title(f\"Finite difference vs tau (n={n}, d={d})\")\n",
    "        plt.grid(True)\n",
    "        plt.legend()\n",
    "        plt.show()\n",
    "\n",
    "    if store_all_n:\n",
    "        arr_logZ = np.array(results_logZ_all_n)[:, 1:, :].transpose(1, 0, 2)\n",
    "        arr_fd = np.array(results_fd_all_n)[:, 1:, :].transpose(1, 0, 2)\n",
    "        \n",
    "        np.save(f\"{save_filename}_logZ_all_n.npy\", arr_logZ)\n",
    "        np.save(f\"{save_filename}_fd_all_n.npy\", arr_fd)\n",
    "        np.save(f\"{save_filename}_taus.npy\", taus)\n",
    "        np.save(f\"{save_filename}_tau_fd.npy\", tau_fd)\n",
    "        np.save(f\"{save_filename}_beads.npy\", np.array(N_list))\n",
    "        \n",
    "        print(f\"Data successfully saved to '{save_filename}_fd_all_n.npy' (shape {arr_fd.shape})\")\n",
    "        print(f\"Axes mapping: n (axis 0), beads (axis 1) -> '{save_filename}_beads.npy', taus -> '{save_filename}_tau_fd.npy'\")\n",
    "    else:\n",
    "        df = pd.DataFrame(results[N_list[0]])\n",
    "        df.to_csv(f\"{save_filename}_results.csv\", index=False)\n",
    "\n",
    "    return results\n"
]

new_source_cell4 = [
    "results = plot_logZ_and_fd_multiN(\n",
    "    tau_start=0.1,\n",
    "    tau_end=10.0,\n",
    "    tau_step=0.5,\n",
    "    n=7,\n",
    "    d=2,\n",
    "    N_list=[1, 4, 8],\n",
    "    show_logZ=False,\n",
    "    show_fd=True,\n",
    "    store_all_n=True,\n",
    "    save_filename=\"plot_data\"\n",
    ")\n"
]

new_source_cell5 = [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import pandas as pd\n",
    "\n",
    "save_filename = \"plot_data\"  \n",
    "\n",
    "try:\n",
    "    arr_logZ = np.load(f\"{save_filename}_logZ_all_n.npy\")\n",
    "    arr_fd = np.load(f\"{save_filename}_fd_all_n.npy\")\n",
    "    taus = np.load(f\"{save_filename}_taus.npy\")\n",
    "    tau_fd = np.load(f\"{save_filename}_tau_fd.npy\")\n",
    "    beads = np.load(f\"{save_filename}_beads.npy\")\n",
    "    \n",
    "    print(f\"Loaded arrays with shapes:\")\n",
    "    print(f\"logZ_all_n: {arr_logZ.shape}\")\n",
    "    print(f\"fd_all_n: {arr_fd.shape}\")\n",
    "    print(f\"taus: {taus.shape}\")\n",
    "    print(f\"tau_fd: {tau_fd.shape}\")\n",
    "    print(f\"beads: {beads.shape}\")\n",
    "    \n",
    "    target_n = 7\n",
    "    target_bead = 4\n",
    "    \n",
    "    n_idx = target_n - 1\n",
    "    bead_idx = np.where(beads == target_bead)[0][0]\n",
    "    \n",
    "    print(f\"\\n Extracting data for n={target_n} (index {n_idx}) and bead={target_bead} (index {bead_idx})\")\n",
    "    \n",
    "    logZ_vals = arr_logZ[n_idx, bead_idx, :]\n",
    "    fd_vals = arr_fd[n_idx, bead_idx, :]\n",
    "    \n",
    "    df_extracted = pd.DataFrame({\n",
    "        \"tau_fd\": tau_fd,\n",
    "        \"Energy_FD\": fd_vals\n",
    "    })\n",
    "    export_filename = f\"extracted_n{target_n}_bead{target_bead}_polyconv.csv\"\n",
    "    df_extracted.to_csv(export_filename, index=False)\n",
    "    print(f\"Successfully exported to {export_filename}\\n \")\n",
    "    \n",
    "    plt.figure(figsize=(8, 5))\n",
    "    plt.plot(tau_fd, fd_vals, label=f\"N={target_bead}\", marker=\".\")\n",
    "    \n",
    "    plt.xlabel(\"tau\")\n",
    "    plt.ylabel(\"(log Z(tau) - log Z(tau + dt)) / dt\")\n",
    "    plt.title(f\"Extracted Finite Difference Data: n={target_n}, Bead N={target_bead}\")\n",
    "    plt.legend()\n",
    "    plt.grid(True)\n",
    "    plt.tight_layout()\n",
    "    plt.show()\n",
    "\n",
    "except FileNotFoundError:\n",
    "    print(f\"Error: Could not find the files starting with '{save_filename}'.\")\n",
    "    print(\"Please make sure you have run plot_logZ_and_fd_multiN with store_all_n=True first!\")\n"
]

# The cells in JonasPolyConv-opt currently:
# 0: colab drive mount
# 1: !ls
# 2: fermion_logZ_numeric
# 3: fermion_logZ_numeric test
# 4: plot_logZ_and_fd_multiN
# 5: results = plot_logZ_and_fd_multiN(...)

for i, cell in enumerate(nb['cells']):
    src = "".join(cell['source'])
    if "def plot_logZ_and_fd_multiN(" in src:
        nb['cells'][i]['source'] = new_source_cell3
    elif "results = plot_logZ_and_fd_multiN(" in src:
        nb['cells'][i]['source'] = new_source_cell4

# Check if we already have 6 cells (the extraction cell)
has_extraction_cell = False
for cell in nb['cells']:
    if "target_bead" in "".join(cell['source']) and "target_n" in "".join(cell['source']):
        has_extraction_cell = True
        break

if not has_extraction_cell:
    new_cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": new_source_cell5
    }
    nb['cells'].append(new_cell)

with open(notebook_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("JonasPolyConv-opt updated with store_all_n and extraction cell successfully.")
