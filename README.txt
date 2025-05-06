# cna-inferer

A Python package to infer copy number alterations (CNAs) from single‑cell RNA‑seq data using sliding‑window aggregation and segmentation.

## 📂 Repository Structure

```

.
├── cna\_inferer/                   # Core Python module
│   ├── **init**.py
│   ├── annotation.py
│   ├── aggregation.py
│   ├── segmentation.py
│   └── call.py                    # Contains process\_and\_call\_cnas
├── data/
│   ├── Homo\_sapiens.GRCh38.104.gtf
│   └── PBMC\_simulated\_cnas\_041025.h5ad
├── notebooks/
│   ├── Task2A\_ReadDepth.ipynb
│   ├── Task2A\_Performance.ipynb
│   └── Task2B\_Simulation.ipynb
├── LICENSE
├── setup.py
└── README.md

````

> **Notes:**  
> - `cna_inferer/`: core library code  
> - `data/`: example GTF annotation & simulated `.h5ad`  
> - `notebooks/`: analysis notebooks for Task 2A/2B  
> - `setup.py`: installation setup  
> - `LICENSE`: MIT license  

---

## ⚙️ Installation

```bash
# Clone the repo
git clone https://github.com/chenming-pu/cna-inferer.git
cd cna-inferer

# (Optional) Create & activate a virtual environment
python3 -m venv venv && source venv/bin/activate

# Install editable package
pip install -e .
````

All dependencies (scanpy, pandas, numpy, pyensembl, gtfparse, polars, ruptures, joblib, etc.) will be installed via pip.

---

## 🚀 Quick Start

### 1. Python API

```python
from cna_inferer.main import process_and_call_cnas

adata = process_and_call_cnas(
    h5_file="data/PBMC_simulated_cnas_041025.h5ad",
    gtf_file="data/Homo_sapiens.GRCh38.104.gtf",
    window_size=150,
    z_thresh=1.5,
    min_bins=3,
    n_bkps=4,
    model="l2",
    n_jobs=4
)

# Inspect results
print("Event table:\n", adata.uns["cna_events"].head())
print("Per-cell calls:\n", adata.obs["cna_calls"][:3])
```

### 2. Command‑Line Interface (future)

> Once you add a `cli.py` and configure `entry_points` in `setup.py`, you could run:
>
> ```bash
> cna-infer --input data/PBMC_simulated_cnas_041025.h5ad \
>           --gtf data/Homo_sapiens.GRCh38.104.gtf \
>           --window 150 \
>           --z-thresh 1.5 \
>           --min-bins 3 \
>           --n-bkps 4 \
>           --model l2 \
>           --output results.h5ad
> ```

---

## 🔍 Examples & Notebooks

* **Task 2A**: `notebooks/Task2A_ReadDepth.ipynb` (read‑depth downsampling analysis)
* **Task 2A Performance**: `notebooks/Task2A_Performance.ipynb`
* **Task 2B**: `notebooks/Task2B_Simulation.ipynb`

---

## 🧪 Testing

```bash
pytest tests/
```

---

## 📄 License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

```
```
