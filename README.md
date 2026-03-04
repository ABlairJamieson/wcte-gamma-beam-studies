
# WCTE Gamma Beam Studies

Analysis tools for studying **300 MeV gamma beam simulations in WCSim** for the Water Cherenkov Test Experiment (WCTE).

This project analyzes simulated events to understand how **pion photoproduction** can be distinguished from pure electromagnetic gamma showers using Cherenkov detector observables.

The code processes WCSim ROOT output and produces diagnostic plots comparing events **with and without pion production**.

---

# Physics Motivation

Gamma interactions in water at ~300 MeV are dominated by electromagnetic processes:

- Pair production
- Compton scattering
- Electromagnetic shower development

However, **photonuclear interactions** can occasionally produce pions. These events are rare but important to understand because they may mimic other physics signals or affect detector calibration.

This analysis explores observables that can separate:

- **pure EM gamma showers**
- **gamma → pion photoproduction events**

using Cherenkov detector information such as:

- total charge
- number of PMT hits
- charge-per-hit
- spatial hit patterns
- forward Cherenkov cone structure

---

# Features

The analysis currently includes several studies:

### Basic event observables
- Total number of PMT hits
- Total charge
- Charge per hit
- Comparison of **pion vs non-pion events**

### 2D distributions
- Charge vs number of hits
- Separate overlays for pion and non-pion events

### Vertex studies
- Truth interaction vertex distribution
- Comparison of shower start vs pion production points

### Cherenkov cone analysis
A study of light inside vs outside the Cherenkov cone (~42° in water).

Variables computed:

- \(N_{in}, N_{out}\)
- \(Q_{in}, Q_{out}\)

Derived observables:

- \(N_{out}/N_{tot}\)
- \(Q_{out}/Q_{tot}\)

2D distributions are also produced for:

- \(Q_{in}\) vs \(N_{in}\)
- \(Q_{out}\) vs \(N_{out}\)
- \(Q_{in}\) vs \(N_{out}\)
- \(Q_{out}\) vs \(N_{in}\)

These help visualize how pion events populate different regions of phase space.

---

# Repository Structure

```

BeamMC_studies/
│
├── main.cc
├── Makefile
│
├── include/
│   ├── Study.h
│   ├── Types.h
│   ├── Utils.h
│   └── Studies_*.h
│
├── src/
│   ├── Run.cc
│   ├── Utils.cc
│   ├── Studies_BasicSpectra.cc
│   ├── Studies_QperHit.cc
│   ├── Studies_QvsN_2D.cc
│   ├── Studies_VtxXZ.cc
│   └── Studies_Cone42.cc
│
└── plots/

```

The code is organized so that **each physics study is implemented in its own module**, making it easy to add additional analyses.

---

# Building

The code requires:

- ROOT
- WCSim
- libWCSimRoot

Inside the analysis container environment:

```

make

```

This produces the executable:

```

beam_mc_studies

```

---

# Running

Example:

```

./beam_mc_studies ../WCSim_full/build/src/wcsim_gamma_300_100000.root

```

Optional configuration:

```

--cone-vertex=truth

```

This controls whether the Cherenkov cone study uses:

- a nominal vertex position, or
- the **truth interaction vertex**.

---

# Output

The program generates PDF plots including:

- NDigi distributions
- charge-per-hit distributions
- Q vs N 2D overlays
- vertex distributions
- Cherenkov cone observables

These are saved in the working directory.

---

# Development Notes

The code is designed to:

- process large WCSim samples efficiently
- allow modular addition of new physics studies
- keep the event loop centralized while delegating analysis tasks to study modules

This structure makes it straightforward to add additional observables or event classifiers.

---

# Acknowledgements

The software structure and many analysis utilities in this repository were developed with **AI-assisted programming support (ChatGPT)** under the direction of the repository author.

Design decisions, physics analysis strategy, and validation were performed by the repository author.

---

# License

For research and educational use.
```
