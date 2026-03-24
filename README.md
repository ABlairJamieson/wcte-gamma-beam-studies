
# MC HAMMER 

**Monte Carlo Histogramming, Analysis, and Modular Modeling for Event Research**

A lightweight, modular framework for analyzing **WCSim ROOT data**.  
Designed for rapid physics studies with minimal boilerplate — just *hammer the data*.

<img width="766" height="303" alt="image" src="https://github.com/user-attachments/assets/2cb8d005-7ad1-42ed-a3d7-694a6ff76ceb" />


---

## Overview

MC HAMMER provides a simple structure for:

- looping over WCSim events
- extracting detector and truth information
- implementing modular analysis “studies”
- producing ROOT-based histograms and plots

It is intended as a **starting point for students and researchers** working with water Cherenkov detector simulations.

---

## Example Application (Current Study)

This repository currently includes an example analysis of:

 **300 MeV gamma beam simulations in WCSim (WCTE context)**

The goal is to distinguish:

- **pure electromagnetic showers**
- **gamma -> pion photoproduction events**

using Cherenkov detector observables.

---

## Physics Motivation

Gamma interactions in water at ~300 MeV are dominated by:

- Pair production
- Compton scattering
- Electromagnetic shower development

However, **photonuclear interactions** can produce pions.

These events are:
- rare
- topologically different
- important for backgrounds and calibration

MC HAMMER is used here to explore observables that separate:

- EM-only events
- events containing pions

---

## Framework Design

### Modular Study System

Each analysis is implemented as a **study module**:

```cpp
struct MyStudy : IStudy {
  void FillEvent(...) override { ... }
  void WritePlots() override { ... }
};
```

Studies are registered in the run loop:

```cpp
studies.emplace_back(std::make_unique<MyStudy>(...));
```

This allows:
- clean separation of analyses
- easy extension
- multiple studies running in a single pass over data

### Core Components

- **Event loop** (centralized, efficient)
- **Geometry cache** (PMT positions, detector info)
- **Truth utilities** (event classification, vertex extraction)
- **Plot helpers** (overlays, normalization)

---

## Included Studies

### Basic observables
- number of PMT hits
- total charge
- charge per hit
- pion vs non-pion comparisons

### 2D distributions
- charge vs number of hits
- overlay comparisons

### Vertex studies
- truth interaction vertex distributions

### Cone-based topology study

Cherenkov cone analysis using a configurable angle.

Variables:
- \(N_{in}, N_{out}\)
- \(Q_{in}, Q_{out}\)

Derived:
- \(N_{out}/N_{tot}\)
- \(Q_{out}/Q_{tot}\)

2D phase space:
- \(Q_{in}\) vs \(N_{in}\)
- \(Q_{out}\) vs \(N_{out}\)
- cross-combinations

Used to identify broader angular light distributions in pion events.

---

## Repository Structure

```text
mc-hammer/
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
│   ├── Studies_*.cc
│
└── plots/
```

---

## Building

Requires:
- ROOT
- WCSim
- libWCSimRoot

Inside your environment:

```bash
make
```

This produces the executable:

```bash
beam_mc_studies
```

---

## Running

Example:

```bash
./beam_mc_studies input.root
```

Optional arguments:

```bash
--cone-angle=42
--cone-vertex=truth
```

---

## Output

Produces PDF plots including:

- 1D distributions
- normalized overlays
- 2D phase space plots
- cone-based observables

These are saved in the working directory.

---

## Adding a New Study

1. Create a new study file, for example `Studies_MyStudy.h/.cc`

```cpp
struct MyStudy : IStudy {
  void FillEvent(...) override {
    // event-level analysis
  }

  void WritePlots() override {
    // plotting
  }
};
```

2. Register it in `Run.cc`:

```cpp
studies.emplace_back(std::make_unique<MyStudy>(...));
```

3. Add the corresponding `.cc` file to the Makefile.

That is all that is required for the new study to run automatically in the main event loop.

---

## Development Notes

The code is designed to:

- process large WCSim samples efficiently
- allow modular addition of new physics studies
- keep the event loop centralized while delegating analysis tasks to study modules

This structure makes it straightforward to add additional observables or event classifiers.

---

## Future Directions

This framework can be extended to:

- particle ID studies
- timing-based analyses (prompt vs delayed light)
- Michel electron tagging
- reconstruction validation
- real detector data

---

## Acknowledgements

The software structure and many utilities were developed with **AI-assisted programming support (ChatGPT)** under the direction of the repository author.

Physics design, validation, and interpretation were performed by the author.

---

## License

For research and educational use.

