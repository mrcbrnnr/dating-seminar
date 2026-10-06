# Dating Team: Shiny app for the seminar “Dating Methods”

Interactive learning environment for ¹⁴C calibration, dendrochronology, coins (TPQ), written sources (TAQ),
stratigraphy and Bayesian modelling. Developed for a seminar at the University of Tübingen.
All example data are **fictional**.

**Open the app in the browser:** https://mrcbrnnr.github.io/dating-seminar/

The web version runs entirely in the browser (Shinylive, R as WebAssembly). No R installation and no server
are needed. The first load takes about 20 to 40 seconds. Recommended: laptop with a current version of
Chrome, Edge, Firefox or Safari.

## What the app can do

| Tab | Content |
|---|---|
| 1 · Calibration & model | Enter samples and layers or load a scenario; model with stratigraphy, TPQ, TAQ and dendro windows; switch assumptions on and off individually; comparison with the deposition from the ¹⁴C data only |
| 2 · Wiggle matching | Fit a simulated ring sequence to the IntCal20 curve, χ² test, fit curve |
| 3 · SPD & combination | Summed curve (different events) or combination (same event, weighted mean with χ² test) |
| 4 · OxCal code | Generates OxCal code from the model for comparison |
| 5 · Old-wood check | How much older is a sample than the youngest sample of its layer (or than the deposition)? |
| 6 · Measurement → age | Convert F14C to ¹⁴C age (and back) and calibrate |
| 7 · Dendro crossdating | Shift a sample against a **simulated** reference chronology (r, t-value, sign agreement), felling-year window from bark edge, sapwood or heartwood |

Scenarios in tab 1: “Fort Musterberg” (with a conflict in layer C and with the conflict resolved),
“Plateau example (Hallstatt period)”, “Old-wood exercise” and “Collective grave” (exercise after M. Hinz 2012).

Years may be entered as text in the tables (“480 BC”, “101 AD”). Internally years are counted astronomically as
in OxCal (year 0 = 1 BC). The selector at the top right switches the display between BC/AD and astronomical counting.

Column names in the tables: `id, layer, type, value1, value2` for samples and `layer, rank, tpq, taq` for layers
(the German column names of the German version are accepted as well).

## Repository

```
app/app.R                          the complete app
.github/workflows/deploy.yml       builds the web version and publishes it via GitHub Pages
README.md                          this file
```

With every change to `app/app.R` on the branch `main`, GitHub Actions rebuilds the site automatically
(a few minutes). Prerequisite: **Settings → Pages → Source: GitHub Actions**.
The IntCal20 calibration curve is downloaded from intcal.org during the build and placed next to the app.

## Run locally

```r
install.packages("shiny")
shiny::runApp("app")
```

The curve is read from `app/intcal20.14c`. If the file is missing, the app tries the package `rcarbon`
or a download from intcal.org locally.

## Limitations

- The model is a didactic Monte Carlo model (draw and reject), not a replacement for OxCal or other MCMC software.
  Results are similar but not identical.
- Not included: reservoir effects, Marine20 and SHCal20, outlier models (only in the OxCal export).
- The reference chronology in tab 7 is simulated and does not replace a real standard chronology.
- The ¹⁴C ages of the case study are chosen for illustration.

## Sources and acknowledgements

- IntCal20: Reimer et al. 2020, *Radiocarbon* 62, 725–757. Please observe the terms of use on intcal.org.
- Sapwood ranges (oak) in tab 7 after Tegel et al. 2022, *Frontiers in Ecology and Evolution* 10.
- Bayesian foundations and OxCal building blocks: Bronk Ramsey 2009, *Radiocarbon* 51.
- Ideas and exercise data (collective grave): materials by Martin Hinz for the seminar “Absolute Chronologie und
  Isotopenforschung” (2012), used with attribution.
