# A dynamic stage–discharge rating for a migrating tropical river: Kalman-filtered rating shifts and their uncertainty at Bweramule, Semliki River, Uganda

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24367492.svg)](https://doi.org/10.5281/zenodo.24367492) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kaitale/semliki-sem9.10-dynamic-rating/blob/main/reproduce_paper.ipynb)

**Charles Kakooza (Ministry of Water and Environment); Monica Nafula (KAI Water & Energy); Geoffry Kasibante (Ministry of Water and Environment); Faizal Waiswa (Ministry of Water and Environment); Benon Zaake (Ministry of Water and Environment)**

Preprint (not peer reviewed), version 2.0, 19 November 2026, prepared to the author guidelines of *Journal of Hydrology: Regional Studies*.

## Abstract

**Study region:** The Semliki River at Bweramule (6,634 km²), on the Uganda–Democratic Republic of Congo border between Lake Edward and Lake Albert: the forecast point of an operational flood forecasting system, on a meander bend that migrates up to 11 m a year.

**Study focus:** A fixed stage–discharge rating is exact on the day of its last gauging and silently wrong afterwards if the control moves. We replace it with a dynamic rating – the installed curve plus a stage shift that follows a random walk – estimated from 176 current-meter gaugings (1950–1975) and 11 ADCP gaugings (2024–2025) with a Kalman filter and smoother, the drift rate and gauging error fitted by maximum likelihood. Its 90 % band widens every day without a gauging. We test it against a fixed rating by one-step-ahead prediction, and measure the loop from unsteady flow with the Jones formula against an unsteady HEC-RAS hindcast.

**New hydrological insights for the region:** The control at Bweramule drifts at 6.8 cm per square-root year. The dynamic rating raises the log-likelihood of the 1950–1975 gaugings from 153 to 253 (AIC −304 to −503), halves the gauging scatter (9.5 to 4.2 %) and predicts each next gauging with a 90 % band 42 % narrower (±0.10 against ±0.17 m) that is equally reliable. The control stood about 0.17 m low in 1950–1956 and 0.17 m high by 1975. In October 2026, 375 days after the last gauging, the stage for a given flow is known to ±0.13 m (about ±6 % in flow); a gauging every six months would hold it at ±0.10 m. Hysteresis is small: at the 99th-percentile rate of stage change the Jones correction is 0.6 %, and HEC-RAS agrees.

## Reproduce the paper

One notebook, top to bottom:

```bash
pip install -r requirements.txt
jupyter lab reproduce_paper.ipynb      # then Run All
```

or open it in [Google Colab](https://colab.research.google.com/github/Kaitale/semliki-sem9.10-dynamic-rating/blob/main/reproduce_paper.ipynb); the notebook downloads anything it needs from this repository.

## Contents

| File | What it is |
|---|---|
| `Kakooza_et_al_2026_Bweramule_dynamic_rating_preprint.pdf` | the paper |
| `reproduce_paper.ipynb` | the one notebook: reproduces Tables 4–5, Fig. 5 and the Jones corrections, and runs the full estimation on a synthetic record or on your own gaugings (seconds) |
| `semliki_dynrating.py` | Kalman filter, RTS smoother, maximum-likelihood drift rate, rating with band, Jones correction |
| `rating_R2026-10_table.csv` | the installed rating R2026-10 at Bweramule (stage on the ADCP staff gauge, m; discharge, m³/s) |
| `requirements.txt` | numpy, matplotlib |

## Data not included

The Bweramule gaugings, telemetry and HEC-RAS model are held by the Ministry of Water and Environment, Uganda, and are available on request; with them, the last cell of the notebook reproduces Table 3 and Figs. 2–3.

## How to cite

Kakooza, C., Nafula, M., Kasibante, G., Waiswa, F., Zaake, B., 2026. A dynamic stage–discharge rating for a migrating tropical river: Kalman-filtered rating shifts and their uncertainty at Bweramule, Semliki River, Uganda. Preprint, version 2.0. Zenodo. https://doi.org/10.5281/zenodo.24367492

## Licence

Manuscript and data: CC BY 4.0. Code: MIT. See `LICENSE`.

## Funding

Carried out under the Semliki Flood Forecasting System of the Ministry of Water and Environment, Uganda, with support from the World Food Programme.

## Version history

* **2.0** (19 Nov 2026): the work as a journal-style preprint with one reproduction notebook.
* 1.x (Oct 2026): the original technical report and analysis notebooks (superseded).
