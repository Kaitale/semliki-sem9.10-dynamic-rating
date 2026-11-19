# The Dynamic Rating at Bweramule: a Stage–Discharge Relation that Moves with the River

**SFFS Technical Report SEM.9.10, version 1.0 (4 October 2026)** — Semliki Flood Forecasting System, Ministry of Water and Environment, Uganda.

The report replaces the fixed rating at Bweramule (Semliki River, DWRM 85205) with a dynamic rating: the installed curve plus a stage shift that drifts as a random walk, estimated from all 191 gaugings (1950–2025) by a Kalman filter, with a 90 % band that widens without gaugings and narrows at each new one, and the Jones correction for unsteady flow.

**Key findings:** a drifting rating halves the gauging scatter of 1950–75 (9.5 % to 4.2 %) and raises the log-likelihood from 153 to 253; the drift rate is 6.8 cm per √year; at the issue the stage for a given flow was known to ±0.13 m (90 %), and an ADCP gauging about every six months would keep it within ±0.10 m; the Jones correction is under 1 % at Bweramule.

## Files

| File | What it is |
|---|---|
| `SEM.9.10.R.1.0.SFFS.Dynamic.Rating.Bweramule.pdf` | The technical report |
| `semliki_dynrating.py` | Kalman filter and RTS smoother for a random-walk stage shift, maximum-likelihood drift rate, rating with a 90 % band, band growth, Jones correction |
| `01_dynamic_rating_kalman.ipynb` | Fit, smoothed shift, fixed-against-drifting comparison, modern shift, band growth |
| `02_jones_loop_correction.ipynb` | The Jones loop correction with the report's slope and width |

Run with Python 3.10+ and `pip install -r requirements.txt`.

**Data not included.** The 191 gaugings at Bweramule (DWRM 85205) are held by the Ministry of Water and Environment, Uganda, and are available from the Ministry on request. The notebooks run on a synthetic example (clearly labelled) with a drift rate and gauging error close to those of the report; set `GAUGINGS` to a CSV with columns `date, stage_m, Q_m3s, era` to use real data.

## Authors

Charles Kakooza, Monica Nafula, Geoffry Kasibante, Faizal Waiswa and Benon Zaake.

Prepared under the Semliki Flood Forecasting System (SFFS) of the Ministry of Water and Environment (MWE), Uganda, with support from the World Food Programme (WFP).

## Licence

Report and data: CC BY 4.0. Code and notebooks: MIT. See `LICENSE`.

## How to cite

See `CITATION.cff`, or the DOI of the Zenodo record of this repository.
