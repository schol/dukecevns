# dukecevns
Simple code for sharing within Duke group for CEvNS rate checks

# Pre-requisition

- https://github.com/nlohmann/json.git
- https://code.ornl.gov/COHERENT/COHERENTProposal2018/tree/master/assumptions (private)

Install json:
```sh
cd /path/to/dukecevns
git clone https://github.com/nlohmann/json.git
```

Grab the following files from ORNL GitLab and put them in `/path/to/dukecevns`:

- `get_flavor_weight.cc`
- `sns_out_BERT_convolved.root` and `sns_out_BERT.root`
- everything in `eff/`, `gs/`, `qf/` and `jsonfiles/`

# Compilation

```sh
cd /path/to/dukecevns
make sns_rates
```

# Usage

```sh
cd /path/to/dukecevns
./sns_rates base-name-of-a-file-in-the-folder-named-jsonfiles
```

Example `json` files can be found in ORNL GitLab repository `COHERENTProposal2018/assumptions/jsonfiles`.

The output files are located in a new folder named `out`.

## Experimental setup

```json
"timewindow": {
  "start": 1400.0,
  "end": 7400.0
},
```

Kate recommends using the convolved flux histogram. The convolved flux histogram is shifted in time; nominal window to use is 1400 to 7400, although that may get optimized if one knows something about the background.

## Detector response

```json
"detectorresponse": {
  "efftype": "qc", <- qc: charge (q) collected, it can also be # of PEs, ADCs, etc.
  "efficiencyfile": "none", <- file in eff/ folder
  "stepthresh": 2.0, <- lower energy threshold
  "upperthresh": 6.0, <- upper limit of energy
  "qftype": "poly",
  "qfname": "nai",
  "qcperkeVee": 30.0, <- ionization (or light) yield
  "qcbinning": 1, <- bin width
  "gsname": "none", <- electron-equivalent energy resolution (file in gs/)
  "qcsmearing": "poisson", <- "qcgsname", "gsname" are not needed if this is set
  "qcgsname": "none" <- energy resolution in terms of Qc
},
```

### gs file format

`polysqrt` means that it's a polynomial in terms of the square root of energy. The second line defines the range, and the third line includes the coefficients of the formula defined here: https://coherent.phy.duke.edu/wiki/Assumptions#NaI_Smearing (private).

### efficiency file format

Column one: energy, # of PEs, or Qc, etc.
Column two: efficiency

### qf file format

Row one: energy range
Row two: quenching factor

If `qftype==poly`, the second row contains polynomial coefficients: C0, C1, C2, etc. so that QF = C0 + C1*E + C2*E^2.

# Output

## sns_diff_rates_quenched-alliso-*.out
Column one: Energy [MeVee]
Column two: Number of events per MeVee
