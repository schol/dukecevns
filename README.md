Simple code to calculate CEvNS rates in various experiments.

# Pre-requisites

## JSON for Modern C++

[JSON for Modern C++](https://github.com/nlohmann/json) is needed to read experiment configuration files in [JSON][] format. It has been added as a [git][] [submodule][] of [dukecevns][], and can be automatically fetched with the `--recurse-submodules` argument for `git clone`:

```sh
git clone --recurse-submodules https://github.com/schol/dukecevns
```

For an existing [dukecevns][] repository, the submodule can be fetched with the following commands:

```sh
cd /path/to/dukecevns
rm -fr json
git pull
git submodule update --init
```

## CERN ROOT

[ROOT][] is needed to compile [dukecevns][]. The [Makefile](Makefile) included in [dukecevns][] can be used to compile [dukecevns][] on Linux and MacOS with [ROOT][] installed. For Windows users, a [Docker][] image [physino/root][] that has the latest [ROOT][] installed in [Fedora][] Linux can be used to compile [dukecevns][]. If [Docker Desktop][] has been installed and is running, open a terminal in the folder of [dukecevns][] by highlighting the address bar of the file explorer and typing `cmd` and return. In the terminal, type

```sh
docker-compose run --rm sh
```

Wait and let [Docker][] do its magic until you see a familiar [SHELL][] [prompt][] showing up:

```sh
root@fedora:~/dukecevns $
```

Now you are in a [Docker][] container running [Fedora][] Linux with [ROOT][] pre-installed in it, where you can simply type `make` to compile [dukecevns][] as in a real Linux machine.

## Optional dependencies

Some calculations rely on additional code hosted in private git repositories, such as
- <https://code.ornl.gov/COHERENT/COHERENTProposal2018>

One can run the included `bootstrap.sh` to grab the required packages automatically:

```sh
./bootstrap.sh
```

# Compilation

```sh
cd /path/to/dukecevns
# create libDiffSpec_1.0.a
make
# create an executable 'sns_rates' (needs private code)
# done automatically if bootstrap.sh is used
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

If `efftype = "qc"`, `qcgsname` or `qcsmearing` should be used. If `efftype="eee"`, `gsname` should be used. They cannot be mixed.

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

[JSON]: https://www.json.org/json-en.html
[git]: https://git-scm.com/book/en/v2/Getting-Started-What-is-Git%3F
[submodule]: https://git-scm.com/book/en/v2/Git-Tools-Submodules
[dukecevns]: https://github.com/schol/dukecevns
[ROOT]: https://root.cern.ch
[Docker]: https://www.docker.com
[physino/root]: https://hub.docker.com/r/physino/root
[Fedora]: https://en.wikipedia.org/wiki/Fedora_Linux
[Docker Desktop]: https://www.docker.com/products/docker-desktop
[SHELL]: https://en.wikipedia.org/wiki/Unix_shell
[prompt]: https://en.wikibooks.org/wiki/Guide_to_Unix/Explanations/Shell_Prompt
