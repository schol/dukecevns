# dukecevns
Simple code for sharing within Duke group for CEvNS rate checks

# Compilation
For sns_rates.cc:  need json installed,  https://github.com/nlohmann/json.git. Simply run

```sh
cd /path/to/dukecevns
git clone https://github.com/nlohmann/json.git
make sns_rates
```

# Usage

```sh
cd /path/to/dukecevns
./sns_rates file-in-folder-named-jsonfiles-without-path-without-dot-json
```

Example `json` files can be found in ORNL GitLab repository `COHERENTProposal2018/assumptions/jsonfiles`.

The output files are located in a new folder named `out`.
