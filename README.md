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
- everything in `qf/` and `jsonfiles/`

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
