# ZFN Design Workbench · User Manual

Version 0.7.0 · 2026-11-18

## 1 What this is

A design tool for zinc-finger nucleases (ZFNs). You give it a DNA sequence and it does four things:

- **S1 DesignScore** — scans candidate sites, scores each design on six weighted components, ranks them, and produces an all-target report.
- **S2 Target filter** — optional. Puts each candidate back on the whole genome; candidates that land in repeats, or occur 3 or more times genome-wide, are removed entirely.
- **S3 Off-target Score** — optional. Runs a whole-genome off-target scan on the top N designs and reports a burden.
- **S4 Fusion design** — attaches effector domains (VP64 / KRAB / FokI) to selected designs, emits the fusion protein and a humanised CDS, plus an order form.

All computation runs locally or on the lab server. Sequences are never sent to a third party.

## 2 Access

### 2.1 Web UI

Open in a browser:

- Local: `http://127.0.0.1:8797/`
- Public example: `https://zfn.mobilome.cn/` (our deployment)

The `EN / 中文` button in the top-left corner switches the whole page language and remembers the choice.

### 2.2 Run it yourself

```bash
python -m zfn_design gui --port 8797 --no-browser
```

`--out-root` sets the result directory; the default is `<project>/out/gui/`.

## 3 Screen map

Inputs on the left, results on the right.

- Input column, top to bottom: DNA box, architecture, spacer length, finger array, parameters, Target filter (collapsed), two buttons.
- The two buttons — **Design Score (S1 + S2)** and **Off-target Score (S3)** — each have their own progress bar.
- Result column: statistic cards, design list, off-target table, S4 panel, file links.

## 4 Five steps

### 4.1 Paste a sequence

Paste DNA into the large box. Plain text or FASTA both work; a FASTA header becomes the sequence name. "Load example" fills in a sequence known to run.

### 4.2 Optional: pick a genome for Target filter

Expand "Target filter · whole-genome repeat filtering". Choose a genome FASTA (e.g. `GRCh38.p14.fna`); the second dropdown optionally takes a RepeatMasker annotation (`.out` or `.bed`).

With a genome selected, filtering is genome-wide:

- any landing that overlaps a RepeatMasker interval removes the whole candidate;
- a candidate occurring 3 or more times genome-wide (reverse complement included) is removed;
- removed sites are **not scored**.

Without a genome, only DesignScore runs and repeats are ignored.

### 4.3 Press Design Score

The result column shows statistics, the candidate table, and design cards sorted by Design_score. A "Download bundle (zip)" link appears when done.

Without a genome this takes seconds. With GRCh38 it takes about 3 minutes (3.3 billion bases scanned).

### 4.4 Optional: press Off-target Score

Precondition: **only runs that went through Target filter may run S3**. A run without filtering is rejected with an explanation.

Enter Top N (default 10) and the per-side mismatch cap (default 2), then press the button. It reads `report.json` from the previous run and scans the genome for the top N designs.

Two numbers come out: `off-target score` (0–100, higher is better) and `burden` (lower is better). They are shown **side by side and are never added to Design_score**.

A genome scan per design takes about 1–2 minutes.

### 4.5 Optional: press S4 to build fusion protein and CDS

In the "④ S4 · Fusion design" panel:

- **Top N** — how many designs to use (default 3);
- tick effector architectures: VP64 (activator), KRAB (repressor), FokI (nuclease); multiple allowed;
- press "Generate fusion protein and CDS".

Outputs: fusion protein sequences, humanised CDS, an order form, and a rebuilt zip. S4 only reads this run's `report.json`; it does not recompute S1–S3.

## 5 Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| Architecture | whatever the finger library supports | Monomer (activator/repressor) or paired nuclease |
| Spacer length | 5–7 | Only meaningful for paired designs; ignored for monomers |
| Finger array | three fingers | Three fingers read 9 bp, six read 18 bp; other counts use the generic scaffold |
| Designs per site | 3 | Designs generated per candidate site |
| Max sites scored | 20 | Candidates kept before scoring |
| Max mismatches per side | 2 | Off-target cutoff (0–3) |
| Top N (S3) | 10 | Designs entering the off-target scan |
| Mismatch cap (S3) | 2 | Mismatch cap during the scan |
| Top N (S4) | 3 | Designs entering fusion design |

## 6 Output files

One directory per run, by default `<project>/out/gui/<timestamp>_<name>_<short-hash>/`.

| File | Content |
| --- | --- |
| `report.html` | Full report, open in a browser |
| `report.json` | Structured result; downstream scripts read this |
| `sites.tsv` | Candidate sites |
| `designs.tsv` | Designs with score components |
| `offtarget_genome.tsv` | S3 results (only if S3 ran) |
| `fusions.json` | S4 entries (only if S4 ran) |
| `proteins.tsv` | Fusion proteins, one chain per row |
| `cds.tsv` | CDS details: start/stop codon, GC, sites |
| `cds.fasta` | CDS sequences, ready for ordering |
| `order_form.csv` | Synthesis order form |
| `zfn_bundle.zip` | Everything above, packed |

## 7 How to read the numbers

**Design_score** (0–100) is a weighted sum of six components: finger evidence strength, triplet quality, GC balance, finger diversity, sequence complexity, self-complementarity risk. A ≥ 85, B ≥ 70, C ≥ 55.

The score means "this site looks plausible under published design rules". It is **not an activity prediction**.

**burden** is the off-target load: every near-match site weighted by its mismatch count (0 → 1.0, 1 → 0.25, 2 → 0.06, 3 → 0.015, 4 → 0.004), summed. `off-target score = 100 − 25 × burden`, floored at 0. When the hit cap (500) is reached the burden is a lower bound and the UI says so.

**CDS status**:

- `CDS ok` — starts with `ATG`, ends with a stop codon, passes GC and site checks;
- `needs review` — sequence produced, but red flags remain (e.g. off-target not computed, unverified part);
- `part missing` — no reliable sequence for some part, no CDS emitted.

Read every red flag. Do not skip them.

## 8 Command line

```bash
# S1 only (optional S2: add --genome / --repeats)
python -m zfn_design run --seq "ACGT..." --out out/run1

# with whole-genome filtering
python -m zfn_design run --seq input.fa --genome refs/GRCh38.p14.fna \
    --repeats refs/GRCh38.p14.rm.out --out out/run2

# S4: fusion design on an existing result directory
python -m zfn_design fusion out/run2 --top 3 --arch single_activator paired_nuclease
```

Subcommands: `run`, `fusion`, `scan`, `library`, `info`, `gui`.

## 9 Reports and downloads

Every run on the page links its files; "Download bundle (zip)" grabs all of them at once.

`report.html` is the full numbered report: input, Target filter, finger library, candidate sites, triplets, designs, Off-target score, fusion, weights.

## 10 Troubleshooting

**Nothing happens after pressing the button.** Watch the progress bar. If it never moves, the request never left the browser; check the console for `POST /api/run`.

**"Pick a genome first".** S3 and Target filter need a genome FASTA. Select one and rerun Design Score.

**"Only candidates that passed Target filter can run Off-target Score".** The previous run skipped filtering. Enable the genome and rerun.

**It is slow.** S1+S2 on GRCh38 takes about 3 minutes; each S3 design takes 1–2 minutes. Everything is single-process; do not press twice.

**Zero designs.** All candidates were filtered out, or the finger library lacks triplets. Check the "missing triplets" count in the report.

**"Part missing" in the CDS.** At least one part in that fusion has no reliable sequence. Supply it in the S4 parts JSON, or pick another architecture.

**Can I order it?** The sequence, yes. Activity, no. Every design must be validated by T7E1 or amplicon deep sequencing before it leaves the bench.

## 11 Limits and caveats

- Scoring follows rules distilled from public literature. No activity prediction, no structure modelling.
- Genome filtering only knows RepeatMasker annotation and copy counts; unannotated low-complexity regions slip through.
- Off-target search is a half-site similarity scan, not a cleavage probability.
- Finger scaffolds derive from a 1AAY-derived template; inter-module linker geometry has no structural evidence.
- Codon optimisation avoids common restriction sites and homopolymers, but does not check splice sites.

## 12 Version and maintenance

Current version 0.7.0. The service directory is chosen by the deployer; restart with `zfn_gui.sh restart`; logs in `gui.log`.

Send issues and suggestions to Song Ji or Song Chengyi.
