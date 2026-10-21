# Testing the liftover package

User is checking two things:

1. **Does it work** on any machine, following only the written
   instructions?
2. **Is it understandable**: could you use it on your own SNP list without
   asking author (Aryana) for guidance?

Please test version **`liftover-v1.0`**.
It takes about 30 minutes (plus 5 optional minutes on Hoffman2 for step 6).

**What the package does:** moves SNP positions from the hg19 genome build to
hg38, refusing any SNP whose position can't be moved unambiguously. Background
is in [`docs/liftover_README.md`](liftover_README.md).

Copy each command exactly. Where something says **Expected**, compare what you
see. If anything differs, note it.

---

## 1. Get the code

You need a GitHub account with access to the repo.

On **Hoffman2**, work in your own folder:

```bash
cd $HOME          # or your own project/scratch folder
git clone https://github.com/aryanasatheesh/alphagenome-project.git
cd alphagenome-project
git checkout liftover-v1.0
```

On a **Mac**, same commands in Terminal from any folder you like.

When git asks for a password, it needs a GitHub personal access token, not your
GitHub password (GitHub → Settings → Developer settings → Personal access
tokens).

**Expected:** `git checkout` says you are in "detached HEAD" state. That is
normal for a tagged version. Please don't commit or push anything; keep notes
for your feedback instead.

## 2. Set up Python

The package needs Python 3.6 or newer and one library, `pyliftover`.

**Mac:**
```bash
python3 -m venv ~/liftover-test-env
source ~/liftover-test-env/bin/activate
pip install pyliftover==0.4.1
```

**Hoffman2:**
```bash
python --version
python -c "import pyliftover; print(pyliftover.__version__)" || pip install --user pyliftover==0.4.1
```
If `python --version` says 2.x or "command not found", run
`module load anaconda3` and try again.

**Expected:** `import pyliftover` works and prints `0.4.1`.

## 3. Download the chain file

The chain file is the hg19→hg38 map. It isn't in the repo because it is a
downloaded data file.

```bash
curl -o phase1/data/hg19ToHg38.over.chain.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
ls -lh phase1/data/hg19ToHg38.over.chain.gz
```

**Expected:** a file of about 222K.

## 4. Run the automatic checks

```bash
python tests/liftover/check_liftover.py
```

This lifts six real PGC3 SNPs whose outcome is already known: two that lift
cleanly to known hg38 positions, and one each that cannot be lifted for a
different reason. It also triggers each of the script's safety checks. It
writes nothing inside the repo.

**Expected:** 22 lines saying `PASS`, ending with `ALL CHECKS PASSED`.
If any line says `FAIL`, copy the whole output into your feedback.

## 5. Reproduce a real result

Lift the 173 schizophrenia lead SNPs and compare with the file I committed:

```bash
python phase1/liftover_to_hg38.py \
    --input  phase1/data/scz_pgc3_lead_snps.tsv \
    --output /tmp/$USER-lead_hg38.tsv
cmp /tmp/$USER-lead_hg38.tsv phase1/data/scz_pgc3_lead_snps_hg38.tsv && echo "OUTPUT IDENTICAL"
cmp /tmp/$USER-lead_hg38.liftover_report.tsv phase1/data/scz_pgc3_lead_snps_hg38.liftover_report.tsv && echo "REPORT IDENTICAL"
```

**Expected:** the script prints `lifted 173`, then both `IDENTICAL` lines appear.
Also open `/tmp/$USER-lead_hg38.liftover_log.txt` and check it records your
command, the chain file and your pyliftover version.

## 6. (Optional, Hoffman2 only) Reproduce the whole-genome survey

This reruns the edge-case count over all 7.66 million PGC3 SNPs, reading the
GWAS file from my folder (read-only; nothing is written there). It needs a
compute node:

```bash
qrsh -l h_data=4G,h_rt=1:00:00
cd ~/alphagenome-project            # wherever you cloned it
python phase1/survey_liftover_edge_cases.py \
    --gwas /u/project/cluo/aryasath/alphagenome_project/data/PGC3_SCZ_european.vcf.tsv.gz \
    --out-prefix /tmp/$USER-pgc3
cmp /tmp/$USER-pgc3.summary.tsv results/liftover_survey/pgc3_european.summary.tsv && echo "SURVEY IDENTICAL"
exit
```

**Expected:** about 70 seconds, then `SURVEY IDENTICAL`. If you get
"Permission denied" on the GWAS file, tell me and skip this step.

## 7. Try it on your own SNPs

This is the most useful part. Make a small tab-separated file with a header, for
example in a spreadsheet saved as tab-delimited text:

```
rsid	chr	pos	my_note
rs...	6	...	anything
```

Use a few hg19 positions you know or can look up (Ensembl GRCh37 or the UCSC
browser set to hg19). Then:

```bash
python phase1/liftover_to_hg38.py --input my_snps.tsv --output my_snps_hg38.tsv
```

Things worth trying:
- write chromosomes as `chr6` in one file and `6` in another
- use different column names (e.g. `CHR`, `BP`) and work out from `--help` how to tell the script
- run it twice on the same file, or with a chromosome or position that doesn't exist
- check one or two lifted positions against Ensembl GRCh38 or UCSC hg38

When something goes wrong, does the message tell you what to do next?

## 8. Send feedback

Open an **Issue** on the repo (Issues tab → New issue), titled
`liftover-v1.0 test: <your name>`, and answer:

1. **Machine:** Mac or Hoffman2? Output of `python --version` and the pyliftover version.
2. **Steps 4–6:** did each give the expected result? Paste any `FAIL` lines or errors.
3. **Where did you get stuck or have to guess?** Even small things.
4. **Which error messages were unclear?** Quote them.
5. **Step 7:** what did you try, and did the output make sense?
6. **README:** was anything in `docs/liftover_README.md` confusing or missing?
7. **Time:** roughly how long did it take?

If you prefer, send the same answers by message. Thank you!