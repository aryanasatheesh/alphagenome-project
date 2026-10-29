#!/usr/bin/env python3
"""
build_track_classes.py
----------------------
Label every AlphaGenome human output channel with what biological sample it
represents: brain tissue, cells from brain, stem-cell-derived neural cells,
tumour cell lines, peripheral nerve, eye, pluripotent lines, or other.

Why
    Comparing "brain" with "non-brain" tracks needs a defensible definition of
    brain. Matching words in sample names is unsafe ("gastrocnemius" contains
    "astro", "Peyer's patch" contains "eye"), so labels come from:
      1. the ontology term each track carries (UBERON anatomy / Cell Ontology),
         used to FIND every candidate neural or nervous-system sample, and
      2. a reviewed table, phase2/metadata/biosample_curation.tsv, that makes the
         final call for each candidate and records the evidence.
    The script stops if the ontology finds a candidate that has not been reviewed,
    so a new metadata release cannot silently change the brain set.

How a track is described (all from the AlphaGenome metadata)
    ontology term   UBERON (tissue), CL (cell type), EFO/CLO (cell line),
                    NTR (ENCODE provisional term)
    biosample type  TISSUE, PRIMARY_CELL, CELL_LINE, IN_VITRO_DIFFERENTIATED_CELLS,
                    ORGANOID
    life stage      embryonic, newborn, child, adult, unknown (set per track)

Columns added per track
    group           CNS_tissue | CNS_primary_cell | CNS_organoid | neural_in_vitro |
                    neural_crest_in_vitro | CNS_non_neural | CNS_tumor_line | PNS |
                    PNS_tumor_line | eye | pluripotent | non_neural | padding
    analysis_set    brain        = CNS_tissue + CNS_primary_cell (cells from human CNS)
                    neural_model = neural_in_vitro + CNS_organoid (made from stem cells)
                    background   = non_neural
                    excluded     = everything else (neither clearly brain nor clearly unrelated)
                    padding      = padding channels
    origin          donor_tissue, donor_cells, ESC_derived, iPSC_derived,
                    stem_cell_derived, tumor_line, pluripotent_line, cell_line
                    For lab-made cells, ESC vs iPSC is INFERRED from the life stage
                    recorded for the parent line's donor: embryonic -> ESC,
                    child/adult -> iPSC.
    life_stage_group  prenatal | perinatal | postnatal | unknown | n/a (cultured lines)
    cell_class, brain_region, evidence   from the curation table

Usage (from the repo root)
    python phase2/build_track_classes.py
Defaults: phase2/metadata/human_output_metadata.textproto,
          phase2/metadata/biosample_curation.tsv,
          phase2/metadata/ontologies/{uberon-basic.obo,cl-basic.obo},
          output phase2/metadata/track_classes.tsv (+ .summary.tsv, .log.txt)

Ontology files (download once; gitignored as large):
    curl -L -o phase2/metadata/ontologies/uberon-basic.obo \\
      https://github.com/obophenotype/uberon/releases/latest/download/uberon-basic.obo
    curl -L -o phase2/metadata/ontologies/cl-basic.obo \\
      https://github.com/obophenotype/cell-ontology/releases/latest/download/cl-basic.obo
"""

import argparse
import csv
import datetime
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
META = HERE / "metadata"

# Ontology roots used to find candidates
CNS, PNS, EYE = "UBERON:0001017", "UBERON:0000010", "UBERON:0000970"
NEURAL_ROOTS = {"CL:0002319",  # neural cell
                "CL:0000540",  # neuron
                "CL:0000125",  # glial cell
                "CL:0011020",  # neural progenitor cell
                "CL:0000047",  # neural stem cell
                "CL:0000133",  # neurectodermal cell
                "CL:0011012"}  # neural crest cell

ANALYSIS_SET = {"CNS_tissue": "brain", "CNS_primary_cell": "brain",
                "neural_in_vitro": "neural_model", "CNS_organoid": "neural_model",
                "non_neural": "background", "padding": "padding"}
ONTO_PREFIX = {"ONTOLOGY_TYPE_CL": "CL", "ONTOLOGY_TYPE_UBERON": "UBERON",
               "ONTOLOGY_TYPE_EFO": "EFO", "ONTOLOGY_TYPE_CLO": "CLO", "ONTOLOGY_TYPE_NTR": "NTR"}
SUMMARY_HEADS = ["ATAC", "DNASE", "CHIP_HISTONE", "RNA_SEQ", "CAGE", "CHIP_TF"]


# ----------------------------------------------------------------------------- parsing
def parse_textproto(path):
    """Tiny parser for the metadata textproto: nested blocks become lists of dicts."""
    root, stack = {}, []
    stack.append(root)
    for raw in open(path, encoding="utf-8"):
        line = raw.strip()
        if not line:
            continue
        if line.endswith("{"):
            node = {}
            stack[-1].setdefault(line[:-1].strip(), []).append(node)
            stack.append(node)
        elif line == "}":
            stack.pop()
        else:
            key, val = line.split(":", 1)
            val = val.strip()
            if val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            stack[-1].setdefault(key.strip(), []).append(val)
    return root


def first(d, key):
    v = d.get(key)
    return v[0] if v else ""


def metadata_tracks(path):
    """One dict per output channel, in channel order, for every output head."""
    rows = []
    for om in parse_textproto(path)["output_metadata"]:
        head = first(om, "output_type").replace("OUTPUT_TYPE_", "")
        idx = 0
        for block in om.get("tracks", []) + om.get("junctions", []):
            for m in block.get("metadata", []):
                onto = (m.get("ontology_term") or [{}])[0]
                bio = (m.get("biosample") or [{}])[0]
                prefix = ONTO_PREFIX.get(first(onto, "ontology_type"), "")
                oid = first(onto, "id")
                if prefix in ("CL", "UBERON", "EFO", "CLO") and oid:
                    curie = f"{prefix}:{int(oid):07d}"
                elif prefix and oid:
                    curie = f"{prefix}:{oid}"
                else:
                    curie = ""
                rows.append(dict(
                    output_type=head, track_index=idx, track_name=first(m, "name"),
                    is_padding=first(m, "name") == "Padding",
                    ontology_curie=curie,
                    biosample_type=first(bio, "type").replace("BIOSAMPLE_TYPE_", ""),
                    biosample_name=first(bio, "name"), stage=first(bio, "stage"),
                    assay=first(m, "assay"), data_source=first(m, "data_source"),
                    gtex_tissue=first(m, "gtex_tissue"), histone_mark=first(m, "histone_mark_code"),
                    tf=first(m, "transcription_factor_code"), strand=first(m, "strand").replace("STRAND_", "")))
                idx += 1
    return rows


def load_obo(*paths):
    """Term names and is_a / part_of parents. Qualifiers like {is_inferred="true"} are dropped."""
    names, parents = {}, defaultdict(set)
    for p in paths:
        term = None
        for line in open(p, encoding="utf-8"):
            line = line.rstrip("\n")
            if line.startswith("["):
                term = None if line != "[Term]" else ""
                continue
            if term is None or not line:
                continue
            if line.startswith("id: "):
                term = line[4:].strip()
            elif not term:
                continue
            elif line.startswith("name: "):
                names.setdefault(term, line[6:].strip())
            elif line.startswith("is_a: "):
                parents[term].add(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                parents[term].add(line.split()[2])
    return names, parents


def ancestors(term, parents):
    seen, stack = set(), [term]
    while stack:
        for p in parents.get(stack.pop(), ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return seen


# ----------------------------------------------------------------------------- labelling
def life_stage_group(stage, btype):
    if btype in ("CELL_LINE", "IN_VITRO_DIFFERENTIATED_CELLS", "ORGANOID"):
        return "n/a"
    parts = {s for s in stage.split(",") if s and s != "unknown"}
    if not parts:
        return "unknown"
    if parts == {"embryonic"}:
        return "prenatal"
    if "newborn" in parts:
        return "perinatal"
    if parts <= {"child", "adult"}:
        return "postnatal"
    return "mixed"


def origin(row, group):
    t = row["biosample_type"]
    if t == "TISSUE":
        return "donor_tissue"
    if t == "PRIMARY_CELL":
        return "donor_cells"
    if t == "ORGANOID":
        return "stem_cell_derived"
    if t == "IN_VITRO_DIFFERENTIATED_CELLS":
        parts = set(row["stage"].split(","))
        if parts == {"embryonic"}:
            return "ESC_derived"
        if parts and parts <= {"child", "adult"}:
            return "iPSC_derived"
        return "stem_cell_derived"
    if group in ("CNS_tumor_line", "PNS_tumor_line") or (group == "eye" and t == "CELL_LINE"):
        return "tumor_line"
    if group == "pluripotent":
        return "pluripotent_line"
    return "cell_line"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description="Label AlphaGenome output tracks by biological sample.")
    ap.add_argument("--metadata", default=str(META / "human_output_metadata.textproto"))
    ap.add_argument("--curation", default=str(META / "biosample_curation.tsv"))
    ap.add_argument("--uberon", default=str(META / "ontologies" / "uberon-basic.obo"))
    ap.add_argument("--cl", default=str(META / "ontologies" / "cl-basic.obo"))
    ap.add_argument("--out", default=str(META / "track_classes.tsv"))
    args = ap.parse_args()
    for p in (args.metadata, args.curation, args.uberon, args.cl):
        if not Path(p).exists():
            sys.exit(f"Missing input: {p}\n(see the docstring for how to download ontology files)")

    tracks = metadata_tracks(args.metadata)
    names, parents = load_obo(args.uberon, args.cl)
    curation = list(csv.DictReader(open(args.curation, newline="", encoding="utf-8"), delimiter="\t"))

    # ---- match curated rows to biosamples; most specific key wins
    biosamples = {(t["ontology_curie"], t["biosample_type"], t["biosample_name"])
                  for t in tracks if not t["is_padding"]}
    by_key = defaultdict(list)
    for i, c in enumerate(curation):
        by_key[(c["ontology_curie"], c["biosample_type"], c["biosample_name"])].append(i)
    match, used, problems = {}, set(), []
    for b in sorted(biosamples):
        curie, btype, name = b
        hits = []
        for key in ((curie, btype, name), (curie, "*", name), ("*", btype, name), ("*", "*", name)):
            hits = by_key.get(key, [])
            if hits:
                break
        if len(hits) > 1:
            problems.append(f"more than one curated row matches {b}")
        if hits:
            match[b] = curation[hits[0]]
            used.add(hits[0])
    for i, c in enumerate(curation):
        if i not in used:
            problems.append(f"curated row matches nothing in the metadata: "
                            f"{c['ontology_curie']} / {c['biosample_type']} / {c['biosample_name']}")

    # ---- every ontology candidate must have been reviewed
    unreviewed = []
    for b in sorted(biosamples):
        curie = b[0]
        if not curie or curie.split(":")[0] not in ("UBERON", "CL"):
            continue
        anc = ancestors(curie, parents) | {curie}
        flags = [f for f, hit in (("CNS", CNS in anc), ("PNS", PNS in anc), ("eye", EYE in anc),
                                  ("neural", bool(NEURAL_ROOTS & anc))) if hit]
        if flags and b not in match:
            unreviewed.append(f"{b} flagged {flags}")
    if unreviewed:
        problems.append("ontology candidates not in the curation table:\n    " + "\n    ".join(unreviewed))
    if problems:
        sys.exit("STOPPED - fix the curation table first:\n  " + "\n  ".join(problems))

    # ---- label every channel
    out_rows = []
    for t in tracks:
        if t["is_padding"]:
            c = {"group": "padding", "cell_class": "", "brain_region": "", "evidence": "padding channel"}
        else:
            c = match.get((t["ontology_curie"], t["biosample_type"], t["biosample_name"]),
                          {"group": "non_neural", "cell_class": "", "brain_region": "", "evidence": "not in curation table"})
        g = c["group"]
        out_rows.append({**{k: t[k] for k in ("output_type", "track_index", "track_name", "is_padding",
                                               "ontology_curie", "biosample_type", "biosample_name", "stage")},
                         "life_stage_group": "" if t["is_padding"] else life_stage_group(t["stage"], t["biosample_type"]),
                         "group": g, "analysis_set": ANALYSIS_SET.get(g, "excluded"),
                         "origin": "" if t["is_padding"] else origin(t, g),
                         "cell_class": c["cell_class"], "brain_region": c["brain_region"],
                         **{k: t[k] for k in ("assay", "histone_mark", "tf", "data_source", "gtex_tissue", "strand")},
                         "evidence": c["evidence"]})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

    # ---- summary
    groups = ["CNS_tissue", "CNS_primary_cell", "CNS_organoid", "neural_in_vitro", "neural_crest_in_vitro",
              "CNS_non_neural", "CNS_tumor_line", "PNS", "PNS_tumor_line", "eye", "pluripotent", "non_neural", "padding"]
    cnt = Counter((r["output_type"], r["group"]) for r in out_rows)
    summary_path = Path(str(out).replace(".tsv", ".summary.tsv"))
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("group\tanalysis_set\t" + "\t".join(SUMMARY_HEADS) + "\n")
        for g in groups:
            f.write(f"{g}\t{ANALYSIS_SET.get(g, 'excluded')}\t" + "\t".join(str(cnt[(h, g)]) for h in SUMMARY_HEADS) + "\n")
    print(f"{'group':<22}{'set':<13}" + "".join(f"{h:>13}" for h in SUMMARY_HEADS))
    for g in groups:
        print(f"{g:<22}{ANALYSIS_SET.get(g, 'excluded'):<13}" + "".join(f"{cnt[(h, g)]:>13}" for h in SUMMARY_HEADS))
    brain_stage = Counter((r["output_type"], r["life_stage_group"]) for r in out_rows if r["analysis_set"] == "brain")
    print("\nbrain set by life stage:")
    for h in SUMMARY_HEADS:
        print(f"  {h:<13}" + "  ".join(f"{s} {brain_stage[(h, s)]}" for s in ("prenatal", "perinatal", "postnatal", "unknown", "mixed") if brain_stage[(h, s)]))

    with open(str(out).replace(".tsv", ".log.txt"), "w", encoding="utf-8") as f:
        f.write(f"date            : {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"metadata        : {Path(args.metadata).name}  md5 {md5(args.metadata)}\n")
        f.write(f"curation        : {Path(args.curation).name}  md5 {md5(args.curation)}  ({len(curation)} rows)\n")
        for p in (args.uberon, args.cl):
            ver = next((l.split(':', 1)[1].strip() for l in open(p, encoding='utf-8') if l.startswith('data-version')), '?')
            f.write(f"ontology        : {Path(p).name}  {ver}\n")
        f.write(f"channels        : {len(out_rows)}  (padding {sum(r['is_padding'] for r in out_rows)})\n")
    print(f"\nOUTPUT : {out}\nSUMMARY: {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())