#!/usr/bin/env python3
"""
assign_tissue_categories.py
---------------------------
Give every real AlphaGenome track a TISSUE CATEGORY (organ system) and a SAMPLE KIND,
so variant effects can be summarised per category instead of "brain vs everything
else". (Choo, Nov 5: "find categorization of what common tissues are".)

Why categories
    "Background" mixes very different samples (DNase: 107 of 305 tracks are cell
    lines; ATAC: 92 of 167 are lymphoblastoid lines). A mean over all non-brain
    tracks is dominated by whatever is most common. Summarising per category,
    with each category weighted equally, avoids that.

How
    1. Neural samples keep the reviewed decision from track_classes.tsv:
         CNS_tissue, CNS_primary_cell, CNS_tumor_line -> CNS
         neural_in_vitro, CNS_organoid                -> neural_stem_cell_derived
         eye                                          -> eye_retina
         PNS, PNS_tumor_line, neural_crest_in_vitro   -> PNS_neural_crest
         pluripotent                                  -> pluripotent
         CNS_non_neural (brain vessels, meninges, pituitary) -> by ontology, step 2;
                                    if nothing matches -> CNS_non_neural
    2. Tissues and primary cells: the first category whose anchor term is an
       ancestor (is_a or part_of, UBERON + Cell Ontology) of the track's term.
       Cell-type anchors that occur in many organs (blood/immune, vessel, stromal)
       are checked before organs, so "fibroblast of lung" is stromal, not lung.
       Samples with no anchor get a reviewed entry in OVERRIDES (eye_anterior =
       cornea, lens, iris, ciliary body, trabecular meshwork: eye, but not retina).
    3. In-vitro-differentiated non-neural cells -> stem_cell_derived_other.
    4. Cell lines -> category of their tissue of origin, from CELL_LINE_ORIGIN below
       (Cellosaurus / ENCODE descriptions; 'low' confidence entries are marked).
    Anchor terms are looked up by exact name in the ontology files, and the script
    stops if any is missing, so a renamed term cannot silently drop a rule.

sample_kind and main_analysis
    sample_kind   tissue | primary_cell | in_vitro | organoid | cell_line
    main_analysis True for tissue, primary_cell, in_vitro, organoid.
                  Cell lines (cancer, immortalised, lymphoblastoid) are kept out of
                  the main comparison, the same rule that kept tumour lines out of
                  "brain"; they can be added back as a sensitivity analysis.

Output
    phase2/metadata/track_categories.tsv           one row per real track
    phase2/metadata/track_categories.summary.tsv   tracks per head x category (main / lines)

Usage (from the repo root)
    python phase2/assign_tissue_categories.py
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_track_classes import load_obo, ancestors  # noqa: E402  (same parser as the label build)

META = Path(__file__).resolve().parent / "metadata"
CLASSES = META / "track_classes.tsv"
OBO = [META / "ontologies" / "uberon-basic.obo", META / "ontologies" / "cl-basic.obo"]

NEURAL_GROUP = {"CNS_tissue": "CNS", "CNS_primary_cell": "CNS", "CNS_tumor_line": "CNS",
                "neural_in_vitro": "neural_stem_cell_derived", "CNS_organoid": "neural_stem_cell_derived",
                "eye": "eye_retina", "PNS": "PNS_neural_crest", "PNS_tumor_line": "PNS_neural_crest",
                "neural_crest_in_vitro": "PNS_neural_crest", "pluripotent": "pluripotent"}

# checked in this order; first match wins
RULES = [
    ("immune_blood", ["hematopoietic cell", "immune system", "hematopoietic system", "blood", "thymus",
                      "spleen", "bone marrow", "tonsil", "lymph node"]),
    ("vascular", ["endothelial cell", "vascular associated smooth muscle cell", "pericyte", "blood vessel",
                  "vasculature", "artery", "vein"]),
    ("stromal_connective_adipose", ["fibroblast", "mesenchymal stem cell", "stromal cell", "osteoblast",
                                    "chondrocyte", "smooth muscle cell", "adipose tissue", "bone tissue",
                                    "cartilage tissue"]),
    ("liver", ["liver", "hepatocyte"]),
    ("lung_airway", ["lung", "respiratory system", "trachea", "epithelial cell of lung"]),
    ("heart", ["heart", "cardiac muscle cell", "cardiac muscle tissue"]),
    ("skeletal_muscle", ["skeletal muscle tissue", "skeletal muscle myoblast", "myoblast", "musculature"]),
    ("kidney_urinary", ["kidney", "renal system", "urinary bladder", "ureter", "urothelium",
                        "kidney epithelial cell"]),
    ("digestive", ["digestive system", "stomach", "intestine", "colon", "esophagus", "pancreas",
                   "oral cavity", "tongue", "gingiva"]),
    ("endocrine", ["thyroid gland", "adrenal gland", "endocrine system", "pituitary gland"]),
    ("embryonic_extraembryonic", ["extraembryonic structure", "placenta", "trophoblast cell"]),
    ("reproductive_breast", ["reproductive system", "breast", "mammary gland", "prostate gland", "testis",
                             "ovary", "uterus", "vagina", "epithelial cell of prostate",
                             "mammary gland epithelial cell"]),
    ("skin", ["skin of body", "integument", "integumental system", "keratinocyte", "melanocyte"]),
    ("embryonic_extraembryonic", ["embryo", "embryonic structure"]),
]

# Cell lines: tissue of origin. (category, confidence, note)
_LCL = ("immune_blood", "high", "lymphoblastoid (EBV-B) line")
CELL_LINE_ORIGIN = {
    **{n: ("immune_blood", "high", "leukaemia / lymphoma / myeloma / immortalised blood line") for n in
       ["BLaER1", "CMK", "DND-41", "DOHH2", "HAP-1", "HL-60", "Jurkat", "Jurkat, Clone E6-1", "K562", "KBM-7",
        "KMS-11", "KOPT-K1", "Karpas-422", "Loucy", "MM.1S", "NAMALWA", "NB4", "NCI-H929", "OCI-LY1", "OCI-LY3",
        "OCI-LY7", "RPMI8226", "Raji", "SU-DHL-6", "THP-1"]},
    **{n: ("liver", "high", "hepatocellular carcinoma / hepatoblastoma") for n in ["HepG2", "HuH-7", "HuH-7.5"]},
    **{n: ("lung_airway", "high", "lung carcinoma") for n in ["A549", "Calu3", "NCI-H226", "NCI-H460", "PC-9"]},
    **{n: ("digestive", "high", "colorectal / pancreatic line") for n in
       ["Caco-2", "HCEC 1CT", "HCT116", "HT-29", "LoVo", "RKO", "SW480", "Panc1"]},
    **{n: ("kidney_urinary", "high", "kidney / urothelial line") for n in
       ["ACHN", "Caki2", "G401", "HEK293", "HEK293T", "HK-2", "RCC", "RCC 7860", "urothelium cell line"]},
    **{n: ("reproductive_breast", "high", "prostate / breast / cervix / endometrium line") for n in
       ["22Rv1", "C4-2B", "DU 145", "LNCAP", "LNCaP clone FGC", "PC-3", "RWPE1", "RWPE2", "VCaP",
        "MCF-7", "T47D", "MCF 10A", "HeLa-S3", "Ishikawa"]},
    "HTR-8/SVneo": ("embryonic_extraembryonic", "high", "trophoblast line"),
    **{n: ("skin", "high", "melanoma") for n in ["A375", "COLO829", "RPMI7951", "SK-MEL-5"]},
    **{n: ("stromal_connective_adipose", "high", "fibroblast / stromal line") for n in
       ["AG04449", "AG04450", "AG08395", "AG08396", "AG09309", "AG09319", "AG10803", "BJ", "GM03348",
        "GM04503", "GM04504", "GM23248", "HFF-Myc", "HFFc6", "IMR-90", "WI38", "skin fibroblast",
        "HS-27A", "HS-5"]},
    **{n: ("stromal_connective_adipose", "high", "bone / soft-tissue sarcoma") for n in
       ["HT1080", "MG63", "SJSA1", "A673"]},
    "SJCRH30": ("skeletal_muscle", "high", "rhabdomyosarcoma"),
    "LHCN-M2": ("skeletal_muscle", "high", "immortalised skeletal myoblast"),
    "AG20443": ("stromal_connective_adipose", "low", "Coriell AG line, assumed fibroblast: check"),
    "ACC112": ("other", "low", "adenoid cystic carcinoma (salivary): check"),
    **{n: ("other", "low", "origin not verified: check ENCODE biosample") for n in ["EH", "EL", "ELR", "ELF-1"]},
}


# Samples the ontology rules leave as "other" (Cell Ontology terms often lack part_of
# links to organs in the basic release). Applied only when the rules found nothing.
OVERRIDES = {
    **{n: "skeletal_muscle" for n in ["skeletal muscle satellite cell", "cell of skeletal muscle",
                                      "skeletal muscle cell", "myotube"]},
    **{n: "heart" for n in ["left ventricle myocardium inferior", "left ventricle myocardium superior",
                            "Right ventricle myocardium inferior", "Right ventricle myocardium superior",
                            "cardiac mesenchymal cell"]},
    **{n: "lung_airway" for n in ["bronchial epithelial cell", "tracheal epithelial cell", "airway epithelial cell",
                                  "epithelial cell of alveolus of lung", "respiratory epithelial cell",
                                  "nasal cavity respiratory epithelium epithelial cell of viscerocranial mucosa",
                                  "olfactory region"]},
    **{n: "embryonic_extraembryonic" for n in ["chorion", "amnion", "placental epithelial cell",
                                               "placental basal plate", "amniotic epithelial cell", "egg chorion"]},
    # Whole prenatal limbs (DNase tracks 301-304; UBERON:8300001-8300004 have no organ-system
    # anchor). Mixed embryonic tissue (muscle, bone, skin, vessels), so grouped with the other
    # embryonic samples rather than one organ system. Decision: Aryana, 2026-11-17.
    **{n: "embryonic_extraembryonic" for n in ["right forelimb", "left forelimb",
                                               "right hindlimb", "left hindlimb"]},
    "epithelial cell of umbilical artery": "vascular",
    "parathyroid adenoma": "endocrine",
    **{n: "stromal_connective_adipose" for n in [
        "fibroblast of skin of scalp", "fibroblast of skin of left biceps", "fibroblast of skin of left quadriceps",
        "fibroblast of skin of right quadriceps", "fibroblast of skin of right biceps", "fibroblast of breast",
        "fibroblast derived cell line", "synovial cell", "nucleus pulposus cell of intervertebral disc",
        "annulus pulposus cell", "cruciate ligament of knee", "femur", "mesothelial cell", "fat cell",
        "subcutaneous fat cell", "adipocyte of omentum tissue", "perirenal adipocyte cell", "adipocyte of breast"]},
    **{n: "eye_anterior" for n in ["non-pigmented ciliary epithelial cell", "iris pigment epithelial cell",
                                   "corneal epithelial cell", "lens epithelial cell", "trabecular meshwork cell"]},
    **{n: "digestive" for n in ["epithelial cell of Malassez", "gingival epithelial cell",
                                "acinar cell of salivary gland", "epithelial cell of esophagus",
                                "intestinal epithelial cell"]},
    **{n: "immune_blood" for n in ["lymphoblast", "bone marrow cell"]},
    **{n: "reproductive_breast" for n in ["mammary stem cell", "Sertoli cell"]},
    **{n: "kidney_urinary" for n in ["urothelial cell", "kidney tubule cell"]},
    **{n: "skin" for n in ["outer root sheath cell", "acinar cell of sebaceous gland", "stem cell of epidermis"]},
}


def line_origin(name):
    if name in CELL_LINE_ORIGIN:
        return CELL_LINE_ORIGIN[name]
    if (name.startswith("GM") or name.startswith("HG")) and name[2:].isdigit():
        return _LCL
    return None


def kind(btype):
    return {"TISSUE": "tissue", "PRIMARY_CELL": "primary_cell", "IN_VITRO_DIFFERENTIATED_CELLS": "in_vitro",
            "ORGANOID": "organoid", "CELL_LINE": "cell_line"}.get(btype, "unknown")


def main():
    for p in [CLASSES, *OBO]:
        if not p.exists():
            sys.exit(f"Not found: {p} (run build_track_classes.py first; ontologies in phase2/metadata/ontologies/)")
    names, parents = load_obo(*OBO)
    by_name = defaultdict(list)
    for tid, n in names.items():
        by_name[n].append(tid)
    rules, missing = [], []
    for cat, anchors in RULES:
        ids = set()
        for a in anchors:
            if len(by_name.get(a, [])) != 1:
                missing.append(f"{a!r} ({len(by_name.get(a, []))} terms)")
            ids.update(by_name.get(a, []))
        rules.append((cat, ids))
    if missing:
        sys.exit("Anchor terms not found exactly once in the ontologies: " + ", ".join(missing))

    rows = [r for r in csv.DictReader(open(CLASSES, encoding="utf-8"), delimiter="\t") if r["is_padding"] != "True"]
    out, unresolved = [], Counter()
    for r in rows:
        k = kind(r["biosample_type"])
        g = r["group"]
        conf, source = "high", ""
        if g in NEURAL_GROUP:
            cat, source = NEURAL_GROUP[g], f"group {g}"
        elif k == "cell_line":
            o = line_origin(r["biosample_name"])
            cat, conf, source = (o[0], o[1], o[2]) if o else ("other", "low", "cell line not in origin table")
        elif k in ("in_vitro", "organoid"):
            cat, source = "stem_cell_derived_other", "in vitro differentiated, non-neural"
        else:
            term = r["ontology_curie"]
            anc = ancestors(term, parents) | {term}
            cat = next((c for c, ids in rules if ids & anc), "other")
            hit = next((ids & anc for c, ids in rules if ids & anc), set())
            source = "ontology: " + ", ".join(sorted(names.get(t, t) for t in hit)) if hit else "no anchor matched"
            if cat == "other" and g == "CNS_non_neural":
                cat, source = "CNS_non_neural", "group CNS_non_neural (choroid plexus, meninges, CSF)"
            elif cat == "other" and r["biosample_name"] in OVERRIDES:
                cat, source = OVERRIDES[r["biosample_name"]], "reviewed override (no ontology anchor)"
        if cat == "other":
            unresolved[(r["biosample_type"], r["biosample_name"], r["ontology_curie"])] += 1
            conf = "low"
        out.append({"output_type": r["output_type"], "track_index": r["track_index"],
                    "biosample_name": r["biosample_name"], "biosample_type": r["biosample_type"],
                    "group": g, "analysis_set": r["analysis_set"], "life_stage_group": r["life_stage_group"],
                    "tissue_category": cat, "sample_kind": k,
                    "main_analysis": str(k in ("tissue", "primary_cell", "in_vitro", "organoid")),
                    "confidence": conf, "category_source": source})

    with open(META / "track_categories.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    heads = ["DNASE", "ATAC", "CHIP_HISTONE", "RNA_SEQ", "CAGE", "CHIP_TF"]
    cats = sorted({o["tissue_category"] for o in out}, key=lambda c: (c != "CNS", c))
    cnt = Counter((o["output_type"], o["tissue_category"], o["main_analysis"]) for o in out)
    with open(META / "track_categories.summary.tsv", "w", encoding="utf-8") as f:
        f.write("tissue_category\t" + "\t".join(f"{h}_main\t{h}_lines" for h in heads) + "\n")
        for c in cats:
            f.write(c + "\t" + "\t".join(f"{cnt[(h, c, 'True')]}\t{cnt[(h, c, 'False')]}" for h in heads) + "\n")
    print(f"{'category':<28}" + "".join(f"{h[:12]:>14}" for h in heads) + "     (main + cell lines)")
    for c in cats:
        print(f"{c:<28}" + "".join(f"{cnt[(h, c, 'True')]:>8} +{cnt[(h, c, 'False')]:>4}" for h in heads))
    print(f"\n{len(out)} real tracks -> phase2/metadata/track_categories.tsv")
    if unresolved:
        print(f"\n{sum(unresolved.values())} tracks in 'other' ({len(unresolved)} biosamples) - review:")
        for (bt, n, cu), c in sorted(unresolved.items(), key=lambda x: -x[1]):
            print(f"  {c:>3}  {bt:<14} {cu:<16} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())