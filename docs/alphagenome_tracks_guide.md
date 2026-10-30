# AlphaGenome Tracks: Structure, Identification, and Brain vs iPSC Use

Oct 28, 2026

## Summary

A track is one output channel of AlphaGenome: the model's prediction of what one specific experiment from its training data would show for your DNA sequence. You identify a track by its row in the metadata table, and you tell brain from iPSC by the **biosample type** field, never by the sample name.

This doc answers three questions from Choo:

1. **What is a track and how is the output structured?** Sequence goes in; each output head returns a matrix of positions × tracks. Metadata row *i* describes column *i*.
2. **How do I tell which track is which (e.g. brain vs iPSC)?** Read four fields: ontology term, biosample type, biosample name, life stage. Real brain is TISSUE or PRIMARY\_CELL from the central nervous system. Stem-cell-derived neural cells are IN\_VITRO\_DIFFERENTIATED\_CELLS, often with the same name as a real brain cell type.
3. **How do we extract and use them?** Join scores to `phase2/metadata/track_classes.tsv` on `(output_type, track_index)`, drop padding, then filter on `analysis_set`.

For the project: keep **brain** (donor tissue and cells) as the primary set, and analyse **neural\_model** (ESC/iPSC-derived) as a separate, secondary set. Never pool the two. Section 5 explains why.

## 1. What a track is

AlphaGenome never takes experimental data as input. It reads about 1 Mb of DNA sequence and predicts thousands of experiments at once; each predicted experiment is a **track**.

**Where tracks come from.** DeepMind trained the model on real experiments from ENCODE, FANTOM5, GTEx and 4D Nucleome. Each experiment is one **biosample** (a tissue, cell type or cell line) measured with one **assay** (DNase-seq, ATAC-seq, RNA-seq, a histone ChIP and so on). The model learned to reproduce every one of those experiments from sequence alone. A track is the output slot reserved for one of them, for example "DNase-seq in putamen".

**How the output is laid out.** Think of each output head as a spreadsheet:

- **One head per assay type.** ATAC, DNase, CAGE, RNA-seq, histone ChIP, TF ChIP, PRO-cap, three splicing heads and contact maps.
- **Rows are positions along the input sequence.** At 1 bp resolution a 1,048,576 bp input gives 1,048,576 rows; at 128 bp resolution it gives 8,192.
- **Columns are tracks.** The DNase head has 384 columns, the histone head 1,152.
- **The metadata is the column header.** Metadata row *i* describes column *i* of that head. Nothing in the numbers themselves says which experiment a column is; only the metadata does.

**What the values mean.** Each value is the predicted signal for that experiment at that position, in the assay's own units (the AlphaGenome paper's Methods give each one). Values are comparable within a track, such as reference vs alternate allele, but not across heads.

**Strand.** Assays that read one DNA strand at a time (CAGE, RNA-seq, PRO-cap) have a `+` and a `-` track per biosample. DNase, ATAC and ChIP are unstranded (`.`).

**Padding.** Most heads have extra all-zero columns at the end so the width is a multiple of 128. They are named `Padding` and have no biosample. DNase has 305 real tracks and 79 padding columns, so 384 in total. Padding must be dropped before any average or ranking.

**What ISM does to this layout.** Our variant score sums each track's signal over a window around the variant, for the reference and the alternate sequence, and takes log2\[(Σ alt + 1) / (Σ ref + 1)\]. Positions collapse, so each variant ends up with **one score per track**: 305 DNase scores, 1,116 histone scores, and so on. Every brain-specificity question in this project is a question about which of those columns to compare.

## 2. Reading a track's metadata

Every real track carries the same handful of fields, and four of them identify the biology: ontology term, biosample type, biosample name and life stage. The exact column names differ slightly between the DeepMind API, the PyTorch port and our parsed table, so all three are given.

| What it tells you | DeepMind API column | PyTorch port | Our `track_classes.tsv` | Example |
| --- | --- | --- | --- | --- |
| Track label | `name` | `track_name` | `track_name` | CL:0000100 ATAC-seq |
| Standard ID of the sample | `ontology_curie` | `ontology_curie` | `ontology_curie` | CL:0000100 |
| **How the sample was obtained** | `biosample_type` | `biosample_type` | `biosample_type` | in\_vitro\_differentiated\_cells |
| Human-readable name | `biosample_name` | `biosample_name` | `biosample_name` | motor neuron |
| **Life stage of the donor** | `biosample_life_stage` | biosample\_life\_stage | `stage` | adult |
| Assay | `Assay title` | assay\_title | `assay` | ATAC-seq |
| Histone mark / TF | `histone_mark`, `transcription_factor` | same | `histone_mark`, `tf` | H3K27ac, CTCF |
| Consortium | `data_source` | same | `data_source` | encode |
| GTEx tissue (RNA-seq, splicing) | `gtex_tissue` | same | `gtex_tissue` | Brain\_Putamen\_basal\_ganglia |
| Strand | `strand` | `strand` | `strand` | + / - / . |

That example row is real: ATAC track 1 is motor neurons made in a dish, from an adult donor's cells. Its name says "neuron"; its type says it is not brain tissue.

**The ontology prefix says what kind of thing the sample is.**

| Prefix | Ontology | Describes | Example |
| --- | --- | --- | --- |
| UBERON | Uberon anatomy | a tissue or organ | UBERON:0001874 putamen |
| CL | Cell Ontology | a cell type | CL:0000127 astrocyte |
| EFO / CLO | Experimental Factor / Cell Line Ontology | a named cell line | SK-N-SH, H9 |
| NTR | ENCODE "new term request" | a provisional ENCODE term | germinal matrix, neurosphere |

**The biosample type says where the cells came from.** This is the field that separates brain from iPSC.

| Biosample type | Meaning | Brain-relevant example |
| --- | --- | --- |
| TISSUE | a piece of a person's organ, post-mortem or surgical | dorsolateral prefrontal cortex |
| PRIMARY\_CELL | cells isolated from a person's tissue and grown briefly | astrocyte of the cerebellum |
| IN\_VITRO\_DIFFERENTIATED\_CELLS | cells **made in a dish** by differentiating stem cells | neural progenitor cell, glutamatergic neuron |
| ORGANOID | 3-D tissue grown from stem cells | brain organoid |
| CELL\_LINE | an immortal line: a tumour line, or an ESC / iPSC line itself | SK-N-SH, H9, WTC11 |

**The life stage is recorded per track, not per sample.** For TISSUE and PRIMARY\_CELL it is the donor's age group: `embryonic` means prenatal. For lab-made cells it is the age of the person whose cells started the parent stem cell line (Section 3 uses this).

## 3. How to tell brain from iPSC

Read the biosample type first, the ontology term second and the life stage third; use the name only to find candidates. The name cannot decide, because stem-cell-derived cells carry the same names as real brain cells ("neural progenitor cell" exists as both).

&#91;embedded content: brain vs stem-cell-derived decision tree · 3 sample origins\]

The sample's origin decides the branch; the ontology term decides brain; for dish-made cells, the life stage separates ESC-derived from iPSC-derived.

**"iPSC" can mean three different kinds of track.** It helps to say which one you mean:

- **The stem cell line itself** (WTC11, iPS DF 19.11, H9). Biosample type CELL\_LINE. Not neural at all.
- **Neural cells made from a stem cell line** (glutamatergic neurons, neural progenitors, motor neurons). Biosample type IN\_VITRO\_DIFFERENTIATED\_CELLS, or ORGANOID for 3-D cultures. This is what people usually mean by "iPSC neurons".
- **ESC-derived vs iPSC-derived.** Both are made in a dish. The metadata does not name the parent line, but the life stage gives it away (below).

**The procedure:**

1. **Drop padding** (name `Padding`, empty biosample).
2. **TISSUE or PRIMARY\_CELL** = material taken from a person. It is **brain** if the ontology term is part of the central nervous system (Uberon) or is a CNS cell type such as an astrocyte. Exception: cells that sit in the brain but are not neural (blood-vessel cells, choroid plexus, meninges, pituitary) are excluded.
3. **IN\_VITRO\_DIFFERENTIATED\_CELLS or ORGANOID** with a neural cell type = **stem-cell-derived neural**. Then read the stage:
   - `embryonic` → made from an **embryonic** stem cell line (H1, H9)
   - `child` or `adult` → made from an **induced** pluripotent stem cell line, reprogrammed from that person's cells
4. **CELL\_LINE** = look up the line. ESC/iPSC lines are pluripotent, not neural. Tumour lines are excluded (glioblastoma, medulloblastoma, and neuroblastoma, which is peripheral, not brain).
5. **For brain tracks, read the stage again** to split prenatal (`embryonic`) from postnatal (`child`, `adult`).

The ESC-vs-iPSC rule is an inference from how ENCODE records donors, confirmed on two tracks: ENCODE's motor neurons are iPSC-derived from a 49-year-old donor (ENCBS745APJ), and its glutamatergic neurons come from the iPSC line WTC11, from a 30-year-old donor (ENCBS229NUR). For any single track, search the biosample on encodeproject.org and read "originated from".

**Worked examples (real tracks):**

| Biosample name | Head | Biosample type | Stage | Verdict |
| --- | --- | --- | --- | --- |
| brain | DNase | TISSUE | embryonic | **brain**, prenatal |
| astrocyte of the cerebellum | DNase | PRIMARY\_CELL | unknown | **brain** (cultured donor cells) |
| glutamatergic neuron | DNase | IN\_VITRO\_DIFFERENTIATED\_CELLS | adult | **iPSC-derived** neural (WTC11) |
| neural progenitor cell | DNase | IN\_VITRO\_DIFFERENTIATED\_CELLS | embryonic | **ESC-derived** neural |
| motor neuron | ATAC | IN\_VITRO\_DIFFERENTIATED\_CELLS | adult | **iPSC-derived** neural |
| WTC11, H9 | various | CELL\_LINE | — | pluripotent line, not neural |
| SK-N-SH | DNase | CELL\_LINE | — | neuroblastoma (peripheral tumour), excluded |
| brain microvascular endothelial cell | DNase | PRIMARY\_CELL | — | in brain but vascular, excluded |

**Traps found while building the labels:**

- **Name matching gives false hits.** "gastrocnemius" contains "astro", "renal cortex" contains "cortex", "Peyer's patch" contains "eye".
- **A cell-type name carries no location.** In the Cell Ontology, "neuron" and "astrocyte" are not linked to the CNS, and "neural progenitor cell" is not even classed as a neural cell.
- **The same name appears at different stages.** "brain" DNase is embryonic, while "brain" histone tracks include adult donors. Read stage per track.
- **Retina is not CNS in Uberon**; it is filed under the eye. We keep it as its own `eye` group.
- **Some cell lines are mislabelled.** SK-N-MC is an Askin tumour, not a neuroblastoma (Cellosaurus CVCL\_0530). H54 is a contaminated copy of A172 (CVCL\_7185).

## 4. Extracting and using tracks

In this project, the labels live in one file, `phase2/metadata/track_classes.tsv`: one row per output channel, keyed by `(output_type, track_index)`. Everything downstream should select tracks through it rather than re-deriving labels. It is built by `phase2/build_track_classes.py` from the DeepMind metadata, the Uberon and Cell Ontology files, and the reviewed decisions in `biosample_curation.tsv`.

**Columns that matter for selection:**

| Column | Values | Use it to |
| --- | --- | --- |
| `analysis_set` | brain, neural\_model, background, excluded, padding | pick the comparison groups |
| `group` | CNS\_tissue, CNS\_primary\_cell, neural\_in\_vitro, CNS\_organoid, eye, PNS, … | finer splits or sensitivity checks |
| `origin` | donor\_tissue, donor\_cells, ESC\_derived, iPSC\_derived, tumor\_line, … | separate iPSC from ESC |
| `life_stage_group` | prenatal, perinatal, postnatal, unknown, n/a | fetal vs adult brain |
| `histone_mark` | H3K27ac, H3K4me3, H3K27me3, … | analyse one mark at a time |
| `brain_region` | cerebral\_cortex, basal\_ganglia, cerebellum, … | regional breakdowns |

**A. Selecting tracks for ISM scores (our pipeline).** Join, check nothing failed to match, drop padding, then filter:

```python
import pandas as pd

tc = pd.read_csv("phase2/metadata/track_classes.tsv", sep="\t")
scores = pd.read_csv("results/pilot_ism_scores.tsv", sep="\t")

# 'modality' in the score file must use the same spelling as 'output_type' (e.g. DNASE)
m = scores.merge(tc, left_on=["modality", "track_index"],
                 right_on=["output_type", "track_index"],
                 how="left", validate="many_to_one")
assert m["analysis_set"].notna().all(), "score rows with no metadata: index mapping broke"
m = m[m["analysis_set"] != "padding"]

brain        = m[m.analysis_set == "brain"]
neural_model = m[m.analysis_set == "neural_model"]
background   = m[m.analysis_set == "background"]

ipsc_neural  = neural_model[neural_model.origin == "iPSC_derived"]
fetal_brain  = brain[brain.life_stage_group == "prenatal"]
brain_k27ac  = brain[(brain.output_type == "CHIP_HISTONE") & (brain.histone_mark == "H3K27ac")]
```

**B. Slicing model output directly (PyTorch port).** Named outputs strip padding by default. Because padding sits at the end of each head, real track *i* keeps index *i*, so our table's indices still apply. Check that the names line up before trusting it:

```python
from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.named_outputs import TrackMetadataCatalog

model = AlphaGenome.from_pretrained("weights.pth", device="cuda")
model.set_track_metadata_catalog(TrackMetadataCatalog.load_builtin("human"))
out = model.predict(dna_onehot, organism_index=0, named_outputs=True)

idx = tc.query("output_type == 'DNASE' and analysis_set == 'brain'").track_index.tolist()
meta = out.dnase.to_dataframe()
assert (meta.iloc[idx].track_name.values ==
        tc.query("output_type == 'DNASE'").set_index("track_index").loc[idx, "track_name"].values).all()

brain_dnase = out.dnase[128].tensor[..., idx]   # positions x 16 brain tracks
```

The port can also filter by metadata, e.g. `out.chip_histone[128].select(histone_mark="H3K27ac")` or `out.select(biosample_type="tissue")`. That is fine for exploring, but `biosample_type` alone does not mean brain, so use `analysis_set` for anything that goes in the thesis.

**C. Browsing metadata with the DeepMind API.** Useful for looking things up, e.g. every stem-cell-derived track:

```python
from alphagenome.models import dna_client

md = dna_model.output_metadata(dna_client.Organism.HOMO_SAPIENS).concatenate()
md[md.biosample_type == "in_vitro_differentiated_cells"][
    ["output_type", "biosample_name", "biosample_life_stage", "ontology_curie"]]
```

**Changing a label.** Edit the row's `group` in `biosample_curation.tsv` and rerun `python phase2/build_track_classes.py`. The build stops if the ontology finds a neural sample nobody has reviewed, and its log records the file checksums and ontology versions, so every result can be traced to the labels it used.

## 5. Using brain vs iPSC tracks in this project

Recommendation: test schizophrenia variants against **brain** tracks as the primary analysis, run the same test on **neural\_model** tracks as a separate secondary analysis, and never pool them. The two sets answer different questions, and the literature gives a specific reason to keep the second one rather than discard it.

**Why keep stem-cell-derived tracks at all.** Schizophrenia common-variant risk is enriched in genes of developing neurons in the second-trimester fetal brain (Cameron et al. 2023). Stem-cell-derived neurons look fetal rather than adult: their DNA methylation predicts an early fetal age (Steg et al. 2021), and their 3-D chromatin resembles fetal neurons (Zagirova et al. 2026). Our donor brain set is mostly postnatal, so the stem-cell set partly covers the developmental window the donor set misses.

**Why not count them as brain.** They are few, they are specific cell types grown under specific protocols, and they keep some chromatin features of the undifferentiated stem cell (Zagirova et al. 2026). Counting them as brain also changes results: in the August list, 4 of the 23 "CNS" DNase tracks were stem-cell-derived, and 3 of rs1198588's top 8 "CNS" tracks were among them.

|  | brain | neural\_model |
| --- | --- | --- |
| What it is | tissue or cells from human CNS donors | neural cells made from ESC or iPSC lines |
| Developmental stage | mostly postnatal | fetal-like |
| Tracks: DNase / histone / ATAC | 16 / 70 / 0 | 4 / 32 / 1 |
| Prenatal tracks: DNase / histone | 2 / 21 | not applicable (made in a dish) |
| Strength | real human brain, with brain regions | developmental proxy, defined neuron types |
| Weakness | few prenatal tracks; bulk tissue mixes cell types | few tracks; culture artefacts |
| Role here | primary test | secondary test, reported separately |

**Analysis plan, per variant and per modality:**

1. **Primary:** brain vs background, using DNase (16 tracks) and histone ChIP **one mark at a time** (start with H3K27ac and H3K4me3). ATAC cannot be used: it has no brain tracks.
2. **Secondary:** neural\_model vs background, the same statistic, in its own column.
3. **Concordance:** does a variant that is brain-specific also stand out in neural\_model? Agreement strengthens the call. Brain-only may point to mature or glial biology; neural\_model-only may point to a developmental effect. Treat both as hypotheses, not findings.
4. **Fetal sensitivity check:** repeat the primary test on prenatal brain histone tracks (21). DNase has only 2 prenatal brain tracks, so report those descriptively.
5. **Origin as a label, not a test:** carry `origin` (iPSC\_derived vs ESC\_derived) in the output tables. With 4 DNase neural\_model tracks, splitting them further would not support a statistical claim.
6. **Same controls for both:** the random-SNP null distribution must be computed separately for brain and for neural\_model, because the two sets differ in size.

The longer-term fix for the fetal gap is the fine-tuning arm on fetal brain ATAC-seq, which would add real prenatal brain tracks rather than proxies.

**Suggested wording for the methods:** "Brain tracks" are TISSUE or PRIMARY\_CELL samples from the human central nervous system. "Stem-cell-derived neural tracks" are IN\_VITRO\_DIFFERENTIATED\_CELLS or ORGANOID samples with a neural cell type, analysed separately as a model of developing neurons.

## 6. Reference counts and open questions

Counts come from the October 21 build of `track_classes.tsv` against the DeepMind human metadata file.

**Channels per head:**

| Head | Channels | Real tracks | Padding |
| --- | --- | --- | --- |
| ATAC | 256 | 167 | 89 |
| DNase | 384 | 305 | 79 |
| CAGE | 640 | 546 | 94 |
| RNA-seq | 768 | 667 | 101 |
| Histone ChIP | 1,152 | 1,116 | 36 |
| TF ChIP | 1,664 | 1,617 | 47 |
| PRO-cap | 128 | 12 | 116 |
| Contact maps | 28 | 28 | 0 |

**Real tracks per analysis set:**

| analysis\_set | ATAC | DNase | Histone | RNA-seq | CAGE | TF ChIP |
| --- | --- | --- | --- | --- | --- | --- |
| brain | 0 | 16 | 70 | 41 | 66 | 4 |
| neural\_model | 1 | 4 | 32 | 9 | 0 | 14 |
| excluded | 3 | 30 | 128 | 43 | 36 | 134 |
| background | 163 | 255 | 886 | 574 | 444 | 1,465 |

Brain tracks by donor stage: DNase 2 prenatal, 10 postnatal, 4 unknown (cultured astrocytes); histone 21 prenatal, 37 postnatal, 12 unknown. The 70 brain histone tracks span 11 marks, including 10 H3K27ac, 15 H3K4me3 and 6 H3K27me3.

**Open questions for the lab:**

1. Should stem-cell-derived neural tracks count as brain? Proposed answer: no; keep them as a separate secondary set (Section 5).
2. Should brain tumour lines stay excluded, or be kept as their own group?
3. Should FANTOM primary neurons and neuronal stem cells, whose location is not given, count as brain? They appear only in CAGE.
4. Should prenatal and postnatal brain be analysed separately, given the fetal-brain focus? DNase has only 2 prenatal brain tracks.
5. Is pituitary gland (mostly endocrine) correctly excluded?

## Sources

**AlphaGenome**

- [Model output metadata](https://google-deepmind-alphagenome.readthedocs-hosted.com/exploring_model_metadata.html), AlphaGenome documentation
- [Navigating data ontologies](https://google-deepmind-alphagenome.readthedocs-hosted.com/colabs/tissue_ontology_mapping.html), AlphaGenome tutorial (metadata columns, track counts)
- [Named Outputs](https://alphagenome-pytorch.readthedocs.io/en/latest/_sources/named_outputs.rst), alphagenome-pytorch documentation (fields, `select`, padding)
- Track metadata file: `OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto`, google-deepmind/alphagenome\_research on GitHub

**Sample identity**

- ENCODE biosamples: [motor neuron ENCBS745APJ](https://www.encodeproject.org/biosamples/ENCBS745APJ/), [glutamatergic neuron ENCBS229NUR](https://www.encodeproject.org/biosamples/ENCBS229NUR/)
- Cellosaurus: [SK-N-MC CVCL\_0530](https://www.cellosaurus.org/CVCL_0530), [H54 CVCL\_7185](https://www.cellosaurus.org/CVCL_7185)
- Uberon (release 2026-10-01) and Cell Ontology (release 2026-06-08), obophenotype on GitHub

**Brain vs stem-cell-derived neurons**

- [Cameron et al. 2023, Biological Psychiatry](https://orca.cardiff.ac.uk/151244/4/PIIS0006322322014044.pdf): schizophrenia risk in developing neurons of the second-trimester fetal brain
- [Steg et al. 2021, Molecular Brain](https://pmc.ncbi.nlm.nih.gov/articles/PMC8236187): iPSC-derived neurons have an early fetal epigenetic age
- [Zagirova et al. 2026, Genome Research](https://genome.cshlp.org/content/early/2026/07/31/gr281488125): iPSC-derived neurons resemble fetal neurons in 3-D chromatin and keep features of undifferentiated cells
