#!/usr/bin/env python3
"""
MetaChlam global report logic
============================

This script is the final reporting step for the MetaChlam Nextflow pipeline.
It aggregates the outputs from the four strain-classification tools
(StrainGE, StrainScan, Sourmash, and LINtax) into:

  1. a per-sample summary CSV
  2. a per-sample evidence TSV containing all tool-level hits

Important reporting principles
------------------------------

1. The evidence table preserves all hits
   Each classifier can report multiple strain/substrain hits for a sample.
   These are all retained in the evidence TSV rather than collapsed early.
   This is especially important for mixed-strain samples.

2. Ct read support comes from the Kraken2 Ct-only classification step
   Because the Kraken2 database used here is Ct-specific, the number of
   classified reads is interpreted as the number of Ct-classified reads.
   For paired-end data, this value is converted to read pairs by dividing by 2.

3. Ct detection and strain classification are reported separately
   - Ct_reads_detected indicates whether Kraken2 found Ct-classified reads.
   - Ct_classification_detected indicates whether any downstream classifier
     reported a recognized Ct strain or substrain.
   - Ct_detected is only set to 1 when both read-level and classifier-level
     evidence are present.

4. Consensus calls use a two-tool support rule
   For public metagenome reporting, MetaChlam conservatively calls:
     - single infection: exactly one strain supported by at least two tools
     - mixed infection: two or more strains each supported by at least two tools
     - low_support: Ct strain signal detected, but no strain reaches two-tool support
     - none: no Ct reads and no classifier calls
     - unclassified: Ct reads exist but no recognized strain call is produced

5. One-tool strain signals are not discarded
   Strains detected by only one classifier are not used to define the main
   infection type, but they are preserved as secondary signals in:
     - All_detected_strains
     - Strain_confidence
     - Secondary_strains
     - Has_secondary_strain_signal

6. Confidence strings report all evidence
   Strain_confidence and Substrain_confidence include every detected label and
   the tools supporting it, for example:
     NP-UA(3/4:LINtax,Sourmash,StrainGE); LGV(1/4:StrainScan)

7. Substrain calls are handled conservatively
   For NP-UA, P-UA, and Ocular calls, substrains are taken from parsed clade
   labels such as NP-UA_clade1 or P-UA_clade3.
   LGV subclasses are a special case: some tools report LGV subclass information
   only in Raw_hit, e.g. LGV__L2b__..., while Detected_substrain is blank.
   Therefore, LGV subclasses are inferred from Raw_hit only for LGV hits.

8. The summary is a report-level interpretation, not a replacement for evidence
   The summary provides conservative consensus calls for quick interpretation.
   The evidence TSV remains the source of truth for all raw tool outputs.

MetaChlam global report

Creates:
  - per-sample summary CSV (1 row per sample)
  - per-sample evidence TSV (multi-row; all hits from all tools)

Ct reads logic:
  - Kraken2 database is Ct-only => "sequences classified" == Ct reads (per-read)
  - paired-end => convert to read-pairs: Ct_reads = classified_reads // 2
  - Ct_reads_pct computed from (classified+unclassified)//2 when available; else fallback to Kraken %.

Tool agreement:
  - 4 tools: StrainGE, StrainScan, Sourmash, LINtax
  - Evidence file preserves all hits from all tools.
  - Summary file reports consensus calls supported by >=2 tools, all detected signals,
    secondary one-tool signals, and discordance/low-support flags.
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

CANONICAL_STRAINS = {
    "ocular": "ocular",
    "np_ua": "NP-UA",
    "p_ua": "P-UA",
    "lgv": "LGV",
}

TOOLS = ["StrainGE", "StrainScan", "Sourmash", "LINtax"]


def normalize_hit(s: str) -> str:
    if s is None:
        return ""
    x = str(s).strip().lower()
    x = x.replace("-", "_")
    x = x.replace("/", "_")
    x = re.sub(r"\s+", " ", x).strip()
    x = x.replace("__clade", "_clade")
    x = x.replace("-clade", "_clade")
    x = re.sub(r"_+", "_", x)
    return x


SUBSTRAIN_RE = re.compile(r"(ocular|np_ua|p_ua|lgv)_clade(\d+)(?=\b|_|$)", re.IGNORECASE)
STRAIN_RE = re.compile(r"\b(ocular|np_ua|p_ua|lgv)\b", re.IGNORECASE)


def infer_lgv_subclass_from_raw_hit(raw_hit: str) -> str:
    """Infer LGV subclass only from raw hits like LGV__L2b__..."""
    raw_hit = str(raw_hit or "").strip()
    m = re.match(r"LGV__(L\d+[a-z]?)__", raw_hit)
    if m:
        return f"LGV_{m.group(1)}"
    return ""


def parse_strain_substrain(norm: str, raw_hit: str = "") -> Tuple[Optional[str], Optional[str], str, str]:
    """
    Returns (strain, substrain, matched_pattern, substrain_source).
    LGV subclasses are inferred from Raw_hit only when standard clade parsing is absent.
    """
    if norm:
        m = SUBSTRAIN_RE.search(norm)
        if m:
            key = m.group(1).lower().replace("-", "_")
            clade_num = m.group(2)
            strain = CANONICAL_STRAINS.get(key, None)
            substrain = f"{strain}_clade{clade_num}" if strain else None
            return strain, substrain, "substrain", "detected_substrain_column"

    lgv_sub = infer_lgv_subclass_from_raw_hit(raw_hit)
    if lgv_sub:
        return "LGV", lgv_sub, "substrain", "raw_hit_lgv_subclass"

    if norm:
        m = STRAIN_RE.search(norm)
        if m:
            key = m.group(1).lower().replace("-", "_")
            strain = CANONICAL_STRAINS.get(key, None)
            return strain, None, "strain", "none"

    return None, None, "none", "none"


'''
CLASSIFIED_RE = re.compile(r"^\s*(\d+)\s+sequences classified\s+\(([\d\.]+)%\)", re.MULTILINE)
UNCLASSIFIED_RE = re.compile(r"^\s*(\d+)\s+sequences unclassified\s+\(([\d\.]+)%\)", re.MULTILINE)


def parse_kraken_log(log_path: str) -> Dict[str, object]:
    out = {"Ct_reads": 0, "Ct_reads_pct": 0.0, "Kraken_note": ""}
    if not log_path:
        out["Kraken_note"] = "kraken_log_missing"
        return out
    p = Path(log_path)
    if not p.exists() or p.stat().st_size == 0:
        out["Kraken_note"] = "kraken_log_empty"
        return out
    txt = p.read_text(errors="replace")
    mC = CLASSIFIED_RE.search(txt)
    if not mC:
        out["Kraken_note"] = "kraken_classified_not_found"
        return out
    classified_reads = int(mC.group(1))
    classified_pct = float(mC.group(2))
    ct_pairs = classified_reads // 2
    out["Ct_reads"] = ct_pairs
    mU = UNCLASSIFIED_RE.search(txt)
    if mU:
        unclassified_reads = int(mU.group(1))
        total_pairs = (classified_reads + unclassified_reads) // 2
        out["Ct_reads_pct"] = (ct_pairs / total_pairs) * 100.0 if total_pairs > 0 else 0.0
    else:
        out["Ct_reads_pct"] = classified_pct
        out["Kraken_note"] = "unclassified_not_found_used_pct_fallback"
    return out
'''
EXTRACT_PRINTED_RE = re.compile(r"^\s*(\d+)\s+reads printed to file", re.MULTILINE)

def parse_extract_read_log(log_path: str) -> Dict[str, object]:
    out = {"Ct_reads": 0, "Extract_note": ""}

    if not log_path:
        out["Extract_note"] = "extract_read_log_missing"
        return out

    p = Path(log_path)
    if not p.exists() or p.stat().st_size == 0:
        out["Extract_note"] = "extract_read_log_empty"
        return out

    txt = p.read_text(errors="replace")
    m = EXTRACT_PRINTED_RE.search(txt)

    if not m:
        out["Extract_note"] = "extract_reads_printed_not_found"
        return out

    out["Ct_reads"] = int(m.group(1))
    return out

def _row(tool: str, raw: str, hit_rank: int, evidence_value, evidence_metric: str) -> Dict[str, object]:
    norm = normalize_hit(raw)
    strain, substrain, pat, substrain_source = parse_strain_substrain(norm, raw_hit=raw)
    strain_source = "detected_strain_column" if strain else "none"
    if strain == "LGV" and substrain_source == "raw_hit_lgv_subclass":
        strain_source = "raw_hit_lgv_subclass"
    return {
        "Tool": tool,
        "Raw_hit": raw,
        "Normalized_hit": norm,
        "Detected_strain": strain or "",
        "Detected_substrain": substrain or "",
        "Detected_strain_source": strain_source,
        "Detected_substrain_source": substrain_source,
        "Hit_rank": hit_rank,
        "Evidence_value": evidence_value,
        "Evidence_metric": evidence_metric,
        "Matched_pattern": pat,
    }


def read_strainge(strainge_file: str) -> pd.DataFrame:
    rows = []
    tool = "StrainGE"
    if not strainge_file:
        return pd.DataFrame(rows)
    p = Path(strainge_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame(rows)
    name_col = next((c for c in ["strain", "reference", "ref", "name", "strain_name"] if c in df.columns), None)
    score_col = next((c for c in ["gcov", "cov", "coverage"] if c in df.columns), None)
    if name_col is None:
        return pd.DataFrame(rows)
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec[name_col])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_strainscan(strainscan_dir: str) -> pd.DataFrame:
    rows = []
    tool = "StrainScan"
    if not strainscan_dir:
        return pd.DataFrame(rows)
    p = Path(strainscan_dir) / "final_report.txt"
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame(rows)
    if "Strain_Name" not in df.columns:
        return pd.DataFrame(rows)
    score_col = "Predicted_Depth (Enet)" if "Predicted_Depth (Enet)" in df.columns else None
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["Strain_Name"])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_sourmash(sourmash_file: str) -> pd.DataFrame:
    rows = []
    tool = "Sourmash"
    if not sourmash_file:
        return pd.DataFrame(rows)
    p = Path(sourmash_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep=",")
    except Exception:
        return pd.DataFrame(rows)
    if "name" not in df.columns:
        return pd.DataFrame(rows)
    score_col = "f_orig_query" if "f_orig_query" in df.columns else None
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["name"])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_lintax(linreport_file: str) -> pd.DataFrame:
    rows = []
    tool = "LINtax"
    if not linreport_file:
        return pd.DataFrame(rows)
    p = Path(linreport_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep=",")
    except Exception:
        return pd.DataFrame(rows)
    required = {"LINgroup_Name", "Assigned_reads"}
    if not required.issubset(set(df.columns)):
        return pd.DataFrame(rows)
    if "Percentage_assigned_reads" in df.columns:
        df["Percentage_assigned_reads"] = pd.to_numeric(df["Percentage_assigned_reads"], errors="coerce").fillna(0.0)
        df = df[(df["Percentage_assigned_reads"] > 0)]
    df["Assigned_reads"] = pd.to_numeric(df["Assigned_reads"], errors="coerce").fillna(0).astype(int)
    df = df.sort_values("Assigned_reads", ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["LINgroup_Name"])
        rows.append(_row(tool, raw, i, int(rec["Assigned_reads"]), "Assigned_reads"))
    return pd.DataFrame(rows)


def detection_bin_level(ct_reads_pairs: int) -> Tuple[str, str]:
    if ct_reads_pairs <= 0:
        return "0", "none"
    if ct_reads_pairs < 500:
        return "<500", "low"
    if ct_reads_pairs < 1000:
        return "500-999", "confident"
    if ct_reads_pairs < 1500:
        return "1000-1499", "strain_conf"
    if ct_reads_pairs < 5000:
        return "1500-4999", "substrain_conf"
    return ">=5000", "very_high"


def build_conf_string(support_map: Dict[str, List[str]]) -> str:
    parts = []
    for label, tools in sorted(support_map.items(), key=lambda x: (-len(set(x[1])), x[0])):
        uniq_tools = sorted(set(tools))
        parts.append(f"{label}({len(uniq_tools)}/4:{','.join(uniq_tools)})")
    return "; ".join(parts) if parts else "-"


def join_or_dash(items: List[str]) -> str:
    return ";".join(items) if items else "-"


def main(extract_read_log: str, strainge_file: str, linreport_file: str, sourmash_file: str,
         strainscan_dir: str, sample_id: str, output_summary: str, output_evidence: str):
    kdat = parse_extract_read_log(extract_read_log)
    ct_reads = int(kdat["Ct_reads"])
    #ct_pct = float(kdat["Ct_reads_pct"])
    extract_note = str(kdat.get("Extract_note", "") or "")
    ct_bin, ct_level = detection_bin_level(ct_reads)
    ct_reads_detected = 1 if ct_reads > 0 else 0

    evidence_frames = [
        read_strainge(strainge_file),
        read_strainscan(strainscan_dir),
        read_sourmash(sourmash_file),
        read_lintax(linreport_file),
    ]
    frames = [f for f in evidence_frames if f is not None and not f.empty]
    if frames:
        ev = pd.concat(frames, ignore_index=True)
    else:
        ev_cols = [
            "Tool", "Raw_hit", "Normalized_hit", "Detected_strain", "Detected_substrain",
            "Detected_strain_source", "Detected_substrain_source", "Hit_rank", "Evidence_value",
            "Evidence_metric", "Matched_pattern"
        ]
        ev = pd.DataFrame({c: [] for c in ev_cols})
    ev.insert(0, "Sample", sample_id)

    tools_reporting = []
    strain_support: Dict[str, List[str]] = {}
    substrain_support: Dict[str, List[str]] = {}

    for t in TOOLS:
        sub = ev[ev["Tool"] == t]
        strains = set([x for x in sub["Detected_strain"].fillna("").tolist() if x])
        substrains = set([x for x in sub["Detected_substrain"].fillna("").tolist() if x])
        if strains or substrains:
            tools_reporting.append(t)
        for s in strains:
            strain_support.setdefault(s, []).append(t)
        for ss in substrains:
            substrain_support.setdefault(ss, []).append(t)

    tools_reporting_ct_n = len(tools_reporting)
    tools_reporting_ct = join_or_dash(tools_reporting)
    ct_classification_detected = 1 if tools_reporting_ct_n > 0 else 0
    ct_detected = 1 if (ct_reads_detected and ct_classification_detected) else 0

    all_detected_strains = sorted(strain_support.keys())
    all_detected_substrains = sorted(substrain_support.keys())
    consensus_strains = sorted([s for s, ts in strain_support.items() if len(set(ts)) >= 2])
    consensus_substrains = sorted([ss for ss, ts in substrain_support.items() if len(set(ts)) >= 2])
    secondary_strains = sorted([s for s, ts in strain_support.items() if len(set(ts)) == 1])

    has_secondary_strain_signal = 1 if secondary_strains else 0

    if ct_reads_detected == 0 and ct_classification_detected == 0:
        infection_type = "none"
    elif len(consensus_strains) == 1:
        infection_type = "single"
    elif len(consensus_strains) >= 2:
        infection_type = "mixed"
    elif len(all_detected_strains) > 0:
        infection_type = "low_support"
    else:
        infection_type = "unclassified"

    has_discordant_tool_calls = 1 if (len(all_detected_strains) >= 2 and len(consensus_strains) == 0) else 0

    notes = []
    if extract_note:
        notes.append(extract_note)
    if ct_reads_detected == 1 and ct_classification_detected == 0:
        notes.append("Ct_reads_present_but_no_tool_calls")
    if ct_reads_detected == 0 and ct_classification_detected == 1:
        notes.append("tools_report_hits_but_Ct_reads_zero_check_inputs")
    if has_secondary_strain_signal:
        notes.append("secondary_strain_signal_present")
    if has_discordant_tool_calls:
        notes.append("discordant_low_support_strain_calls")
    notes_str = ";".join(notes) if notes else "-"

    summary = pd.DataFrame([{
        "Sample": sample_id,
        "Ct_detected": ct_detected,
        "Ct_reads_detected": ct_reads_detected,
        "Ct_classification_detected": ct_classification_detected,
        "Ct_reads": ct_reads,
       # "Ct_reads_pct": round(ct_pct, 4),
        "Ct_detection_bin": ct_bin,
        "Ct_detection_level": ct_level,
        "Tools_reporting_Ct_n": f"{tools_reporting_ct_n}/4",
        "Tools_reporting_Ct": tools_reporting_ct,
        "Ct_infection_type": infection_type,
        "Consensus_strains": join_or_dash(consensus_strains),
        "Consensus_strain_count": len(consensus_strains),
        "All_detected_strains": join_or_dash(all_detected_strains),
        "Strain_confidence": build_conf_string(strain_support),
        "Consensus_substrains": join_or_dash(consensus_substrains),
        "Consensus_substrain_count": len(consensus_substrains),
        "All_detected_substrains": join_or_dash(all_detected_substrains),
        "Substrain_confidence": build_conf_string(substrain_support),
        "Has_secondary_strain_signal": has_secondary_strain_signal,
        "Secondary_strains": join_or_dash(secondary_strains),
        "Has_discordant_tool_calls": has_discordant_tool_calls,
        "Extract_note": extract_note if extract_note else "-",
        "Notes": notes_str,
    }])

    Path(output_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(output_evidence).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_summary, index=False)
    ev.to_csv(output_evidence, sep="\t", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MetaChlam: summarize tool outputs with transparent consensus reporting.")
    ap.add_argument("-el", "--extract_read_log", required=False, help="Path to Extract_Kraken_Reads.py log (stdout redirection).")
    ap.add_argument("-st", "--strainge", required=False, help="Path to StrainGE strains TSV.")
    ap.add_argument("-lr", "--linreport", required=False, help="Path to LINtax report CSV.")
    ap.add_argument("-sm", "--sourmash", required=False, help="Path to sourmash gather CSV.")
    ap.add_argument("-sc", "--strainscan", required=True, help="Path to StrainScan output directory (contains final_report.txt).")
    ap.add_argument("-s", "--sample", required=True, help="Sample ID.")
    ap.add_argument("-o", "--output", required=True, help="Summary CSV output path.")
    ap.add_argument("-e", "--evidence", required=True, help="Evidence TSV output path.")
    args = ap.parse_args()
    main(args.extract_read_log, args.strainge, args.linreport, args.sourmash, args.strainscan,
         args.sample, args.output, args.evidence)
