"""Self-contained HTML report for RT-RubyRed (plotly.js inlined, works offline)."""
import datetime as dt
import glob
import gzip
import html
import json
import os
import re

RANKS = {"K": "kingdom", "P": "phylum", "C": "class", "O": "order", "F": "family",
         "G": "genus", "S": "species"}


def taxon_label(lineage):
    """Species name if resolved, otherwise the lowest named rank."""
    named = []
    for part in lineage.split(";"):
        m = re.match(r"\s*([A-Za-z])__(.*)", part)
        if m and m.group(2).strip():
            named.append((m.group(1).upper(), m.group(2).strip().replace("_", " ")))
    if not named:
        return "Unassigned"
    rank, name = named[-1]
    if rank == "S":
        return name
    if rank == "G":
        return f"{name} sp."
    return f"{name} ({RANKS.get(rank, rank)})"


# RubyRed logs "Number of sequences <reason>: N" for every filtering step
FATE_KEYS = [("primer_length", r"primer|incorrect length|amplicon length"),
             ("quality", r"quality"),
             ("subsampling", r"subsampl|read number"),
             ("chimera", r"chimer"),
             ("unmatched", r"reorient|matched to database")]


def read_fate(logs_dir, tags):
    """Sum RubyRed's per-step removal counts over the batch logs (logs/batchNN.log)."""
    out = {}
    for tag in tags:
        try:
            with open(os.path.join(logs_dir, f"{tag}.log"), errors="replace") as fh:
                text = fh.read()
        except OSError:
            continue
        for reason, n in re.findall(r"Number of sequences (.+?):\s*(\d+)", text):
            key = next((k for k, rx in FATE_KEYS if re.search(rx, reason, re.I)), reason.strip())
            out[key] = out.get(key, 0) + int(n)
    return out


def minknow_info(watch_dir):
    """Sequencing details from the MinKNOW run folder (report_*.json, written at the end of
    the run) or, while the run is still going, from the first FASTQ header."""
    info = {}
    d = os.path.realpath(watch_dir or ".")
    run_dir = None
    for _ in range(4):
        if glob.glob(os.path.join(d, "report_*.json")) or glob.glob(os.path.join(d, "final_summary_*.txt")):
            run_dir = d
            break
        d = os.path.dirname(d)
    if run_dir:
        try:
            with open(sorted(glob.glob(os.path.join(run_dir, "report_*.json")))[0]) as fh:
                j = json.load(fh)
            pri = j.get("protocol_run_info", {})
            sw = pri.get("software_versions", {})
            fc = pri.get("flow_cell", {})
            dev = pri.get("device", {})
            acq = (j.get("acquisitions") or [{}])[-1]
            bc = acq.get("basecaller_config", {})
            info.update({
                "minknow": (sw.get("minknow") or {}).get("full", ""),
                "distribution": sw.get("distribution_version", ""),
                "dorado": (sw.get("basecaller_connected_version") or "").split("+")[0],
                "model": (bc.get("model_names") or {}).get("simplex_model", ""),
                "min_qscore": (bc.get("read_filtering") or {}).get("min_qscore", ""),
                "kits": ", ".join((bc.get("barcoding_configuration") or {}).get("barcoding_kits", [])),
                "flow_cell_id": fc.get("flow_cell_id", ""),
                "flow_cell": fc.get("product_code", ""),
                "device": dev.get("device_type", ""),
                "device_id": dev.get("device_id", ""),
            })
        except (OSError, ValueError, IndexError, AttributeError):
            pass
    if not info.get("model"):
        files = sorted(glob.glob(os.path.join(watch_dir or ".", "*.fastq.gz")))[:1]
        try:
            with gzip.open(files[0], "rt") as fh:
                head = dict(kv.split("=", 1) for kv in fh.readline().split()[1:] if "=" in kv)
            info.setdefault("model", head.get("basecall_model_version_id", ""))
            info.setdefault("flow_cell_id", head.get("flow_cell_id", ""))
        except (OSError, IndexError, ValueError, EOFError):
            pass
    return {k: v for k, v in info.items() if v not in ("", None)}


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def build(out_path, table, stats, meta, tool_dir, batch_tables=None):
    """batch_tables: paths of the per-batch feature tables; their text is embedded
    so the report can offer them (and the merged table) as downloads."""
    labels, seen = [], {}
    for lin in table.index:
        lab = taxon_label(lin)
        if lab in seen:
            seen[lab] += 1
            lab = f"{lab} [{seen[lab]}]"
        else:
            seen[lab] = 1
        labels.append(lab)

    batches = [b for b in stats.get("batches", [])]
    raw = {}
    for b in batches:
        for s, n in b.get("raw_reads", {}).items():
            raw.setdefault(s, {})[b["tag"]] = n
    samples = sorted(set(table.columns) | set(raw), key=natural_key)
    t = table.reindex(columns=samples, fill_value=0)

    data = {
        "samples": samples,
        "taxa": labels,
        "lineages": list(table.index),
        "counts": t.values.astype(int).tolist(),
        "raw": raw,
        "batches": [b["tag"] for b in batches],
        "batch_classified": {b["tag"]: int(b.get("classified_added", 0) or 0) for b in batches},
        "meta": {
            "run": meta.get("run_name", ""),
            "dir": meta.get("dir_name", ""),
            "watch": meta.get("watch_dir", ""),
            "out": meta.get("output_dir", ""),
            "started": meta.get("started", ""),
            "updated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "files": sum(b.get("files", 0) for b in batches),
            "raw_total": sum(b.get("raw_total", 0) for b in batches),
            "version": meta.get("version", ""),
            "merged_name": os.path.basename(out_path).replace("_report.html", "_feature_table.tsv"),
        },
        "tables": {},
        "fate": read_fate(os.path.join(os.path.dirname(os.path.abspath(out_path)), "logs"),
                          [b["tag"] for b in batches if b.get("status", "ok") == "ok"]),
        "failed_raw": sum(b.get("raw_total", 0) for b in batches if b.get("status", "ok") != "ok"),
        "pipeline": meta.get("pipeline", {}),
        "sequencing": minknow_info(meta.get("watch_dir", "")),
    }
    for path in batch_tables or []:
        try:
            with open(path, encoding="utf-8") as fh:
                data["tables"][os.path.basename(path)] = fh.read()
        except OSError:
            pass

    assets = os.path.join(tool_dir, "assets")

    def asset(name):
        with open(os.path.join(assets, name), encoding="utf-8") as fh:
            return fh.read()

    page = (asset("report_template.html")
            .replace("<!--__RUBYRED_LOGO__-->", strip_xml(asset("rubyred_logo.svg")))
            .replace("<!--__KU_LOGO__-->", strip_xml(asset("ku_seal.svg")))
            .replace("__TITLE__", html.escape(f"RT-RubyRed · {data['meta']['run']} · {data['meta']['dir']}"))
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
            .replace("/*__PLOTLY__*/", asset("plotly.min.js").replace("</script", "<\\/script")))
    tmp = out_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(page)
    os.replace(tmp, out_path)


def strip_xml(svg):
    return re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>|<!--.*?-->", "", svg, flags=re.S).strip()
