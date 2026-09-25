"""Self-contained HTML report for RT-RubyRed (plotly.js inlined, works offline)."""
import datetime as dt
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


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def build(out_path, table, stats, meta, tool_dir):
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
        "meta": {
            "run": meta.get("run_name", ""),
            "dir": meta.get("dir_name", ""),
            "watch": meta.get("watch_dir", ""),
            "out": meta.get("output_dir", ""),
            "started": meta.get("started", ""),
            "updated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "files": sum(b.get("files", 0) for b in batches),
            "raw_total": sum(b.get("raw_total", 0) for b in batches),
        },
    }

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
