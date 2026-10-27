#!/usr/bin/env python3
"""Write the gene functional table for one 13_final gene set.

A final transcript whose protein is identical to a step-09 protein reuses that
protein's eggNOG, InterProScan, SwissProt, ZH13 v3 and Wm82 a6 results. Proteins
changed in step 10 (extension, merge) are annotated again with the same tools
and databases. Gene rows follow 13_final; removed and merged IDs do not appear.
"""

import argparse
import os
import re
import subprocess
from collections import defaultdict

from antisense_qc import attr, load_fasta, read_table, write_table

ENV_ANNP = "/home/zzh/anaconda3/envs/annp/bin"
ENV_EGGNOG = "/home/zzh/anaconda3/envs/eggnog213/bin"
EMAPPER = "/tiandata2/zzh/tools/eggnog-mapper-2.1.13/emapper.py"
EGGNOG_DATA = "/tiandata2/zzh/tools/eggnog-mapper-2.1.13/data"
INTERPROSCAN = "/home/zzh/software/interproscan/interproscan-5.77-108.0/interproscan.sh"
SWISSPROT_DMND = "/tiandata2/zzh/tools/uniprot_sprot.dmnd"
REFS = (
    ("ZH13v3", "ZH13v3.0_mRNA_nuclear_20250729_pro",
     "/tiandata2/zzh/tools/ref_diamond/ZH13v3.0_mRNA_nuclear_20250729_pro.dmnd"),
    ("Wm82a6", "glyma.Wm82.gnm6.ann1.PKSW.protein",
     "/tiandata2/zzh/tools/ref_diamond/glyma.Wm82.gnm6.ann1.PKSW.protein.dmnd"),
)
STAMP = "功能注释表已写入"


def clean(value):
    value = (value or "").strip()
    return "" if value in ("", "-", "NA") else value


def go_terms(text):
    return {item.split("(")[0] for item in re.split(r"[|,]", text or "") if item.strip().startswith("GO:")}


def swissprot_product(stitle):
    text = stitle.split(" ", 1)[1] if " " in stitle else ""
    return text.split(" OS=")[0].strip()


def ready(path, source):
    return os.path.isfile(path) and os.path.getsize(path) > 0 and os.path.getmtime(path) >= os.path.getmtime(source)


def write_if_changed(path, records):
    text = "".join(f">{name}\n{seq}\n" for name, seq in records)
    if os.path.isfile(path):
        with open(path) as handle:
            if handle.read() == text:
                return
    with open(path, "w") as handle:
        handle.write(text)


EGGNOG_SUFFIXES = (".emapper.annotations", ".emapper.hits", ".emapper.seed_orthologs")


def tool_paths(downstream, prefix):
    paths = {"swissprot": os.path.join(downstream, "05_diamond", f"{prefix}.final_novel.vs_swissprot.tsv")}
    for label, _name, _db in REFS:
        paths[label] = os.path.join(downstream, "05_diamond", f"{prefix}.final_novel.vs_{label}.tsv")
    paths["eggnog"] = os.path.join(downstream, "06_eggnog", f"{prefix}_final_novel.emapper.annotations")
    paths["interpro"] = os.path.join(downstream, "07_interpro", f"{prefix}.final_novel.interproscan.tsv")
    return paths


def run_tools(query, paths, threads):
    diamond_format = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore"]
    jobs = [("swissprot", SWISSPROT_DMND, diamond_format + ["stitle"])]
    jobs += [(label, db, diamond_format) for label, _name, db in REFS]
    for key, db, fields in jobs:
        if not ready(paths[key], query):
            subprocess.run([
                f"{ENV_ANNP}/diamond", "blastp", "-q", query, "-d", db, "-o", paths[key], "-f", "6", *fields,
                "--evalue", "1e-5", "--max-target-seqs", "5", "--threads", str(threads), "--quiet",
            ], check=True)
    if not ready(paths["eggnog"], query):
        out_dir = os.path.dirname(paths["eggnog"])
        name = os.path.basename(paths["eggnog"])[:-len(EGGNOG_SUFFIXES[0])]
        env = dict(os.environ, PATH=f"{ENV_EGGNOG}:{ENV_ANNP}:{os.environ['PATH']}")
        subprocess.run([
            EMAPPER, "-i", query, "--itype", "proteins", "--data_dir", EGGNOG_DATA,
            "-o", name, "--output_dir", out_dir, "--temp_dir", out_dir,
            "--cpu", str(threads), "--override",
        ], check=True, env=env)
    if not ready(paths["interpro"], query):
        # FunFam searches one superfamily at a time on a single CPU, so one
        # InterProScan over every novel protein takes days. Splitting the FASTA
        # runs the same analyses side by side. The TSV is concatenated.
        run_interpro(query, paths["interpro"], threads)


def run_interpro(query, dest, threads):
    records = []
    name = None
    buf = []
    with open(query) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            records.append((name, "".join(buf)))
    n = min(32, threads, len(records))
    work = dest + ".chunks"
    os.makedirs(work, exist_ok=True)
    groups = [[] for _ in range(n)]
    for index, record in enumerate(records):
        groups[index % n].append(record)
    procs = []
    for index, group in enumerate(groups):
        fasta = os.path.join(work, f"chunk_{index:02d}.faa")
        with open(fasta, "w") as handle:
            for rec_name, seq in group:
                handle.write(f">{rec_name}\n{seq}\n")
        base = os.path.join(work, f"chunk_{index:02d}")
        os.makedirs(base + ".tmp", exist_ok=True)
        log = open(os.path.join(work, f"chunk_{index:02d}.log"), "w")
        procs.append((index, subprocess.Popen(
            [INTERPROSCAN, "-i", fasta, "-f", "tsv", "-dp", "-goterms", "-pa",
             "-cpu", "2", "-b", base, "-T", base + ".tmp"],
            stdout=log, stderr=subprocess.STDOUT,
        ), log))
    print(f"interproscan chunks {len(procs)}", flush=True)
    failed = []
    for index, proc, log in procs:
        code = proc.wait()
        log.close()
        if code:
            failed.append(index)
    if failed:
        raise SystemExit(f"InterProScan failed for chunks {failed}; logs in {work}")
    with open(dest + ".tmp", "w") as out:
        for index, _proc, _log in procs:
            tsv = os.path.join(work, f"chunk_{index:02d}.tsv")
            if os.path.isfile(tsv):
                with open(tsv) as handle:
                    out.write(handle.read())
    os.replace(dest + ".tmp", dest)
    subprocess.run(["rm", "-rf", work], check=True)


def split_outputs(combined, targets):
    """Route each result line to the sample whose prefix starts the query ID."""
    files = [(key, "") for key in ("swissprot", *[label for label, _n, _d in REFS], "interpro")]
    files += [("eggnog", suffix) for suffix in EGGNOG_SUFFIXES]
    for key, suffix in files:
        src = combined[key][:-len(EGGNOG_SUFFIXES[0])] + suffix if suffix else combined[key]
        handles = {}
        for prefix, paths in targets.items():
            dest = paths[key][:-len(EGGNOG_SUFFIXES[0])] + suffix if suffix else paths[key]
            handles[prefix] = open(dest + ".tmp", "w")
        with open(src) as handle:
            for line in handle:
                if line.startswith("#"):
                    for out in handles.values():
                        out.write(line)
                    continue
                prefix = line.split("\t", 1)[0].split("_", 1)[0]
                if prefix in handles:
                    handles[prefix].write(line)
        for out in handles.values():
            out.close()
            os.replace(out.name, out.name[:-len(".tmp")])


def first_hits(path, wanted, stitle=False):
    hits = {}
    if not os.path.isfile(path):
        return hits
    with open(path) as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if p[0] in wanted and p[0] not in hits and len(p) >= 12:
                hits[p[0]] = (p[1], p[2], p[12] if stitle and len(p) > 12 else "")
    return hits


def best_hit_table(path, wanted):
    hits = {}
    with open(path) as handle:
        handle.readline()
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if p[0] in wanted:
                hits[p[0]] = (p[1], p[2], "")
    return hits


def eggnog_table(path, wanted):
    rows = {}
    header = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            p = line.rstrip("\n").split("\t")
            if header and len(p) == len(header) and p[0] in wanted:
                rows[p[0]] = dict(zip(header, p))
    return rows


def interpro_table(path, wanted):
    rows = defaultdict(lambda: {"ipr": set(), "pfam": set(), "go": set()})
    with open(path) as handle:
        for line in handle:
            q = line.split("\t", 1)[0]
            if q not in wanted:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 11 and clean(p[11]):
                rows[q]["ipr"].add(p[11])
            if len(p) > 4 and p[3] == "Pfam" and clean(p[4]):
                rows[q]["pfam"].add(p[4])
            if len(p) > 13:
                rows[q]["go"] |= go_terms(p[13])
    return rows


def te_overlap(bed_path, genes):
    by_chrom = defaultdict(list)
    for gene_id, (chrom, start, end) in genes.items():
        by_chrom[chrom].append((start - 1, end, gene_id))
    total = defaultdict(int)
    classes = defaultdict(set)
    with open(bed_path) as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if p[0] not in by_chrom:
                continue
            s, e = int(p[1]), int(p[2])
            for gs, ge, gene_id in by_chrom[p[0]]:
                bp = min(e, ge) - max(s, gs)
                if bp > 0:
                    total[gene_id] += bp
                    classes[gene_id].add(p[4])
    out = {}
    for gene_id, (chrom, start, end) in genes.items():
        length = end - start + 1
        out[gene_id] = (f"{total[gene_id] / length:.4f}", ",".join(sorted(classes[gene_id])) or "NA")
    return out


def prepare(downstream, prefix):
    final = os.path.join(downstream, "13_final")
    genes = {}
    order = []
    tx_of = defaultdict(list)
    with open(os.path.join(final, f"{prefix}.final.gff3")) as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] == "gene":
                gene_id = attr(p[8])["ID"]
                genes[gene_id] = (p[0], int(p[3]), int(p[4]), p[6])
                order.append(gene_id)
            elif p[2] == "mRNA":
                a = attr(p[8])
                tx_of[a["Parent"]].append(a["ID"])
    proteins = {k: v.rstrip("*") for k, v in load_fasta(os.path.join(final, f"{prefix}.final.proteins.fa")).items()}
    old = load_fasta(os.path.join(downstream, "02_sequences", f"{prefix}.protein.fa"))
    old_by_seq = {}
    for tid, seq in old.items():
        old_by_seq.setdefault(seq.rstrip("*."), tid)
    source = {}
    novel = []
    for gene_id in order:
        for tid in tx_of[gene_id]:
            hit = old_by_seq.get(proteins[tid])
            if hit:
                source[tid] = hit
            else:
                source[tid] = tid
                novel.append((tid, proteins[tid]))
    novel_fa = os.path.join(downstream, "05_diamond", f"{prefix}.final_novel.faa")
    write_if_changed(novel_fa, novel)
    paths = tool_paths(downstream, prefix)
    return {
        "downstream": downstream, "prefix": prefix, "final": final, "genes": genes,
        "order": order, "tx_of": tx_of, "proteins": proteins, "source": source,
        "novel": novel, "novel_fa": novel_fa, "paths": paths,
        "stale": bool(novel) and not all(ready(path, novel_fa) for path in paths.values()),
    }


def annotate_novel(states, threads, work):
    """Run the tools once on the novel proteins of every sample that needs them."""
    stale = [state for state in states if state["stale"]]
    if not stale:
        return
    if len(stale) == 1:
        run_tools(stale[0]["novel_fa"], stale[0]["paths"], threads)
        return
    os.makedirs(work, exist_ok=True)
    query = os.path.join(work, "final_novel.faa")
    with open(query, "w") as handle:
        for state in stale:
            for name, seq in state["novel"]:
                handle.write(f">{name}\n{seq}\n")
    combined = tool_paths(work, "combined")
    for path in combined.values():
        os.makedirs(os.path.dirname(path), exist_ok=True)
    run_tools(query, combined, threads)
    split_outputs(combined, {state["prefix"]: state["paths"] for state in stale})
    subprocess.run(["rm", "-rf", work], check=True)


def build_table(state):
    downstream, prefix, final = state["downstream"], state["prefix"], state["final"]
    genes, order, tx_of = state["genes"], state["order"], state["tx_of"]
    proteins, source, tools = state["proteins"], state["source"], state["paths"]
    old_ids = {v for k, v in source.items() if v != k}
    new_ids = {k for k, v in source.items() if v == k}
    down = os.path.join(downstream, "05_diamond")
    swiss = first_hits(os.path.join(down, f"{prefix}.vs_swissprot.tsv"), old_ids, stitle=True)
    refs = {}
    for label, name, _db in REFS:
        refs[label] = best_hit_table(os.path.join(down, f"{prefix}.vs_{name}.best_hit.tsv"), old_ids)
    eggnog = eggnog_table(os.path.join(downstream, "06_eggnog", f"{prefix}_eggnog.emapper.annotations"), old_ids)
    interpro = interpro_table(os.path.join(downstream, "07_interpro", f"{prefix}.interproscan.tsv"), old_ids)
    if new_ids:
        swiss.update(first_hits(tools["swissprot"], new_ids, stitle=True))
        for label, _name, _db in REFS:
            refs[label].update(first_hits(tools[label], new_ids))
        eggnog.update(eggnog_table(tools["eggnog"], new_ids))
        interpro.update(interpro_table(tools["interpro"], new_ids))

    classes = read_table(os.path.join(final, f"{prefix}.final.gene_class.tsv"))
    summary = read_table(os.path.join(downstream, "12_renamed", f"{prefix}.renamed.annotation.summary.tsv"))
    moved = {
        gene_id: genes[gene_id][:3] for gene_id in order
        if (int(summary[gene_id]["start"]), int(summary[gene_id]["end"])) != genes[gene_id][1:3]
    }
    te = te_overlap(os.path.join(downstream, "08_TE_overlap", f"{prefix}.RepeatMasker.TE.bed"), moved)

    rows = []
    for gene_id in order:
        chrom, start, end, strand = genes[gene_id]
        tids = sorted(tx_of[gene_id], key=lambda t: (-len(proteins[t]), t))
        keys = [source[t] for t in tids]

        def first(table):
            for key in keys:
                if key in table:
                    return table[key]
            return None

        eg = first(eggnog) or {}
        sw = first(swiss)
        z13 = first(refs["ZH13v3"])
        w82 = first(refs["Wm82a6"])
        ipr, pfam, ipr_go = set(), set(), set()
        for key in keys:
            if key in interpro:
                ipr |= interpro[key]["ipr"]
                pfam |= interpro[key]["pfam"]
                ipr_go |= interpro[key]["go"]
        name = clean(eg.get("Preferred_name"))
        product = clean(eg.get("Description")) or (swissprot_product(sw[2]) if sw else "") or "hypothetical protein"
        dbxref = []
        if sw:
            dbxref.append(f"UniProtKB:{sw[0].split('|')[1] if sw[0].count('|') >= 2 else sw[0]}")
        dbxref += [f"InterPro:{x}" for x in sorted(ipr)]
        dbxref += [f"Pfam:{x}" for x in sorted(pfam) if x.startswith("PF")]
        go = go_terms(eg.get("GOs", "")) | ipr_go
        info = summary[gene_id]
        if gene_id in moved:
            te_ratio, te_classes = te[gene_id]
        else:
            te_ratio, te_classes = info.get("TE_overlap_ratio", "NA"), info.get("TE_repeat_classes", "NA")
        rows.append({
            "gene_id": gene_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "gene_class": classes[gene_id]["class"],
            "transcript_count": len(tids),
            "representative_transcript": tids[0],
            "annotation_source": "修订蛋白重新注释" if any(t in new_ids for t in tids) else "沿用09结果",
            "gene_name": name or "NA",
            "product": product,
            "Dbxref": ",".join(dbxref) or "NA",
            "GO_terms": ",".join(sorted(go)) or "NA",
            "KEGG_ko": clean(eg.get("KEGG_ko")) or "NA",
            "KEGG_pathway": clean(eg.get("KEGG_Pathway")) or "NA",
            "EC": clean(eg.get("EC")) or "NA",
            "eggNOG_PFAM": clean(eg.get("PFAMs")) or "NA",
            "InterPro_IDs": ",".join(sorted(ipr)) or "NA",
            "InterPro_Pfam": ",".join(sorted(pfam)) or "NA",
            "InterPro_GO": ",".join(sorted(ipr_go)) or "NA",
            "SwissProt_hit": sw[0] if sw else "NA",
            "SwissProt_identity": sw[1] if sw else "NA",
            "SwissProt_description": swissprot_product(sw[2]) if sw else "NA",
            "ZH13v3_hit": z13[0] if z13 else "NA",
            "ZH13v3_identity": z13[1] if z13 else "NA",
            "Wm82a6_hit": w82[0] if w82 else "NA",
            "Wm82a6_identity": w82[1] if w82 else "NA",
            "TE_overlap_ratio": te_ratio,
            "TE_repeat_classes": te_classes,
            "TE_recomputed": "yes" if gene_id in moved else "no",
            "max_TPM": info.get("max_TPM", "NA"),
            "expression_supported": classes[gene_id]["expression_supported"],
        })
    header = list(rows[0])
    write_table(os.path.join(final, f"{prefix}.final.functional_annotation.tsv"), header, rows)

    n_new = sum(1 for r in rows if r["annotation_source"] == "修订蛋白重新注释")
    n_func = sum(1 for r in rows if r["product"] != "hypothetical protein" or r["InterPro_IDs"] != "NA")
    criteria_path = os.path.join(final, f"{prefix}.final.criteria.txt")
    with open(criteria_path) as handle:
        criteria = handle.read()
    if STAMP in criteria:
        criteria = criteria[:criteria.index(STAMP)].rstrip() + "\n"
    note = (
        f"{STAMP}\n"
        "final.functional_annotation.tsv 每行一个 13_final 基因，顺序与 GFF 相同；删除和合并掉的基因号不出现。\n"
        "蛋白与 09 步某条蛋白完全相同的转录本，沿用那条蛋白的 eggNOG、InterProScan、SwissProt、ZH13 v3 和 Wm82 a6 结果；"
        "10 步延长或合并后的新蛋白用同样的工具、数据库和参数重新注释，结果在 05_diamond、06_eggnog、07_interpro 的 final_novel 文件。\n"
        "基因的 eggNOG、SwissProt 和参考蛋白命中取蛋白最长、且有结果的转录本；InterPro 和 GO 取全部转录本的并集。"
        "product 依次取 eggNOG 描述、SwissProt 描述，否则为 hypothetical protein。\n"
        "坐标改变的基因按 09 步的方法重新计算 TE 重叠；max_TPM 和 expression_supported 沿用 09 步对该基因号的定量。\n"
        f"共 {len(rows)} 个基因，其中重新注释 {n_new} 个，有 eggNOG/SwissProt 描述或 InterPro 结构域的 {n_func} 个。\n"
    )
    with open(criteria_path, "w") as handle:
        handle.write(criteria + "\n" + note)
    print(f"{prefix}\tgenes={len(rows)}\treannotated={n_new}\tte_recomputed={len(moved)}\twith_function={n_func}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downstream", nargs="+", required=True)
    parser.add_argument("--prefix", nargs="+", required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--work", help="scratch directory for the combined run; removed afterwards")
    args = parser.parse_args()
    if len(args.downstream) != len(args.prefix):
        parser.error("--downstream and --prefix need the same number of values")
    states = [prepare(down, prefix) for down, prefix in zip(args.downstream, args.prefix)]
    for state in states:
        print(f"{state['prefix']}\tnovel_proteins={len(state['novel'])}\tneeds_annotation={state['stale']}", flush=True)
    work = args.work or os.path.join(os.path.dirname(os.path.abspath(args.downstream[0])), "final_function_work")
    annotate_novel(states, args.threads, work)
    for state in states:
        build_table(state)


if __name__ == "__main__":
    main()
