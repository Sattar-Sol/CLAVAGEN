"""Minimal desktop GUI (tkinter) for clinicians: paste/load sequences, run, open the report."""

from __future__ import annotations

import os
import tempfile
import webbrowser


def main():
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, scrolledtext
    except ImportError:
        raise SystemExit("tkinter is not available. Install python3-tk or use 'clavagen run' / 'wizard'.")
    from .pipeline import Pipeline
    from .report import write_all

    root = tk.Tk()
    root.title("CLAVAGEN - clinical variant calling")
    root.geometry("980x760")

    def box(label, row, h=8):
        tk.Label(root, text=label, anchor="w", font=("Helvetica", 10, "bold")).grid(row=row, column=0, sticky="w", padx=8)
        t = scrolledtext.ScrolledText(root, height=h, width=110, font=("Courier", 9))
        t.grid(row=row + 1, column=0, columnspan=3, padx=8, sticky="we")
        def load():
            p = filedialog.askopenfilename(filetypes=[("FASTA/text", "*.fa *.fasta *.fna *.faa *.txt *.seq"), ("All", "*")])
            if p:
                t.delete("1.0", "end")
                t.insert("1.0", open(p).read())
        tk.Button(root, text="Load file...", command=load).grid(row=row, column=2, sticky="e", padx=8)
        return t

    q = box("Query sequence(s) - DNA or protein, FASTA or raw (IUPAC codes = heterozygous)", 0)
    r = box("Reference sequence(s) - cDNA/genomic DNA or protein", 2)
    fields = {}
    labels = [("Annotation JSON", "ann"), ("Gene symbol", "gene"), ("CDS range (e.g. 1-444)", "cds"),
              ("Phenotypes (comma-separated)", "ph"), ("Reviewer evidence (e.g. c.20A>T=PS3; *=PP1)", "ev"),
              ("Patient ID", "pid"), ("Output folder", "out")]
    for i, (lab, key) in enumerate(labels):
        tk.Label(root, text=lab, anchor="w").grid(row=4 + i, column=0, sticky="w", padx=8)
        e = tk.Entry(root, width=80)
        e.grid(row=4 + i, column=1, sticky="w")
        fields[key] = e
    fields["out"].insert(0, os.path.abspath("clavagen_results"))
    tk.Button(root, text="Browse...", command=lambda: (fields["ann"].delete(0, "end"),
              fields["ann"].insert(0, filedialog.askopenfilename(filetypes=[("JSON", "*.json")])))).grid(row=4, column=2)
    strict = tk.BooleanVar(value=False)
    tk.Checkbutton(root, text="Strict mode (no heuristic evidence)", variable=strict).grid(row=11, column=0, sticky="w", padx=8)
    log = scrolledtext.ScrolledText(root, height=10, width=110, font=("Courier", 9))
    log.grid(row=13, column=0, columnspan=3, padx=8, pady=6)

    def tmp(text, role):
        fd, p = tempfile.mkstemp(suffix=".fa", prefix=role)
        with os.fdopen(fd, "w") as fh:
            fh.write(text if text.lstrip().startswith(">") else f">{role}1\n{text}\n")
        return p

    def run():
        try:
            qt, rt = q.get("1.0", "end").strip(), r.get("1.0", "end").strip()
            if not qt or not rt:
                messagebox.showerror("Input", "Query and reference sequences are required.")
                return
            f = {k: e.get().strip() for k, e in fields.items()}
            pl = Pipeline([tmp(qt, "query")], [tmp(rt, "reference")], annotation=f["ann"] or None,
                          gene=f["gene"] or None, cds=f["cds"] or None,
                          evidence=[x.strip() for x in f["ev"].split(";") if x.strip()],
                          phenotypes=[x.strip() for x in f["ph"].split(",") if x.strip()],
                          strict=strict.get(), online=True, patient={"id": f["pid"]} if f["pid"] else {})
            res = pl.run()
            write_all(res, f["out"])
            log.delete("1.0", "end")
            for s in res["samples"]:
                log.insert("end", f"{s['query'].id} vs {s['reference'].id}: {s['clinical']['status']}\n")
                for v in s["variants"]:
                    log.insert("end", f"  {v['id']} {v['hgvs_c_short']} {v['hgvs_p']} {v['consequence_label']} "
                                      f"-> {v['classification']} [{', '.join(v['criteria_met'])}]\n")
            log.insert("end", f"\nReport: {os.path.join(f['out'], 'report.html')}\n")
            webbrowser.open("file://" + os.path.abspath(os.path.join(f["out"], "report.html")))
        except Exception as e:  # show any error to the user
            messagebox.showerror("Error", str(e))

    tk.Button(root, text="Run analysis", command=run, bg="#13294b", fg="white", width=20).grid(row=12, column=0, sticky="w", padx=8, pady=6)
    root.mainloop()


if __name__ == "__main__":
    main()
