# Report

This directory contains the final research report and its rendered versions:

| File | Description |
|---|---|
| `TER_report_final.md` | **Authoritative source** (Markdown). |
| `TER_report_final.tex` | LaTeX rendering (editable). |
| `TER_report_final.pdf` | PDF rendering. |
| `TER_report_final.docx` | Word rendering. |
| `figures/` | All report figures (PNG), referenced by the Markdown. |

## Rebuilding the renders

With [pandoc](https://pandoc.org) and a LaTeX engine (xelatex):

```bash
pandoc TER_report_final.md -o TER_report_final.tex --standalone
pandoc TER_report_final.md -o TER_report_final.docx
pandoc TER_report_final.md -o TER_report_final.pdf --pdf-engine=xelatex
```

Notes:
- The Markdown references figures as `figures/<name>.png` (relative to this directory).
- A math-safe variant of the Markdown is used for the LaTeX/PDF render (superscripts and
  symbols are converted to LaTeX math) so that PDF builds work with standard fonts.
- Numbers in the report are generated from the analysis outputs; if you re-run the pipeline,
  regenerate the report's tables from `results/` before rebuilding renders.
