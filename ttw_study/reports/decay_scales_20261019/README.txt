On-shell ttW+ production/decay scale report, 19 October 2026

Open report.pdf for the 16-page report with nine figures.
Rebuild the LaTeX document in this directory with: make
Alternatively run pdflatex on report.tex twice.

The source ZIP contains report.tex, four generated table includes, all nine
standalone PDF figures, Makefile, README.txt, and the compiled report PDF.
The directory also retains the collected arrays, numerical JSON, archived
input cards, precision inventory, and SHA256 provenance. Those data files
are not required merely to rebuild the PDF from the source ZIP.

This report concerns one specified W+ direct-lepton assignment at 13 TeV.
It compares S/Pi, varying/fixed widths, and explicit correlated/independent
production and top/antitop scale ensembles. Scale bands and Monte Carlo
errors are separate. No new integrations were launched.

Reduction scripts, from the repository root:
  python ttw_study/scripts/collect_decay_scale_report.py
  python ttw_study/scripts/render_decay_scale_report.py
The numerical inputs and prose should be reviewed together if runs change.
