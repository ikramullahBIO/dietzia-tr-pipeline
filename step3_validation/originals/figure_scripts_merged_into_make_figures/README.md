# The seven figure scripts that make_figures.py merged

`../../make_figures.py` is a merge of these. Each figure's code sits inside it
as a function, unchanged. The merge was verified by regenerating all seven
figures from the merged script and confirming the output files are
byte-identical to the ones these scripts produce.

    08_make_figure2.py   -> make_figures.py --only 2   (superseded by 19_)
    19_make_figure2_v2.py-> make_figures.py --only 2   (the version used)
    13_make_figure3.py   -> make_figures.py --only 3
    10_make_figure5.py   -> make_figures.py --only 5
    17_make_figure6.py   -> make_figures.py --only 6
    16_make_figure7.py   -> make_figures.py --only 7
    09_make_figure8.py   -> make_figures.py --only 8
    14_make_figure9.py   -> make_figures.py --only 9

Figure 4 is not here: it shells out to FastANI and stays separate, as
`step3_validation.py ani`.
