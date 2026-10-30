// Shared setup + helpers, and Figures 1-5.
const pptxgen = require("pptxgenjs");
const fs = require("fs"), path = require("path"), os = require("os");
const OUT = path.join(os.homedir(),
  "mnt/Dietzia_analysis/results_visuallization/Dietzia_all_report_figures_editable.pptx");

const BLUE="2A78D6", ORANGE="EB6834", GREY="B8B8B4", RED="C0392B", GREEN="1F7A3D";
const INK="0B0B0B", INK2="52514E", PANEL="F4F4F1", EDGE="C9C9C4", TINT="F2F6FC";
const IDENT="E6F6E6", IDENT_E="63B563", SUBST="FBD2A8", SUBST_E="E08C3A",
      GAPC="D9D9D9", GAP_E="B0B0B0";
const F="Calibri";

const pptx = new pptxgen();
pptx.defineLayout({ name:"W", width:10, height:5.625 });
pptx.layout = "W";
pptx.author = "Ikram Ullah";
pptx.title  = "Dietzia report - every figure, editable";

// ---------- helpers ----------
const S = () => { const s = pptx.addSlide(); s.background = { color:"FFFFFF" }; return s; };
const head = (s, title, sub) => {
  s.addText(title, { x:0.42, y:0.16, w:9.16, h:0.36, fontFace:F, fontSize:16, bold:true,
                     color:INK, align:"left", valign:"middle" });
  if (sub) s.addText(sub, { x:0.42, y:0.52, w:9.16, h:0.30, fontFace:F, fontSize:10,
                     italic:true, color:INK2, align:"left", valign:"middle" });
};
const foot = (s, t) => s.addText(t, { x:0.42, y:5.26, w:9.16, h:0.26, fontFace:F,
  fontSize:8, italic:true, color:INK2, align:"left", valign:"middle" });
const rect = (s,x,y,w,h,fill,line,dash) => s.addShape(pptx.ShapeType.rect, {
  x,y,w,h, fill:{color:fill}, line: line ? {color:line, width:0.75, dashType: dash||"solid"} : {type:"none"} });
const rrect = (s,x,y,w,h,fill,line) => s.addShape(pptx.ShapeType.roundRect, {
  x,y,w,h, rectRadius:0.06, fill:{color:fill}, line:{color:line, width:1} });
const txt = (s,x,y,w,h,t,o={}) => s.addText(t, Object.assign(
  { x,y,w,h, fontFace:F, fontSize:10, color:INK, align:"center", valign:"middle", margin:0 }, o));
const arrow = (s,x,y,w,h,rot) => s.addShape(pptx.ShapeType.rightArrow, {
  x,y,w,h, fill:{color:INK2}, line:{type:"none"}, rotate: rot||0 });

// ============================================================ FIGURE 1
{
const s = S();
head(s, "Figure 1  —  MLVA panel design and validation workflow",
        "Redrawn as editable shapes. The published figure is a BioRender illustration; this version keeps the same six stages and wording.");
const BOXW = 2.55, BOXH = 0.78;
const stages = [
  {x:0.45, y:1.05, t:"Reference genome\n+ comparison assemblies", d:"GCA_045016705.1 and the\n83-genome genomospecies-A panel"},
  {x:3.72, y:1.05, t:"TRF scanning &\nmotif identification", d:"Tandem Repeats Finder,\nstringent 2-7-7, score ≥ 50"},
  {x:6.99, y:1.05, t:"Candidate TR loci", d:"772 loci on the reference,\nprimer-flanked and de-duplicated"},
  {x:6.99, y:2.62, t:"In-silico PCR &\nlocus testing", d:"MLVA_finder2_ipcress across\nthe 83-genome panel"},
  {x:3.72, y:2.62, t:"Filtering & selection", d:"≥ 5 binned alleles, ≤ 10% missing,\n3% size tolerance"},
  {x:0.45, y:2.62, t:"Ranking by\ndiscriminatory power", d:"Hunter–Gaston index,\ncumulative, locus by locus"},
];
stages.forEach(b => {
  rrect(s, b.x, b.y, BOXW, BOXH, TINT, BLUE);
  txt(s, b.x, b.y+0.04, BOXW, 0.44, b.t, {fontSize:10.5, bold:true, color:INK});
  txt(s, b.x, b.y+0.44, BOXW, 0.32, b.d, {fontSize:7.8, color:INK2});
});
// arrows: 1->2->3 (right), 3->4 (down), 4->5->6 (left)
arrow(s, 3.06, 1.36, 0.58, 0.16);
arrow(s, 6.33, 1.36, 0.58, 0.16);
arrow(s, 8.12, 1.90, 0.58, 0.16, 90);
arrow(s, 6.33, 2.93, 0.58, 0.16, 180);
arrow(s, 3.06, 2.93, 0.58, 0.16, 180);
arrow(s, 1.44, 3.48, 0.58, 0.16, 90);
// final panel
rrect(s, 2.30, 3.95, 5.40, 1.10, "EAF2E8", GREEN);
txt(s, 2.30, 3.99, 5.40, 0.30, "Final in-silico MLVA panel for wet-lab validation",
    {fontSize:11, bold:true, color:GREEN});
txt(s, 2.30, 4.28, 5.40, 0.28, "9 loci  ·  repeat units 3–18 bp  ·  83 genomospecies-A genomes",
    {fontSize:9.5, color:INK2});
txt(s, 2.30, 4.56, 5.40, 0.32, "Combined HGDI = 0.9880   (61 distinct profiles)",
    {fontSize:11.5, bold:true, color:INK});
foot(s, "Numbers are those reported in Section 3.1 and Table 1. The published figure quotes HGDI 0.989, the 3%-binned value; 0.9880 is the raw-allele value used in the text.");
s.addNotes("Figure 1 redrawn from the report's stage sequence. The original is a BioRender illustration "
 + "(results_visuallization/Report figure 1.png) and cannot be reproduced shape-for-shape; this version is "
 + "structurally equivalent and fully editable.\n"
 + "Note the two HGDI values: 0.9880 raw alleles (Section 3.1) vs 0.9903 binned at 3% by the independent "
 + "selection script. The published figure's 0.989 is the binned one.");
}

// ============================================================ FIGURE 2
{
const B = require(path.join(os.homedir(),
  "mnt/Dietzia_analysis/Project_scripts/Step3_validation_and_figures/fig2_blocks.js"));
const KIND = { cds:[BLUE,"1B4F9C","solid"], rep:[ORANGE,"B4521F","solid"],
               out:[GREY,"6E6E6A","dash"],  is:[RED,"8E2A1E","solid"],
               dup:[ORANGE,"FFFFFF","dash"] };
const s = S();
head(s, "Figure 2  —  CDS length-difference classification workflow",
        "Five criteria, applied in order; the first match wins. Each pair of bars is one orthologous CDS drawn to scale.");
const X0 = 1.35, SC = 4.20/89, ROWH = 0.17, NX = 5.72, NW = 3.95;
B.forEach((b, i) => {
  const yT = 0.98 + i*0.84;
  txt(s, 0.42, yT, 0.28, 0.22, b.n, {fontSize:12, bold:true, color:ORANGE, align:"left"});
  txt(s, 0.72, yT, 5.00, 0.22, b.title, {fontSize:10, bold:true, color:INK, align:"left"});
  txt(s, 0.72, yT+0.20, 5.00, 0.18, b.rule, {fontSize:7.5, color:INK2, align:"left"});
  [["reference", b.ref, yT+0.40], ["neighbour", b.ngb, yT+0.58]].forEach(([lbl, segs, y]) => {
    txt(s, X0-1.02, y, 0.97, ROWH, lbl, {fontSize:7.5, color:INK2, align:"right"});
    segs.forEach(([st, ln, kind]) => {
      const [fc, ec, ds] = KIND[kind];
      rect(s, X0 + st*SC, y, ln*SC, ROWH, fc, ec, ds);
    });
  });
  txt(s, NX, yT+0.18, NW, 0.20, b.note, {fontSize:7.5, italic:true, color:INK2, align:"left"});
  txt(s, NX, yT+0.38, NW, 0.36, b.stat, {fontSize:7, color:INK2, align:"left"});
});
foot(s, "Counts and percentages are those of Table 2. Blue, coding sequence shared by both genomes. Orange, repeat unit or duplicated block. Grey dashed, DNA present in both genomes but outside the called CDS. Dark red, transposase or integrase.");
}
// ============================================================ FIGURE 3
{
const s = S();
head(s, "Figure 3  —  Cumulative discriminatory power of the 9-locus MLVA panel",
        "83 genomospecies-A genomes, raw allele sizes, failure to amplify scored as its own state. 61 distinct profiles.");
const HG = [0.0943,0.6582,0.9192,0.9600,0.9733,0.9853,0.9874,0.9877,0.9880];
s.addChart(pptx.ChartType.line, [
  { name:"Cumulative HGDI", labels: HG.map((_,i)=>String(i+1)), values: HG },
  { name:"0.95 threshold",  labels: HG.map((_,i)=>String(i+1)), values: HG.map(()=>0.95) },
], {
  x:0.55, y:0.95, w:8.9, h:4.10,
  chartColors:[BLUE, ORANGE], lineDataSymbol:["circle","none"],
  lineSize:2, lineDataSymbolSize:7,
  lineDash:["solid","dash"],
  showLegend:true, legendPos:"b", legendFontSize:9, legendFontFace:F,
  showValue:true, dataLabelFontSize:8, dataLabelFontFace:F, dataLabelColor:INK2,
  dataLabelFormatCode:"0.000", dataLabelPosition:"t",
  catAxisTitle:"Number of loci included (cumulative, in panel-selection order)",
  showCatAxisTitle:true, catAxisTitleFontSize:9, catAxisTitleFontFace:F,
  valAxisTitle:"Hunter–Gaston Discriminatory Index",
  showValAxisTitle:true, valAxisTitleFontSize:9, valAxisTitleFontFace:F,
  valAxisMinVal:0, valAxisMaxVal:1.05, valAxisMajorUnit:0.2,
  catAxisLabelFontSize:9, valAxisLabelFontSize:9,
  catAxisLabelFontFace:F, valAxisLabelFontFace:F,
});
foot(s, "Native PowerPoint chart: right-click → Edit Data to change any value. The panel crosses 0.95 at the fourth locus and is flat from the sixth.");
s.addNotes("Values recomputed from MLVA_analysis_genomes.csv with the Hunter–Gaston formula, "
 + "excluding GCA_024732245.1 and the two IBK-p2-16D assemblies (other genomospecies): "
 + "0.0943, 0.6582, 0.9192, 0.9600, 0.9733, 0.9853, 0.9874, 0.9877, 0.9880 for N = 83, 61 profiles.");
}

// ============================================================ FIGURE 4  (2 slides)
{
const meta = fs.readFileSync(path.join(os.homedir(),
  "mnt/Dietzia_analysis/Project_scripts/Step3_validation_and_figures/ani_metadata.csv"),"utf8")
  .trim().split("\n").slice(1).map(l=>{
    const c=[]; let cur="", q=false;
    for (const ch of l){ if(ch==='"'){q=!q;} else if(ch===','&&!q){c.push(cur);cur="";} else cur+=ch; }
    c.push(cur); return c;
  });
const GROUPS = fs.readFileSync(path.join(os.homedir(),
  "mnt/Dietzia_analysis/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/03_three_genome_comparison/figures/fig4_ANI_clustered.groups.tsv"),"utf8")
  .trim().split("\n").slice(1).map(l=>l.split("\t"));
const MX = fs.readFileSync(path.join(os.homedir(),
  "mnt/Dietzia_analysis/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/03_three_genome_comparison/figures/fig4_ANI_clustered.matrix.csv"),"utf8")
  .trim().split("\n").map(l=>l.split(","));
const order = MX[0].slice(1);                       // clustered order
const mrow = Object.fromEntries(meta.map(r=>[r[0], r]));
const grow = Object.fromEntries(GROUPS.map(r=>[r[0], r]));
const REF = "GCA_045016705.1";

// ---- slide 4a: metadata table
{
const s = S();
head(s, "Figure 4 (1 of 2)  —  Pairwise ANI, 12-genome discovery panel: genome metadata",
        "Rows ordered by UPGMA clustering of (100 − ANI); the same row numbers are used in the matrix on the next slide.");
const rows = [[ "#", "Accession", "Strain", "Species", "Country", "Year", "Clade", "ANI to ref.", "Cluster" ]
  .map(t=>({ text:t, options:{ bold:true, color:"FFFFFF", fill:{color:INK2}, fontSize:8.5 } }))];
order.forEach((acc, i) => {
  const m = mrow[acc] || [], g = grow[acc] || [];
  const isRef = acc === REF;
  const fill = isRef ? { color:TINT } : (i % 2 ? { color:"F7F7F5" } : { color:"FFFFFF" });
  rows.push([ String(i+1), acc, m[1]||"", m[2]||"", m[3]||"", m[4]||"", m[8]||"—",
              isRef ? "— (reference)" : (g[2]||""), g[1]||"" ]
    .map((t,ci)=>({ text:t, options:{ fontSize:8, bold:isRef, italic: ci===3,
                    fill, color: isRef ? BLUE : INK } })));
});
s.addTable(rows, { x:0.35, y:0.95, w:9.30, colW:[0.32,1.32,1.15,1.72,1.02,0.48,0.62,1.20,0.47],
  border:{ type:"solid", color:EDGE, pt:0.5 }, fontFace:F, valign:"middle", rowH:0.26, autoPage:false });
foot(s, "Metadata from the NCBI assembly report. The clade column gives the assignment of Dos Santos et al. (2023) where one exists; an em dash means unassigned.");
}

// ---- slide 4b: ANI matrix
{
const s = S();
head(s, "Figure 4 (2 of 2)  —  Pairwise ANI, 12-genome discovery panel: matrix",
        "FastANI, all 66 pairs. Shading marks the ≥ 95% ANI species clusters; bold entries are pairs inside one cluster.");
const rows = [[{ text:"", options:{ fill:{color:INK2} } }].concat(
  order.map((_,i)=>({ text:String(i+1), options:{ bold:true, color:"FFFFFF", fill:{color:INK2}, fontSize:8, align:"center" } })))];
order.forEach((acc, i) => {
  const cl = (grow[acc]||[])[1];
  const cells = [{ text:`${i+1}  ${(mrow[acc]||[])[1]||acc}`,
    options:{ bold:true, fontSize:7.5, align:"left",
              fill:{ color: acc===REF ? TINT : "F7F7F5" }, color: acc===REF ? BLUE : INK } }];
  order.forEach((acc2, j) => {
    const v = parseFloat(MX[i+1][j+1]);
    const same = cl && cl === (grow[acc2]||[])[1];
    const diag = i === j;
    cells.push({ text: diag ? "—" : v.toFixed(1),
      options:{ fontSize:7.5, align:"center", bold: same && !diag,
                fill:{ color: diag ? "E8E8E4" : (same ? "DCEAF8" : "FFFFFF") },
                color: same && !diag ? "1B4F9C" : INK2 } });
  });
  rows.push(cells);
});
s.addTable(rows, { x:0.30, y:0.95, w:9.40, colW:[1.96].concat(Array(12).fill(0.62)),
  border:{ type:"solid", color:EDGE, pt:0.5 }, fontFace:F, valign:"middle", rowH:0.27, autoPage:false });
foot(s, "Values are percent identity to one decimal; the full three-decimal matrix is in ANI_analysis/results/ani_matrix_FastANI_real.csv. Row 6 (BIM B-2035) is the reference.");
s.addNotes("Cluster 5 is the genomospecies-A group: BIM B-2035 (reference), NBRC 15801, MNB45, 313 and "
 + "IITR165. GCA_055568915.1 / NBRC 15801 is the Dietzia maris type-strain assembly and sits at 96.884% "
 + "to the reference, inside the 95% boundary — see Section 3.2.");
}
}

// ============================================================ FIGURE 5
{
const s = S();
head(s, "Figure 5  —  CDS length-difference attribution, percentage of pairs per category",
        "Categories are mutually exclusive and applied in the order of Section 2.4, so each pair is counted once.");
// counts
const ROWS_A = [
  ["Identical length",        2383, 2264],
  ["Tandem repeat",            204,  202],
  ["Start shift",               71,   68],
  ["Stop shift",                23,   26],
  ["IS element",                10,    8],
  ["Large indel (unresolved)",  39,   30],
  ["Other internal indel",     411,  411],
];
const ROWS_B = ROWS_A.slice(1);
const pc = (n, d) => Math.round(1000*n/d)/10;        // one rounding step only

// legend, centred at the top
const lgx = 3.05;
rect(s, lgx, 0.86, 0.17, 0.13, BLUE, "1B4F9C");
txt(s, lgx+0.21, 0.83, 1.60, 0.20, "N1  GCA_051994785.1",
    {fontSize:9, align:"left", color:INK});
rect(s, lgx+1.90, 0.86, 0.17, 0.13, ORANGE, "B4521F");
txt(s, lgx+2.11, 0.83, 1.60, 0.20, "N2  GCA_045689975.1",
    {fontSize:9, align:"left", color:INK});

const panel = (rows, den1, den2, title, LX, BX, BW, vmax, step) => {
  txt(s, LX, 1.12, BX+BW-LX, 0.22, title, {fontSize:9.5, bold:true, align:"left", color:INK});
  const TOP = 1.44, ROWH = 3.28/rows.length, BARH = Math.min(0.15, ROWH*0.30);
  const BOT = TOP + rows.length*ROWH;
  // gridlines + axis ticks
  for (let v = 0; v <= vmax + 1e-9; v += step) {
    const x = BX + (v/vmax)*BW;
    s.addShape(pptx.ShapeType.line, { x, y:TOP-0.06, w:0, h:(BOT-TOP)+0.06,
      line:{ color: v===0 ? "9A9A96" : "E4E4E0", width: v===0 ? 0.9 : 0.6 } });
    txt(s, x-0.22, BOT+0.02, 0.44, 0.18, String(v), {fontSize:7.5, color:INK2});
  }
  rows.forEach(([name, n1, n2], i) => {
    const y = TOP + i*ROWH + (ROWH - 2*BARH - 0.03)/2;
    txt(s, LX, y, BX-LX-0.06, 2*BARH+0.03, name,
        {fontSize:8, color:INK2, align:"right"});
    [[pc(n1, den1), BLUE, "1B4F9C", y],
     [pc(n2, den2), ORANGE, "B4521F", y+BARH+0.03]].forEach(([v, fc, ec, yy]) => {
      const w = Math.max(0.006, (v/vmax)*BW);
      rect(s, BX, yy, w, BARH, fc, ec);
      txt(s, BX+w+0.04, yy-0.01, 0.62, BARH+0.02, v.toFixed(1) + "%",
          {fontSize:7.5, color:INK2, align:"left"});
    });
  });
  txt(s, BX, BOT+0.20, BW, 0.20, "% of pairs", {fontSize:8, color:INK2});
};
panel(ROWS_A, 3141, 3009, "(a)  share of all orthologous pairs",
      0.32, 1.92, 2.52, 80, 20);
panel(ROWS_B,  758,  745, "(b)  share of the pairs that differ in length",
      5.12, 6.72, 2.52, 60, 10);
foot(s, "Denominators: 3,141 (N1) and 3,009 (N2) orthologous pairs in panel (a); 758 and 745 differing pairs in panel (b). Underlying counts in Table 2.");
s.addNotes("Counts: identical 2,383 / 2,264; tandem repeat 204 / 202; start shift 71 / 68; stop shift "
 + "23 / 26; IS element 10 / 8; large indel 39 / 30; other internal indel 411 / 411. "
 + "Denominators 3,141 and 3,009 (all pairs), 758 and 745 (differing pairs).\n"
 + "Drawn as shapes rather than a native chart: each bar is a rectangle and each value a text box, so "
 + "nothing depends on PowerPoint re-reading an embedded workbook.");
}
// ============================================================ FIGURE 6
{
const s = S();
head(s, "Figure 6  —  The one large-indel pair that Swiss-Prot resolves",
        "CP155437.1_3055 (reference) against its N2 ortholog CP143297.1_2802 (GCA_045689975.1), a difference of +648 bp.");
const X0 = 1.75, SC = 6.9/313, H = 0.44;
[["reference", 2.05, 97, true], ["N2 ortholog", 3.00, 313, false]].forEach(([lbl, y, aa, frag]) => {
  txt(s, 0.35, y, 1.30, H, lbl, {fontSize:10, color:INK2, align:"right"});
  if (frag) {
    rect(s, X0, y, (313-97)*SC, H, "FFFFFF", GAP_E, "dash");
    rect(s, X0 + (313-97)*SC, y, 97*SC, H, IDENT, IDENT_E);
    txt(s, X0, y, (313-97)*SC, H, "absent from the reference", {fontSize:8, color:INK2, italic:true});
  } else {
    rect(s, X0, y, (313-97)*SC, H, SUBST, SUBST_E);
    rect(s, X0 + (313-97)*SC, y, 97*SC, H, IDENT, IDENT_E);
  }
  txt(s, X0 + 313*SC + 0.08, y, 1.55, H, aa + " aa", {fontSize:9.5, bold:true, align:"left", color:INK});
});
txt(s, X0, 1.72, 6.9, 0.26, "←   216 residues present only in the neighbour   →          ←  shared C-terminal third  →",
    {fontSize:8, color:INK2, align:"left"});
rrect(s, 1.75, 3.75, 7.0, 1.20, PANEL, EDGE);
txt(s, 1.90, 3.84, 6.7, 0.30, "Swiss-Prot search", {fontSize:10.5, bold:true, align:"left", color:INK});
txt(s, 1.90, 4.12, 6.7, 0.32,
  "N2 copy, 313 aa, annotated SDR family oxidoreductase   —   e = 2 × 10⁻²⁴",
  {fontSize:9.5, align:"left", color:INK2});
txt(s, 1.90, 4.42, 6.7, 0.32,
  "Reference copy, 97 aa, annotated hypothetical protein   —   e = 0.002",
  {fontSize:9.5, align:"left", color:INK2});
txt(s, 1.90, 4.70, 6.7, 0.22,
  "The reference scores weakly because the fragment is short, not because it has diverged.",
  {fontSize:8.5, italic:true, align:"left", color:INK2});
// key
txt(s, 0.42, 5.02, 1.2, 0.20, "Key:", {fontSize:8, align:"left", color:INK2});
rect(s, 0.85, 5.05, 0.18, 0.14, IDENT, IDENT_E);
txt(s, 1.08, 5.02, 1.5, 0.20, "in both copies", {fontSize:8, align:"left", color:INK2});
rect(s, 2.55, 5.05, 0.18, 0.14, SUBST, SUBST_E);
txt(s, 2.78, 5.02, 2.0, 0.20, "neighbour only", {fontSize:8, align:"left", color:INK2});
foot(s, "Bars drawn to scale in amino acids. Section 3.3.");
}

// ============================================================ FIGURE 7
{
const s = S();
head(s, "Figure 7  —  Tandem-repeat genomic context",
        "Intragenic fraction of the relaxed census, compared with raw genomic coding density. Descriptive; no statistical test attached.");
txt(s, 0.45, 0.88, 4.40, 0.24, "(a)  relaxed TRF census, 9,914 loci",
    {fontSize:10, bold:true, align:"left", color:INK});
s.addChart(pptx.ChartType.bar, [
  { name:"Relaxed census (2-3-5)", labels:["Intragenic","Intergenic"], values:[94.07, 5.93] }
], { x:0.45, y:1.12, w:4.40, h:3.90, barDir:"col", chartColors:[BLUE],
     showLegend:false, showValue:true, dataLabelFontSize:9, dataLabelFontFace:F,
     dataLabelColor:INK2, dataLabelFormatCode:'0.0"%"',
     valAxisTitle:"% of tandem-repeat loci", showValAxisTitle:true,
     valAxisTitleFontSize:8.5, valAxisTitleFontFace:F, valAxisMaxVal:100,
     catAxisLabelFontSize:9, valAxisLabelFontSize:8.5,
     catAxisLabelFontFace:F, valAxisLabelFontFace:F, barGapWidthPct:110 });
txt(s, 5.15, 0.88, 4.40, 0.24, "(b)  raw genomic coding density",
    {fontSize:10, bold:true, align:"left", color:INK});
s.addChart(pptx.ChartType.bar, [
  { name:"Prodigal coding density", labels:["Reference","N1","N2"], values:[90.52, 90.39, 91.12] }
], { x:5.15, y:1.12, w:4.40, h:3.90, barDir:"col", chartColors:[ORANGE],
     showLegend:false, showValue:true, dataLabelFontSize:9, dataLabelFontFace:F,
     dataLabelColor:INK2, dataLabelFormatCode:'0.0"%"',
     valAxisTitle:"% of assembly covered by CDS", showValAxisTitle:true,
     valAxisTitleFontSize:8.5, valAxisTitleFontFace:F,
     valAxisMinVal:80, valAxisMaxVal:100,
     catAxisLabelFontSize:9, valAxisLabelFontSize:8.5,
     catAxisLabelFontFace:F, valAxisLabelFontFace:F, barGapWidthPct:110 });
foot(s, "Relaxed census: 9,326 of 9,914 loci intragenic (94.07%), 588 intergenic. Coding density from the same Prodigal calls used throughout. In the stringent 772-locus set the intragenic fraction is 76.3% — see Section 3.5.");
s.addNotes("The comparison asks whether tandem repeats sit inside genes more often than coding sequence "
 + "alone would predict. Relaxed set: 94.1% vs 90.4% coding density, modestly above. Stringent set: 76.3%, "
 + "well below. The two sets run in opposite directions; Section 3.5 and 4.6 discuss this.");
}

// ============================================================ FIGURE 8
{
const D = JSON.parse(fs.readFileSync(path.join(os.homedir(),"fig8_data.json"),"utf8"));
const s = S();
head(s, "Figure 8  —  Two mechanisms behind CDS length differences (reference vs strain 313, GCA_024125715.1)",
        "Needleman–Wunsch, BLOSUM62. Reference on the upper row of each block, neighbour below. Only the windows around each gap block are shown.");
const cellH = 0.20, LEFT = 1.28, RIGHT = 9.62;
let y = 0.98;
const drawWindow = (win, label) => {
  const n = win.ref.length;
  const cw = Math.min(0.185, (RIGHT-LEFT)/n);
  txt(s, 0.42, y, 8.9, 0.20, label, {fontSize:8.5, bold:true, align:"left", color:INK2});
  y += 0.23;
  [["reference", win.ref], ["strain 313", win.ngb]].forEach(([lbl, seq], ri) => {
    txt(s, 0.30, y, 0.92, cellH, lbl, {fontSize:7.5, color:INK2, align:"right"});
    for (let i = 0; i < n; i++) {
      const a = win.ref[i], b = win.ngb[i];
      let fc = IDENT, ec = IDENT_E;
      if (a === "-" || b === "-") { fc = GAPC; ec = GAP_E; }
      else if (a !== b) { fc = SUBST; ec = SUBST_E; }
      const ch = seq[i];
      rect(s, LEFT + i*cw, y, cw, cellH, fc, ec);
      txt(s, LEFT + i*cw, y, cw, cellH, ch === "-" ? "" : ch,
          {fontSize:7, color:INK, bold:false});
    }
    y += cellH;
  });
  y += 0.16;
};
txt(s, 0.42, y, 9.0, 0.22,
  "(A)  N-terminal truncation via alternative start-codon usage — CP155437.1_2816 (427 aa) vs CP099712.1_2753 (402 aa), −75 bp",
  {fontSize:9.5, bold:true, align:"left", color:INK}); y += 0.25;
drawWindow(D.A.windows[0], "alignment columns 1–45  ·  25-residue (75 bp) extension present only in the reference");
txt(s, 0.42, y, 9.0, 0.20,
  "Outside this window the two proteins are 99.7% identical (381 of 382 aligned residues). The DNA is present in both genomes — see Figure 9A.",
  {fontSize:8, italic:true, align:"left", color:INK2}); y += 0.32;
txt(s, 0.42, y, 9.0, 0.22,
  "(B)  In-frame duplication — CP155437.1_367 (373 aa) vs CP099712.1_295 (394 aa), +63 bp",
  {fontSize:9.5, bold:true, align:"left", color:INK}); y += 0.25;
drawWindow(D.B.windows[0], "alignment columns 1–41  ·  duplication 1, 15 residues (45 bp), near the N-terminus");
drawWindow(D.B.windows[1], "alignment columns 201–246  ·  duplication 2, 6 residues (18 bp), central region");
// key
txt(s, 0.42, 5.28, 0.5, 0.18, "Key:", {fontSize:7.5, align:"left", color:INK2});
let kx = 0.80;
[[IDENT, IDENT_E, "identical"], [SUBST, SUBST_E, "substitution"], [GAPC, GAP_E, "gap"]].forEach(([f,e,l])=>{
  rect(s, kx, 5.30, 0.15, 0.13, f, e); txt(s, kx+0.18, 5.27, 1.1, 0.18, l,
    {fontSize:7.5, align:"left", color:INK2}); kx += 1.25;
});
txt(s, kx, 5.27, 4.2, 0.18,
  "Outside the shown windows, panel B is 98.7% identical (303 of 307).",
  {fontSize:7.5, italic:true, align:"left", color:INK2});
s.addNotes("Windows and identities recomputed with align_utils.nw_align, the same function the pipeline "
 + "uses: panel A window 1–45, identity outside 381/382 = 99.7%; panel B windows 1–41 and 201–246, "
 + "identity outside 303/307 = 98.7%. These match fig7_worked_examples.stats.tsv exactly.\n"
 + "Panel B caveat: the pair is assigned to the tandem-repeat category because two repeat loci inside the "
 + "CDS differ in copy number, but neither sits at the 45 bp duplication, so it fails the strict "
 + "position-and-magnitude test. It is an illustration of in-frame duplication, not a validated tandem repeat.");
}
// ============================================================ FIGURE 9
{
const D = JSON.parse(fs.readFileSync(path.join(os.homedir(),"fig9_data.json"),"utf8"));
const s = S();
head(s, "Figure 9  —  The two loci from Figure 8 examined further",
        "(A) the start-codon locus at nucleotide level; (B) copies of the 45 bp repeat unit across the five complete genomospecies-A genomes.");
// ---- panel A
txt(s, 0.42, 0.88, 9.0, 0.22,
  "(A)  Start-codon locus, first 26 codons — both genomes carry the same GTG at codon 1 and the same ATG at codon 26",
  {fontSize:9.5, bold:true, align:"left", color:INK});
const NC = 26, cw = 0.325, LEFT = 0.98;
let y = 1.18;
txt(s, 0.10, y, 0.80, 0.18, "codon", {fontSize:7, color:INK2, align:"right"});
for (let i=0;i<NC;i++) txt(s, LEFT+i*cw, y, cw, 0.18, String(i+1), {fontSize:6.5, color:INK2});
y += 0.20;
[["reference", D.A_ref, D.aa_ref], ["strain 313", D.A_ngb, D.aa_ngb]].forEach(([lbl, nt, aa]) => {
  txt(s, 0.10, y, 0.80, 0.22, lbl, {fontSize:7.5, color:INK2, align:"right"});
  for (let i=0;i<NC;i++){
    const cod = nt.slice(i*3, i*3+3).toUpperCase();
    const diff = D.codon_diff_idx.indexOf(i) >= 0;
    rect(s, LEFT+i*cw, y, cw, 0.22, diff?SUBST:IDENT, diff?SUBST_E:IDENT_E);
    txt(s, LEFT+i*cw, y, cw, 0.22, cod, {fontSize:6.5, color:INK, bold:diff});
  }
  y += 0.22;
  txt(s, 0.10, y, 0.80, 0.18, "residue", {fontSize:7, color:INK2, align:"right"});
  for (let i=0;i<NC;i++) txt(s, LEFT+i*cw, y, cw, 0.18, aa[i]||"", {fontSize:7, color:INK2});
  y += 0.24;
});
rect(s, LEFT, 1.38, cw, 0.70, "FFFFFF", BLUE);         // codon 1 outline
rect(s, LEFT+25*cw, 1.38, cw, 0.70, "FFFFFF", GREEN);  // codon 26 outline
txt(s, 0.42, y, 9.0, 0.20,
  "The only difference across the 78 bp is codon 7, GCC → GCT: synonymous, both encode alanine. The 75 bp length difference is a gene-calling decision, not a difference in the DNA.",
  {fontSize:8, italic:true, align:"left", color:INK2});
// ---- panel B
txt(s, 0.42, 2.72, 9.0, 0.22,
  "(B)  Copies of the 45 bp repeat unit of the Figure 8B locus, counted directly in the DNA",
  {fontSize:9.5, bold:true, align:"left", color:INK});
const GEN = [["BIM B-2035  (reference)",1],["313",2],["IITR165",1],["MNB45",1],["NBRC 15801",1]];
const UX = 3.05, UW = 0.95, UH = 0.30;
let gy = 3.02;
GEN.forEach(([name, n]) => {
  txt(s, 0.42, gy, 2.50, UH, name, {fontSize:9, color: n===2?INK:INK2, bold:n===2, align:"right"});
  for (let k=0;k<n;k++){
    rect(s, UX + k*(UW+0.05), gy, UW, UH, ORANGE, "B4521F");
    txt(s, UX + k*(UW+0.05), gy, UW, UH, "45 bp", {fontSize:8, color:"FFFFFF", bold:true});
  }
  txt(s, UX + 2*(UW+0.05) + 0.15, gy, 3.6, UH,
      n===2 ? "2 copies — the minority allele" : "1 copy",
      {fontSize:8.5, color: n===2?INK:INK2, italic:n!==2, align:"left"});
  gy += UH + 0.09;
});
txt(s, 0.42, 5.02, 9.1, 0.22,
  "Four genomes carry one copy; only strain 313 carries two. With a single copy in the reference, TRF cannot call an array here, which is why the pair fails the strict position-and-magnitude test.",
  {fontSize:8, italic:true, align:"left", color:INK2});
foot(s, "Repeat unit 5'-GGCGCGGGCCGCCCGGGCGGTGCGGGCGACACGCCCGGGTCGGAC-3'.");
s.addNotes("Panel A nucleotides extracted from the assemblies with the same code as make_figures.figure_9: "
 + "the only codon difference across the first 26 codons is codon 7 (GCC vs GCT), synonymous.\n"
 + "Panel B: an exact 45-bp match returns 1 copy for the reference and NBRC 15801, 2 for strain 313, and 0 "
 + "for IITR165 and MNB45 — those two carry a single copy that differs by one residue (T→S), so the "
 + "exact search misses it. The figure follows the report and shows them as one copy each.");
}

// ============================================================ FIGURE 10
{
const s = S();
head(s, "Figure 10  —  Three further mechanisms (reference vs N1, GCA_051994785.1)",
        "Coding sequences drawn to scale. Green, sequence present in both genomes; orange, the length difference itself; dashed outline, sequence absent from that genome.");
const X0 = 2.25, MAXAA = 951, SC = 6.55/MAXAA;
const panels = [
  { tag:"A", y:1.00, title:"Tandem repeat  —  DNA polymerase III subunit gamma/tau",
    ref:{k:"CP155437.1_425", aa:918}, ngb:{k:"CP163397.1_388", aa:951}, d:"+99 bp",
    note:"Array of a 6 bp unit, 5'-CCGGCT-3', encoding Pro–Ala in a low-complexity linker. Reference 47 copies, neighbour 63. Sixteen extra copies predict +96 bp against +99 bp observed, so the pair passes the strict test.",
    extraAtEnd:true },
  { tag:"B", y:2.42, title:"IS element  —  IS21 family transposase",
    ref:{k:"CP155437.1_3243", aa:524}, ngb:{k:"CP163397.1_2742", aa:266}, d:"−774 bp",
    note:"The reference gene is the intact transposase; the neighbour keeps both ends but has lost an internal block of 261 residues.",
    internal:[130, 261] },
  { tag:"C", y:3.84, title:"Unresolved internal indel  —  SRPBCC family protein",
    ref:{k:"CP155437.1_2935", aa:327}, ngb:{k:"CP163397.1_2994", aa:287}, d:"−120 bp",
    note:"A single contiguous block of 40 residues, in frame, coinciding with no tandem-repeat locus and no IS annotation — the residual category.",
    internal:[150, 40] },
];
panels.forEach(p => {
  txt(s, 0.28, p.y, 0.40, 0.22, "(" + p.tag + ")", {fontSize:10, bold:true, color:ORANGE, align:"left"});
  txt(s, 0.62, p.y, 6.2, 0.22, p.title, {fontSize:9.5, bold:true, color:INK, align:"left"});
  txt(s, 6.85, p.y, 2.75, 0.22, p.d, {fontSize:10, bold:true, color:ORANGE, align:"right"});
  const rows = [["reference", p.ref, p.y+0.26], ["N1", p.ngb, p.y+0.50]];
  rows.forEach(([lbl, g, yy]) => {
    txt(s, 0.30, yy, 1.85, 0.20, lbl + "  " + g.k, {fontSize:7.5, color:INK2, align:"right"});
    const other = (g === p.ref) ? p.ngb : p.ref;
    const shorter = Math.min(p.ref.aa, p.ngb.aa), longer = Math.max(p.ref.aa, p.ngb.aa);
    const isLong = g.aa === longer;
    if (p.extraAtEnd) {
      rect(s, X0, yy, shorter*SC, 0.20, IDENT, IDENT_E);
      if (isLong) rect(s, X0 + shorter*SC, yy, (longer-shorter)*SC, 0.20, SUBST, SUBST_E);
      else        rect(s, X0 + shorter*SC, yy, (longer-shorter)*SC, 0.20, "FFFFFF", GAP_E, "dash");
    } else {
      const [st, len] = p.internal;
      rect(s, X0, yy, st*SC, 0.20, IDENT, IDENT_E);
      if (isLong) rect(s, X0 + st*SC, yy, len*SC, 0.20, SUBST, SUBST_E);
      else        rect(s, X0 + st*SC, yy, len*SC, 0.20, "FFFFFF", GAP_E, "dash");
      rect(s, X0 + (st+len)*SC, yy, (longer-st-len)*SC, 0.20, IDENT, IDENT_E);
    }
    txt(s, X0 + longer*SC + 0.06, yy, 0.95, 0.20, g.aa + " aa",
        {fontSize:8, bold:true, color:INK, align:"left"});
  });
  txt(s, 0.58, p.y+0.74, 9.0, 0.30, p.note, {fontSize:7.5, italic:true, color:INK2, align:"left"});
});
foot(s, "All three genomes called ab initio with Prodigal, on both sides of every pair, so none of these differences is an artefact of mixing annotation sources (Section 4.4).");
s.addNotes("Bar geometry is schematic within each pair: the lengths and the size of the difference are to "
 + "scale, but the position of the internal block in panels B and C is indicative. Panel A's extra copies "
 + "are shown at the C-terminal end for clarity; in the gene they lie in the low-complexity linker.");
}
// ============================================================ PROPOSED VENN FIGURE
{
const circle = (s,cx,cy,d,col) => s.addShape(pptx.ShapeType.ellipse, {
  x:cx-d/2, y:cy-d/2, w:d, h:d,
  fill:{color:col, transparency:80}, line:{color:col, width:1.5} });
const num = (s,cx,cy,t,sz) => s.addText(String(t), {
  x:cx-0.6, y:cy-0.22, w:1.2, h:0.44, fontFace:F, fontSize:sz||20, bold:true,
  color:INK, align:"center", valign:"middle" });
const lab = (s,cx,cy,t,col,sz,w) => s.addText(t, {
  x:cx-(w||2.4)/2, y:cy-0.24, w:w||2.4, h:0.48, fontFace:F, fontSize:sz||11,
  bold:true, color:col, align:"center", valign:"middle" });
const pan = (s,x,y,w,h) => s.addShape(pptx.ShapeType.roundRect, {
  x,y,w,h, rectRadius:0.04, fill:{color:PANEL}, line:{color:EDGE, width:1} });

// ---- A
{ const s = S();
head(s, "Proposed Figure 6 (A)  —  Overlap between the mechanistic categories",
        "1,503 orthologous CDS pairs that differ in length (N1 + N2 pooled). Circles are not area-proportional.");
pan(s, 0.55, 1.12, 8.90, 3.98);
const D = 2.50;
circle(s, 3.85, 2.60, D, BLUE); circle(s, 5.75, 2.60, D, ORANGE); circle(s, 4.80, 3.75, D, RED);
num(s,3.05,2.32,370); num(s,6.55,2.32,373); num(s,4.80,4.32,4);
num(s,4.80,2.22,41); num(s,5.52,3.38,14,16); num(s,4.08,3.38,0,16); num(s,4.80,2.98,0,16);
lab(s,2.55,0.99,"Tandem repeat differs\nin copy number",BLUE,10.5,2.5);
lab(s,7.10,0.99,"Alternative start\nor stop codon",ORANGE,10.5,2.5);
lab(s,4.80,4.72,"IS element",RED,10.5,1.8);
s.addText("701",{x:0.80,y:1.30,w:1.3,h:0.44,fontFace:F,fontSize:22,bold:true,color:INK2,align:"left",valign:"middle"});
s.addText("no mechanism detected",{x:0.80,y:1.72,w:2.6,h:0.30,fontFace:F,fontSize:10,italic:true,color:INK2,align:"left",valign:"middle"});
foot(s,"Circle totals reproduce Table 2: tandem repeat 411 = 207 (N1) + 204 (N2); IS 18 = 10 + 8. Source: cds_offsets_IS.tsv, both neighbours.");
s.addNotes("N1 / N2 / pooled: TR only 184/186/370 | codon only 187/186/373 | IS only 2/2/4 | "
 + "TR+codon 23/18/41 | codon+IS 8/6/14 | TR+IS 0/0/0 | all three 0/0/0 | none 354/347/701.\n"
 + "Caption caveat: membership of the codon-shift circle means a non-zero start or stop offset, not that "
 + "the offset accounts for the whole length difference.");
}
// ---- B
{ const s = S();
head(s, "Proposed Figure 6 (B)  —  Large-indel pool, shared between the two neighbours",
        "56 distinct reference CDS in the 69-pair large-indel pool. Circles are not area-proportional.");
const D = 2.90;
circle(s,4.20,3.05,D,BLUE); circle(s,6.10,3.05,D,ORANGE);
num(s,3.35,2.90,26,24); num(s,5.15,2.90,13,24); num(s,6.95,2.90,17,24);
lab(s,3.00,1.28,"vs N1\nGCA_051994785.1",BLUE,10.5,2.5);
lab(s,7.30,1.28,"vs N2\nGCA_045689975.1",ORANGE,10.5,2.5);
s.addText("10 of the 13 shared genes have an identical length difference on both sides,\nwhich places the event on the reference lineage",
  {x:2.40,y:4.62,w:5.60,h:0.52,fontFace:F,fontSize:10,italic:true,color:INK2,align:"center",valign:"middle"});
foot(s,"Reproduces the counts already stated in Section 3.3. Source: cds_offsets_IS.tsv, rows with IS_class = large_indel?.");
}
// ---- C
{ const s = S();
head(s, "Proposed Figure 6 (C)  —  Orthologue coverage of the reference gene set",
        "3,830 coding sequences called on the reference genome with Prodigal. Circles are not area-proportional.");
pan(s, 0.55, 1.12, 8.90, 3.98);
const D = 2.70;
circle(s,4.35,3.05,D,BLUE); circle(s,6.00,3.05,D,ORANGE);
num(s,3.50,3.05,"200"); num(s,5.17,3.05,"2,941"); num(s,6.85,3.05,"68");
lab(s,3.20,1.42,"ortholog in N1",BLUE,11,2.2);
lab(s,7.15,1.42,"ortholog in N2",ORANGE,11,2.2);
s.addText("621",{x:0.80,y:1.30,w:1.3,h:0.44,fontFace:F,fontSize:22,bold:true,color:INK2,align:"left",valign:"middle"});
s.addText("no ortholog in either neighbour",{x:0.80,y:1.72,w:3.2,h:0.30,fontFace:F,fontSize:10,italic:true,color:INK2,align:"left",valign:"middle"});
foot(s,"2,941 + 200 = 3,141 = the N1 pair count; 2,941 + 68 = 3,009 = the N2 pair count. 621 of 3,830 (16.2%) never enter the comparison.");
}
}

pptx.writeFile({ fileName: OUT }).then(() => console.log("saved:", OUT));
