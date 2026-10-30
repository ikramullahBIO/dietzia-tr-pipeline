const pptxgen = require('pptxgenjs');

// geometry copied from dietzia_rerun/08_make_figure2.py
const BLUE="2A78D6", ORANGE="EB6834", GREY="B8B8B4", RED="C0392B";
const INK="0B0B0B", INK2="52514E", TINT="F2F6FC";
const BLOCKS=[
 {n:"1", title:"Tandem repeat copy-number change",
  rule:"the length difference sits on a repeat array whose copy number differs",
  ref:[[0,26,'cds'],[26,9,'rep'],[35,9,'rep'],[44,9,'rep'],[53,26,'cds']],
  ngb:[[0,26,'cds'],[26,9,'rep'],[35,26,'cds']],
  note:"3 copies vs 1  →  shorter by 2 × the repeat unit",
  stat:"204 pairs against N1 and 202 against N2  —  26.9–27.1% of the pairs that differ in length, and 11.3–12.3% of all base-pair variation",
  notes:"The category the project is actually about. A repeat array inside the gene carries a different number of copies in the two genomes, so the length difference is a multiple of the repeat unit. The two independently compared neighbours agree closely, which is hard to attribute to chance. The stricter position-and-magnitude re-test confirms 321 of 406 of these pairs."},
 {n:"2", title:"Alternative start codon",
  rule:"the difference is entirely at the 5′ end; the DNA is present in both genomes",
  ref:[[0,14,'out'],[14,60,'cds']],
  ngb:[[14,60,'cds']],
  note:"start codon moved 3′; the shared body is unchanged",
  stat:"71 pairs against N1 and 68 against N2  —  8.5–10.5% of all base-pair variation",
  notes:"Not a difference in the DNA at all. Both genomes carry the same sequence and the gene caller simply picked a different start codon. Figure 9A shows this at nucleotide level for a real pair, where both genomes carry the same GTG at codon 1 and the same ATG at codon 26."},
 {n:"3", title:"Alternative stop codon",
  rule:"the difference is entirely at the 3′ end",
  ref:[[0,60,'cds'],[60,16,'out']],
  ngb:[[0,60,'cds']],
  note:"stop codon moved 5′; the shared body is unchanged",
  stat:"23 pairs against N1 and 26 against N2  —  0.6–0.9% of all base-pair variation",
  notes:"The mirror image of category 2, at the other end of the gene. Same logic: the shared body of the protein is unchanged, only the terminal codon call moved."},
 {n:"4", title:"Insertion-sequence element",
  rule:"a large block matching a transposase or integrase interrupts the gene",
  ref:[[0,74,'cds']],
  ngb:[[0,22,'cds'],[22,30,'is'],[52,22,'cds']],
  note:"large foreign block inserted",
  stat:"only 10 pairs against N1 and 8 against N2  —  but 7.3–8.0% of all base-pair variation",
  notes:"Rare by count but large by base pairs. Detected from the product annotation alone, which is the main limitation of this category: a curated IS database such as ISfinder would make the calls firmer."},
 {n:"5", title:"Internal duplication, otherwise an unresolved internal indel",
  rule:"an internal block ≥10 aa matching its own adjacent flank at ≥70% identity",
  ref:[[0,20,'cds'],[20,15,'rep'],[35,39,'cds']],
  ngb:[[0,20,'cds'],[20,15,'rep'],[35,15,'dup'],[50,39,'cds']],
  note:"block duplicated immediately next to its source",
  stat:"a tested negative: 0 of 353 pairs reached 70% flank identity, and the count stays at zero for every threshold at and above 55%",
  notes:"The last test in the list, and it returned a clean negative. Rather than test one cut-off, the screen was repeated from 40 to 95 percent; the highest flank identity seen anywhere is 50 percent. Anything reaching this point without matching is recorded as an unresolved internal indel, which is the largest residual category and the study's main open question."},
];

function style(k){
  switch(k){
    case 'cds': return {fill:{color:BLUE},   line:{color:"FFFFFF", width:1}};
    case 'rep': return {fill:{color:ORANGE}, line:{color:"FFFFFF", width:1}};
    case 'dup': return {fill:{color:ORANGE}, line:{color:"FFFFFF", width:1.5, dashType:"dash"}};
    case 'out': return {fill:{color:GREY},   line:{color:"6E6E6A", width:1,   dashType:"dash"}};
    case 'is':  return {fill:{color:RED},    line:{color:"FFFFFF", width:1.5, dashType:"dash"}};
  }
}
const pres = new pptxgen();
pres.layout='LAYOUT_WIDE';
pres.author='Ikram Ullah';
pres.title='Figure 2 - CDS length-difference classification';

function drawBlock(slide,b,x,y,w,o){
  const bx=x+o.labW;
  slide.addText([{text:b.n+"   ",options:{color:ORANGE,bold:true}},{text:b.title,options:{color:INK,bold:true}}],
    {x, y, w, h:o.titleH, fontSize:o.titleSz, fontFace:'Calibri', isTextBox:true, margin:0, valign:'top'});
  slide.addText(b.rule,{x, y:y+o.titleH, w, h:o.ruleH, fontSize:o.ruleSz, color:INK2,
    fontFace:'Calibri', isTextBox:true, margin:0, valign:'top'});
  const yRef=y+o.titleH+o.ruleH+o.preBars, yNgb=yRef+o.barH+o.gap;
  for (const [lab,yy] of [['reference',yRef],['neighbour',yNgb]])
    slide.addText(lab,{x, y:yy-0.03, w:o.labW-0.14, h:o.barH+0.06, fontSize:o.labSz, color:INK2,
      fontFace:'Calibri', align:'right', valign:'middle', isTextBox:true, margin:0});
  for (const [x0,wd,k] of b.ref)
    slide.addShape(pres.ShapeType.rect, Object.assign({x:bx+x0*o.scale, y:yRef, w:wd*o.scale, h:o.barH}, style(k)));
  for (const [x0,wd,k] of b.ngb)
    slide.addShape(pres.ShapeType.rect, Object.assign({x:bx+x0*o.scale, y:yNgb, w:wd*o.scale, h:o.barH}, style(k)));
  const yNote=yNgb+o.barH+o.postBars;
  slide.addText(b.note,{x:bx, y:yNote, w:w-o.labW, h:o.noteH, fontSize:o.noteSz, color:INK2,
    italic:true, fontFace:'Calibri', isTextBox:true, margin:0, valign:'top'});
  return (yNote+o.noteH)-y;
}
function drawLegend(slide,x,y,w,sz,twoCol,short){
  const items=[['cds','coding sequence, shared'],['rep','tandem-repeat unit / duplicated block'],
               ['out', short ? 'present in the DNA, outside the CDS' : 'present in the DNA, outside the called CDS'],['is','IS-element / transposase block']];
  const colW=twoCol?w/2:w, rowH=sz>11?0.34:0.30;
  items.forEach((it,i)=>{
    const cx=x+(twoCol?(i%2)*colW:0), cy=y+(twoCol?Math.floor(i/2):i)*rowH;
    slide.addShape(pres.ShapeType.rect, Object.assign({x:cx, y:cy+0.05, w:0.34, h:0.16}, style(it[0])));
    slide.addText(it[1],{x:cx+0.46, y:cy, w:colW-0.55, h:0.26, fontSize:sz, color:INK,
      fontFace:'Calibri', isTextBox:true, margin:0, valign:'middle'});
  });
}

// ---------- slide 1: whole figure ----------
const s1=pres.addSlide(); s1.background={color:'FFFFFF'};
s1.addText('How a coding-sequence length difference is assigned to a category',
  {x:0.5,y:0.32,w:12.4,h:0.46,fontSize:25,bold:true,color:INK,fontFace:'Calibri',isTextBox:true,margin:0});
s1.addText('Each pair of bars is one orthologous CDS, drawn to scale. Categories are tested in the order shown; the first match wins.',
  {x:0.5,y:0.80,w:12.4,h:0.30,fontSize:12.5,color:INK2,fontFace:'Calibri',isTextBox:true,margin:0});
const O1={scale:0.056,barH:0.26,gap:0.14,labW:0.95,titleH:0.30,ruleH:0.26,preBars:0.10,postBars:0.12,
          noteH:0.26,titleSz:13.5,ruleSz:10.5,labSz:10,noteSz:10};
let y=1.32; for(const b of BLOCKS.slice(0,3)) y+=drawBlock(s1,b,0.5,y,6.05,O1)+0.34;
let y2=1.32; for(const b of BLOCKS.slice(3)) y2+=drawBlock(s1,b,6.85,y2,6.05,O1)+0.34;
drawLegend(s1,6.85,6.42,6.05,10.5,true,true);
s1.addNotes('The classification is a decision list, not five independent tests: each orthologous CDS pair is checked against the categories in order and assigned to the first that matches, so every pair is counted exactly once. Three quarters of pairs never reach the list at all, because the two genes are identical in length.');

// ---------- slides 2-6: one category each ----------
const O2={scale:0.118,barH:0.55,gap:0.34,labW:1.65,titleH:0.62,ruleH:0.40,preBars:0.22,postBars:0.20,
          noteH:0.34,titleSz:26,ruleSz:15,labSz:13,noteSz:14};
for(const b of BLOCKS){
  const s=pres.addSlide(); s.background={color:'FFFFFF'};
  s.addText('Category '+b.n+' of 5',{x:0.7,y:0.40,w:6,h:0.28,fontSize:13,bold:true,color:ORANGE,
    fontFace:'Calibri',isTextBox:true,margin:0});
  drawBlock(s,b,0.7,0.90,11.9,O2);
  s.addShape(pres.ShapeType.rect,{x:0.7,y:4.55,w:11.9,h:0.92,fill:{color:TINT},line:{color:TINT}});
  s.addText(b.stat,{x:1.0,y:4.55,w:11.3,h:0.92,fontSize:15,color:INK,fontFace:'Calibri',
    isTextBox:true,margin:0,valign:'middle'});
  drawLegend(s,0.7,6.05,11.9,12,true);
  s.addNotes(b.notes);
}
pres.writeFile({fileName:'Figure2_CDS_classification.pptx'}).then(f=>console.log('wrote',f));
