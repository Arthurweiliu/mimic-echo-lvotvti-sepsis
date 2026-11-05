#!/usr/bin/env node
/* Portable: set H1_WORK (this repository) and MIMIC_ECHO_DIR (your derived tables). */
/** md2docx_h1.mjs — convert the H1 manuscript markdown to a Word file.
 *  Supports: #/##/### headings, **bold**, *italic*, pipe tables, bullet lists,
 *  pipe-free paragraphs. BMC-style A4, 11 pt, 1.5 line spacing. */
import fs from 'node:fs';
import path from 'node:path';
import { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
         HeadingLevel, AlignmentType, WidthType, BorderStyle, PageBreak } from 'docx';

const src = process.argv[2] || file.path(WORK, "MANUSCRIPT_BMC_ANESTHESIOLOGY.md");
const out = process.argv[3] || file.path(WORK, "H1_BMC_Anesthesiology_manuscript.docx");
const md = fs.readFileSync(src, 'utf8').replace(/\r\n/g, '\n');
const lines = md.split('\n');

const runs = (text, opts = {}) => {
  const out = [];
  const parts = text.split(/(\*\*[^*]+\*\*|\*[^*]+\*)/g).filter(s => s !== '');
  for (const p of parts) {
    if (p.startsWith('**') && p.endsWith('**') && p.length > 4)
      out.push(new TextRun({ text: p.slice(2, -2), bold: true, size: opts.size ?? 22, font: 'Times New Roman' }));
    else if (p.startsWith('*') && p.endsWith('*') && p.length > 2)
      out.push(new TextRun({ text: p.slice(1, -1), italics: true, size: opts.size ?? 22, font: 'Times New Roman' }));
    else out.push(new TextRun({ text: p, size: opts.size ?? 22, font: 'Times New Roman', bold: !!opts.bold }));
  }
  return out;
};

const children = [];
let i = 0;
let tableCount = 0;
while (i < lines.length) {
  const line = lines[i];
  if (/^\s*$/.test(line)) { i++; continue; }

  // pipe table
  if (/^\|/.test(line) && /^\|[\s:|-]+\|?\s*$/.test((lines[i + 1] || '').trim())) {
    const header = line.split('|').slice(1, -1).map(s => s.trim());
    const rows = [];
    i += 2;
    while (i < lines.length && /^\|/.test(lines[i])) {
      rows.push(lines[i].split('|').slice(1, -1).map(s => s.trim()));
      i++;
    }
    tableCount++;
    const cellPara = (t, bold, align) => new Paragraph({
      children: [new TextRun({ text: t, bold, size: 16, font: 'Times New Roman' })],
      alignment: align, spacing: { before: 20, after: 20 }
    });
    children.push(new Table({
      width: { size: 100, type: WidthType.PERCENTAGE },
      borders: {
        top: { style: BorderStyle.SINGLE, size: 4, color: '000000' },
        bottom: { style: BorderStyle.SINGLE, size: 4, color: '000000' },
        left: { style: BorderStyle.SINGLE, size: 4, color: '000000' },
        right: { style: BorderStyle.SINGLE, size: 4, color: '000000' },
        insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: '666666' },
        insideVertical: { style: BorderStyle.SINGLE, size: 2, color: '666666' }
      },
      rows: [
        new TableRow({
          tableHeader: true,
          children: header.map(h => new TableCell({
            shading: { fill: 'EFEFEF' },
            children: [cellPara(h, true, AlignmentType.CENTER)]
          }))
        }),
        ...rows.map(r => new TableRow({
          children: r.map((c, ci) => new TableCell({
            children: [cellPara(c, false, ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER)]
          }))
        }))
      ]
    }));
    children.push(new Paragraph({ text: '', spacing: { after: 120 } }));
    continue;
  }

  if (/^# /.test(line)) {
    children.push(new Paragraph({
      children: [new TextRun({ text: line.slice(2), bold: true, size: 28, font: 'Times New Roman' })],
      spacing: { before: 120, after: 200 }, alignment: AlignmentType.LEFT
    }));
    i++; continue;
  }
  if (/^## /.test(line)) {
    children.push(new Paragraph({
      children: [new TextRun({ text: line.slice(3), bold: true, size: 24, font: 'Times New Roman' })],
      spacing: { before: 260, after: 120 },
      pageBreakBefore: /^(References|Tables|Figure legends|Declarations|List of abbreviations)$/.test(line.slice(3).trim())
    }));
    i++; continue;
  }
  if (/^### /.test(line)) {
    children.push(new Paragraph({
      children: [new TextRun({ text: line.slice(4), bold: true, size: 22, font: 'Times New Roman' })],
      spacing: { before: 200, after: 100 }
    }));
    i++; continue;
  }
  if (/^---+\s*$/.test(line)) { i++; continue; }
  if (/^[-*] /.test(line)) {
    children.push(new Paragraph({
      children: runs(line.slice(2)),
      bullet: { level: 0 },
      spacing: { after: 60, line: 300 }
    }));
    i++; continue;
  }
  if (/^\d+\. /.test(line)) {
    children.push(new Paragraph({
      children: runs(line.replace(/^\d+\.\s*/, '')),
      numbering: undefined,
      indent: { left: 360 },
      spacing: { after: 60, line: 300 }
    }));
    i++; continue;
  }
  children.push(new Paragraph({
    children: runs(line),
    spacing: { after: 140, line: 320 },
    alignment: AlignmentType.JUSTIFIED
  }));
  i++;
}

const doc = new Document({
  creator: 'W. Liu',
  title: 'LVOT-VTI and LVEF in sepsis',
  styles: { default: { document: { run: { font: 'Times New Roman', size: 22 } } } },
  sections: [{
    properties: {
      page: {
        size: { width: 11906, height: 16838 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
      }
    },
    children
  }]
});

const buf = await Packer.toBuffer(doc);
fs.writeFileSync(out, buf);
console.log(`wrote ${out} (${children.length} blocks, ${tableCount} tables, ${(buf.length / 1024).toFixed(0)} KB)`);
