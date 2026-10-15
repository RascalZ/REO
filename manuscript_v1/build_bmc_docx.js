const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  WidthType, BorderStyle, AlignmentType, HeadingLevel, Footer,
  PageNumber, LineNumberRestartFormat
} = require('docx');

const root = __dirname;
const input = path.join(root, 'Manuscript_BMC_Genomics_draft.md');
const output = path.join(root, 'Manuscript_BMC_Genomics_draft.docx');
const markdown = fs.readFileSync(input, 'utf8');
const lines = markdown.split(/\r?\n/);
const children = [];
const PAGE_WIDTH = 11906;
const PAGE_HEIGHT = 16838;
const TABLE_WIDTH = PAGE_WIDTH - 2 * 1440;
let tableIndex = 0;
let currentSection = '';

function runsFor(text, options = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0;
  for (const match of text.matchAll(re)) {
    if (match.index > last) out.push(new TextRun({ text: text.slice(last, match.index), ...options }));
    const raw = match[0];
    if (raw.startsWith('**')) out.push(new TextRun({ text: raw.slice(2, -2), bold: true, ...options }));
    else out.push(new TextRun({ text: raw.slice(1, -1), italics: true, ...options }));
    last = match.index + raw.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), ...options }));
  if (out.length === 0) out.push(new TextRun({ text: '', ...options }));
  return out;
}
function p(text, opts = {}) {
  const isPlaceholder = text.includes('[AUTHOR INPUT REQUIRED') || text.includes('[NOT FOR SUBMISSION');
  let runs;
  if (isPlaceholder) {
    runs = [new TextRun({ text: text.replace(/\*\*/g, '').replace(/`/g, ''), color: 'A00000', bold: true, font: 'Times New Roman', size: 22 })];
  } else {
    runs = runsFor(text, { font: 'Times New Roman', size: 24 });
  }
  return new Paragraph({
    children: runs,
    spacing: { line: 480, after: 120, before: opts.before ?? 0 },
    alignment: opts.align,
    heading: opts.heading,
    keepNext: opts.keepNext ?? false,
    indent: opts.indent,
    ...opts.paragraphOptions,
  });
}
function cellPara(text, header = false) {
  return new Paragraph({
    children: runsFor(text, { font: 'Times New Roman', size: 16, bold: header }),
    spacing: { line: 220, before: 20, after: 20 },
  });
}
function cell(text, width, header = false) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    margins: { top: 60, bottom: 60, left: 70, right: 70 },
    children: [cellPara(text, header)],
  });
}
function makeTable(rawRows) {
  tableIndex += 1;
  const widths = tableIndex === 1 ? [1500, 900, 2050, 900, TABLE_WIDTH - 5350] : [1400, 900, 2750, 1150, TABLE_WIDTH - 6200];
  if (widths.reduce((a,b)=>a+b,0) !== TABLE_WIDTH) throw new Error('Table widths mismatch');
  const rows = rawRows.map((row, ri) => new TableRow({
    tableHeader: ri === 0,
    children: row.map((value, ci) => cell(value, widths[ci], ri === 0)),
  }));
  const border = { style: BorderStyle.SINGLE, size: 2, color: '808080' };
  const borders = { top: border, bottom: border, left: border, right: border, insideHorizontal: border, insideVertical: border };
  return new Table({
    width: { size: TABLE_WIDTH, type: WidthType.DXA },
    columnWidths: widths,
    borders,
    rows,
  });
}
function cellsFromLine(line) {
  return line.trim().replace(/^\|/, '').replace(/\|$/, '').split('|').map(x => x.trim());
}

// A visible status line prevents accidental use of this author-placeholder draft as a final submission.
children.push(p('DRAFT — AUTHOR INPUT REQUIRED — NOT FOR SUBMISSION', { align: AlignmentType.CENTER, keepNext: true }));

for (let i = 0; i < lines.length;) {
  const line = lines[i].trimEnd();
  if (!line.trim()) { i++; continue; }
  if (line.startsWith('# BMC Genomics Research Article')) { i++; continue; }
  if (line.startsWith('**Draft status:**')) {
    children.push(p(line.replace(/\*\*/g, ''), { align: AlignmentType.CENTER, keepNext: true }));
    i++; continue;
  }
  if (line.startsWith('|')) {
    const rawRows = [];
    while (i < lines.length && lines[i].trim().startsWith('|')) {
      const row = cellsFromLine(lines[i]);
      if (!row.every(x => /^:?-{2,}:?$/.test(x))) rawRows.push(row);
      i++;
    }
    if (rawRows.length) children.push(makeTable(rawRows));
    continue;
  }
  const h1 = line.match(/^##\s+(.*)$/);
  if (h1) {
    if (h1[1] === 'Title page') { i++; continue; }
    currentSection = h1[1];
    children.push(new Paragraph({
      children: [new TextRun({ text: h1[1], font: 'Times New Roman', size: 26, bold: true, color: '000000' })], heading: HeadingLevel.HEADING_1,
      spacing: { line: 480, before: 260, after: 140 }, keepNext: true,
      run: { font: 'Times New Roman', size: 26, bold: true },
    }));
    i++; continue;
  }
  const h2 = line.match(/^###\s+(.*)$/);
  if (h2) {
    children.push(new Paragraph({
      children: [new TextRun({ text: h2[1], font: 'Times New Roman', size: 24, bold: true, color: '000000' })], heading: HeadingLevel.HEADING_2,
      spacing: { line: 480, before: 180, after: 80 }, keepNext: true,
      run: { font: 'Times New Roman', size: 24, bold: true },
    }));
    i++; continue;
  }
  const text = line.replace(/\\([*_])/g, '$1').replace(/`/g, '');
  let paragraph = p(text, {
    keepNext: /^\*\*(?:Table|Figure) \d+\./.test(text),
    paragraphOptions: /^\d+\. /.test(text) ? { indent: { left: 360, hanging: 360 } } : {},
  });
  // Additional-file descriptions and figure legends are compact but remain
  // close to double-spaced manuscript text; this avoids orphaning the final
  // legend on a nearly empty page in the rendered Word copy.
  if (currentSection === 'Additional files' || currentSection === 'Figure legends') {
    paragraph = new Paragraph({
      children: runsFor(text, { font: 'Times New Roman', size: 20 }),
      spacing: { line: currentSection === 'Figure legends' ? 340 : 360, after: 80 },
      keepNext: /^\*\*Additional file \d+\./.test(text) || (currentSection === 'Figure legends' && /^Figure [1-3]\./.test(text)),
      keepLines: currentSection === 'Figure legends',
    });
  }
  // The title is the first non-heading line beginning with the explicit title label.
  if (text.startsWith('**Title:**')) paragraph = new Paragraph({ children: runsFor(text, {font:'Times New Roman',size:28}), spacing:{line:480,after:180}, keepNext:true });
  if (currentSection === 'Figure legends' && /^Figure \d+\./.test(text)) {
    paragraph = new Paragraph({
      children: runsFor(text, { font: 'Times New Roman', size: 20 }),
      spacing: { line: 340, after: 80, before: 0 },
      keepLines: true,
      keepNext: /^Figure [1-3]\./.test(text),
    });
  } else if (/^Figure \d+\./.test(text) || /^Table \d+\./.test(text)) {
    paragraph = new Paragraph({ children: runsFor(text,{font:'Times New Roman',size:20}), spacing:{line:360,after:160,before:100}, keepNext:/^Table \d+\./.test(text) });
  }
  if (currentSection === 'Conclusions') {
    paragraph = new Paragraph({
      children: runsFor(text, { font: 'Times New Roman', size: 24 }),
      spacing: { line: 480, after: 120 },
      keepLines: true,
    });
  }
  if (text.startsWith('*Note:*')) paragraph = new Paragraph({ children:runsFor(text,{font:'Times New Roman',size:18}), spacing:{line:300,after:220} });
  children.push(paragraph);
  i++;
}

const footer = new Footer({ children: [new Paragraph({
  alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: 'Page ', font: 'Times New Roman', size: 18 }), PageNumber.CURRENT],
})] });
const doc = new Document({
  creator: 'Project manuscript build',
  title: 'BMC Genomics Research Article draft — author input required',
  description: 'Evidence-bounded BMC Genomics manuscript draft; not submission-ready.',
  styles: {
    default: { document: { run: { font: 'Times New Roman', size: 24 }, paragraph: { spacing: { line: 480, after: 120 } } } },
  },
  sections: [{
    properties: {
      page: { size: { width: PAGE_WIDTH, height: PAGE_HEIGHT }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } },
      lineNumbers: { countBy: 1, start: 1, restart: LineNumberRestartFormat.NEW_PAGE, distance: 360 },
    },
    footers: { default: footer },
    children,
  }],
});
Packer.toBuffer(doc).then(buf => { fs.writeFileSync(output, buf); console.log(output); }).catch(err => { console.error(err); process.exit(1); });
